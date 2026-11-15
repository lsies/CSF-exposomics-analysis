from pathlib import Path
import csv
import hashlib
import json
import math
import zipfile
import numpy as np
import pandas as pd
import networkx as nx

def read_enrichment(root, prefix, family):
    with zipfile.ZipFile(root/'raw/public'/f'{prefix}all.zip') as archive:
        enrichment=pd.read_csv(archive.open('Enrichment_GO/GO_AllLists.csv'))
    workbook=pd.read_excel(root/'raw/public'/f'{prefix}metascape_result.xlsx',sheet_name='Enrichment')
    categories={'GO_BP':'GO Biological Processes','GO_CC':'GO Cellular Components','GO_MF':'GO Molecular Functions','KEGG':'KEGG Pathway'}
    if family in categories:
        assert set(enrichment.Category)=={categories[family]},f'Wrong category in archived CSV: {family}'
        assert set(workbook.Category.dropna())=={categories[family]},f'Wrong category in archived XLSX: {family}'
    exact=workbook[workbook.GroupID.str.endswith('_Member')].drop_duplicates(['Category','Term'])
    exact=exact.rename(columns={'Term':'GO','LogP':'precise_LogP','Log(q-value)':'precise_LogQ'})
    enrichment=enrichment.merge(exact[['Category','GO','precise_LogP','precise_LogQ']],on=['Category','GO'],how='left',validate='one_to_one')
    enrichment['numeric_precision']=np.where(enrichment.precise_LogP.notna(),'full_precision_XLSX','rounded_original_CSV')
    enrichment['log10_p']=enrichment.precise_LogP.fillna(enrichment.LogP)
    enrichment['log10_q']=enrichment.precise_LogQ.fillna(enrichment['Log(q-value)'])
    enrichment['p_value']=10.0**enrichment.log10_p
    enrichment['q_value']=10.0**enrichment.log10_q
    enrichment['significant_FDR05']=enrichment.q_value.lt(.05)
    enrichment['threshold_status']=np.where(enrichment.significant_FDR05,'FDR_below_0.05','FDR_not_below_0.05')
    boundary=enrichment.numeric_precision.eq('rounded_original_CSV') & enrichment.log10_q.between(-1.35,-1.25)
    enrichment.loc[boundary,'threshold_status']='unresolved_at_exported_precision'
    enrichment.loc[boundary,'significant_FDR05']=False
    enrichment['analysis_family']=family
    return enrichment

def main(root):
    destination=root/'results'
    frozen=json.loads((root/'data/frozen_analysis_inputs.json').read_text(encoding='utf-8'))
    for relative,digest in frozen.items():
        assert hashlib.sha256((root/relative).read_bytes()).hexdigest()==digest,f'External analysis must be rerun for changed input: {relative}'
    table=pd.read_csv(root/'raw/cytoscape/node_table.csv')
    with (root/'raw/cytoscape/cytoHubba_scores.csv').open(encoding='utf-8-sig',newline='') as handle:
        rows=list(csv.reader(handle))
    hub=pd.DataFrame([row[:len(rows[0])] for row in rows[1:]],columns=rows[0]).rename(columns={'node_name':'gene_symbol'})
    for column in hub.columns[1:]:
        hub[column]=pd.to_numeric(hub[column])
    hub.to_csv(destination/'cytoHubba_all_scores.csv',index=False)
    nodes=table[['gene_symbol','Degree','Betweenness','Closeness','SUID']].merge(hub[['gene_symbol','MCC']],on='gene_symbol',validate='one_to_one')
    medians=nodes[['Degree','Betweenness','Closeness']].median()
    nodes['CytoNCA_above_three_medians']=nodes[medians.index].gt(medians).all(axis=1)
    nodes['MCC_rank']=nodes.MCC.rank(method='min',ascending=False).astype(int)
    nodes['MCC_top10_including_ties']=nodes.MCC_rank.le(10)
    clusters=json.loads((root/'raw/cytoscape/MCODE_results.json').read_text())['data']['clusters']
    names=dict(zip(nodes.SUID,nodes.gene_symbol))
    memberships=[]
    summaries=[]
    for cluster in clusters:
        genes=sorted(names[x] for x in cluster['nodes'])
        summaries.append({'module':cluster['rank'],'score':cluster['score'],'size':len(genes),'seed':names[cluster['seedNode']],'genes':';'.join(genes)})
        memberships.extend({'module':cluster['rank'],'gene_symbol':gene} for gene in genes)
    modules=pd.DataFrame(memberships)
    modules.to_csv(destination/'MCODE_memberships.csv',index=False)
    pd.DataFrame(summaries).to_csv(destination/'MCODE_modules.csv',index=False)
    nodes['MCODE_module']=nodes.gene_symbol.map(modules.groupby('gene_symbol').module.apply(lambda x:';'.join(map(str,sorted(x))))).fillna('')
    nodes['consensus_hub']=nodes.CytoNCA_above_three_medians&nodes.MCC_top10_including_ties&nodes.MCODE_module.ne('')
    kp=pd.read_csv(root/'data/kp_gene_reference.csv')
    nodes['KP_core']=nodes.gene_symbol.isin(kp.loc[kp.is_kp_core,'gene_symbol'])
    response=pd.read_csv(destination/'metabolite_gene_expression_evidence.csv')
    nodes['exposure_responsive']=nodes.gene_symbol.isin(response.loc[response.responsive,'gene_symbol'])
    nodes.sort_values(['MCC','gene_symbol'],ascending=[False,True]).to_csv(destination/'network_node_statistics.csv',index=False)
    edges=pd.read_csv(destination/'STRING_PPI_edges.csv')
    graph=nx.from_pandas_edgelist(edges,'gene_A','gene_B',edge_attr='score')
    degree=dict(graph.degree())
    close=nx.closeness_centrality(graph)
    between=nx.betweenness_centrality(graph,normalized=False)
    mcc={gene:0 for gene in graph}
    for clique in nx.find_cliques(graph):
        weight=math.factorial(len(clique)-1)
        for gene in clique:
            mcc[gene]+=weight
    assert all(nodes.set_index('gene_symbol').Degree.eq(pd.Series(degree)))
    assert np.allclose(nodes.Closeness,[close[x] for x in nodes.gene_symbol],rtol=1e-6)
    assert np.allclose(nodes.Betweenness,[2*between[x] for x in nodes.gene_symbol],rtol=1e-6)
    assert np.allclose(nodes.MCC,[mcc[x] for x in nodes.gene_symbol],rtol=1e-5)
    nx.set_node_attributes(graph,nodes.set_index('gene_symbol').drop(columns='SUID').to_dict('index'))
    nx.write_graphml(graph,destination/'STRING_annotated_network.graphml')
    proteins=pd.read_csv(root/'data/protein_background.csv')
    mapping=pd.read_csv(destination/'STRING_identifier_mapping.csv')
    coverage=proteins[proteins.selected].copy()
    coverage['STRING_mapped']=coverage.gene_symbol.isin(mapping.queryItem)
    coverage['connected_network']=coverage.gene_symbol.isin(graph)
    coverage.to_csv(destination/'STRING_gene_coverage.csv',index=False)
    original=read_enrichment(root,'Metascape_','original_combined_category_run')
    original.to_csv(destination/'Metascape_original_combined_categories.csv',index=False)
    enrichment=pd.concat([read_enrichment(root,f'Metascape_{category}_',category) for category in ['GO_BP','GO_CC','GO_MF','KEGG']],ignore_index=True)
    columns=['Category','GO','Description','p_value','q_value','significant_FDR05','numeric_precision','threshold_status']
    comparison=original[columns].merge(enrichment[columns],on=['Category','GO','Description'],how='outer',suffixes=('_original','_separate'),validate='one_to_one')
    comparison.to_csv(destination/'FDR_category_comparison.csv',index=False)
    comparison[comparison.Description.str.contains('tryptophan|kynuren',case=False)].to_csv(destination/'KP_FDR_comparison.csv',index=False)
    enrichment.to_csv(destination/'Metascape_enrichment_all_returned_terms.csv',index=False)
    enrichment[enrichment.significant_FDR05].to_csv(destination/'Metascape_enrichment_FDR05.csv',index=False)
    terms=enrichment[enrichment.Category.eq('KEGG Pathway')&enrichment.significant_FDR05].sort_values(['q_value','GO']).head(3)
    pathway_edges=[{'gene_symbol':gene,'pathway_id':row.GO,'pathway_name':row.Description,'q_value':row.q_value} for row in terms.itertuples() for gene in row.Hits.split('|')]
    pd.DataFrame(pathway_edges).to_csv(destination/'KEGG_display_gene_pathway_edges.csv',index=False)
    tryptophan=enrichment[enrichment.GO.eq('hsa00380')].iloc[0]
    assert 'KYNU' in tryptophan.Hits.split('|')
    pd.DataFrame([{'gene_symbol':gene,'pathway_id':'hsa00380','pathway_name':'Tryptophan metabolism','edge_type':'KEGG_annotation','pathway_q_value':tryptophan.q_value} for gene in tryptophan.Hits.split('|')]).to_csv(destination/'KP_KEGG_annotation_edges.csv',index=False)
    primary=response[response.network_primary].copy()
    selected_genes=['KYNU','SLC7A5','SLC7A8','ABCG2','ALAS1','GPAT3','MBOAT1']
    sets=pd.read_csv(root/'data/step5_eligible_analysis_sets.csv')
    primary_sets=sets.loc[sets.network_primary,'dataset_id'].tolist()
    rows=[]
    for gene in selected_genes:
        for dataset in primary_sets:
            match=primary[(primary.gene_symbol==gene)&(primary.dataset_id==dataset)]
            record={'gene_symbol':gene,'dataset_id':dataset,'status':'not_testable','LRT_FDR':np.nan,'max_abs_log2_fold_change':np.nan}
            if len(match):
                row=match.iloc[0]
                testable=pd.notna(row.LRT_FDR) and row.LRT_full_converged and row.LRT_reduced_converged
                record.update(status='responsive' if row.responsive else 'not_responsive' if testable else 'not_testable',LRT_FDR=row.LRT_FDR,max_abs_log2_fold_change=row.max_abs_log2_fold_change)
            rows.append(record)
    pd.DataFrame(rows).to_csv(destination/'KP_expression_evidence_matrix.csv',index=False)
    mapped=pd.read_csv(destination/'differential_metabolite_reaction_protein_edges.csv')
    mapped[mapped.gene_symbol.isin(selected_genes)].drop_duplicates(['FeaturesPeak','gene_symbol']).to_csv(destination/'KP_expression_bridge_mapping.csv',index=False)
    summary={'connected_nodes':len(graph),'STRING_edges':graph.number_of_edges(),'MCODE_modules':len(clusters),'CytoNCA_medians':medians.to_dict(),'CytoNCA_selected':int(nodes.CytoNCA_above_three_medians.sum()),'MCC_top10_including_ties':int(nodes.MCC_top10_including_ties.sum()),'consensus_hubs':sorted(nodes.loc[nodes.consensus_hub,'gene_symbol']),'enrichment_FDR05_by_category':enrichment[enrichment.significant_FDR05].Category.value_counts().to_dict(),'independent_topology_verification':'passed','CytoNCA_betweenness_convention':'ordered pairs; twice NetworkX unnormalized undirected betweenness'}
    (destination/'final_network_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
