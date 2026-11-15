from pathlib import Path
import json
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def main(root):
    proteins=pd.read_csv(root/'data/protein_background.csv')
    universe=set(proteins.gene_symbol)
    selected=set(proteins.loc[proteins.selected,'gene_symbol'])
    response=pd.read_csv(root/'data/step5_all_gene_response.csv.gz')
    metadata=pd.read_csv(root/'data/step5_eligible_analysis_sets.csv')
    metadata=metadata[['dataset_id','network_primary']]
    records=[]
    for dataset,group in response.groupby('dataset_id'):
        tested=set(group.loc[group.LRT_FDR.notna()&group.LRT_full_converged&group.LRT_reduced_converged,'gene_symbol'])
        background=tested&universe
        query=selected&background
        responsive=set(group.loc[group.responsive,'gene_symbol'])&background
        intersection=query&responsive
        contingency=[[len(intersection),len(query-responsive)],[len(responsive-query),len(background-query-responsive)]]
        result=fisher_exact(contingency,alternative='greater')
        records.append({'dataset_id':dataset,'cell_line':group.cell_line.iloc[0],'structure_id':group.structure_id.iloc[0],'background_genes':len(background),'metabolite_genes_tested':len(query),'responsive_background_genes':len(responsive),'overlap_genes':len(intersection),'gene_symbols':';'.join(sorted(intersection)),'odds_ratio':result.statistic,'p_value':result.pvalue})
    stats=pd.DataFrame(records)
    stats['FDR']=multipletests(stats.p_value,method='fdr_bh')[1]
    stats=stats.merge(metadata,on='dataset_id',validate='one_to_one')
    stats.to_csv(root/'results/step5_expression_overlap_tests.csv',index=False)
    selected_response=response[response.gene_symbol.isin(selected)].merge(metadata,on='dataset_id')
    selected_response.to_csv(root/'results/metabolite_gene_expression_evidence.csv',index=False)
    edges=pd.read_csv(root/'results/differential_metabolite_reaction_protein_edges.csv')
    bridges=edges.merge(selected_response[selected_response.responsive],on='gene_symbol')
    bridges.to_csv(root/'results/metabolite_expression_bridge_edges.csv',index=False)
    kp=pd.read_csv(root/'data/kp_gene_reference.csv')
    kp_genes=set(kp.loc[kp.is_kp_core,'gene_symbol'])
    edges[edges.gene_symbol.isin(kp_genes)].to_csv(root/'results/KP_core_metabolite_edges.csv',index=False)
    summary={'responsive_overlap_all_experiments':sorted(set(bridges.gene_symbol)),'responsive_overlap_primary':sorted(set(bridges.loc[bridges.network_primary,'gene_symbol'])),'experiments_FDR05':stats.loc[stats.FDR.lt(.05),'dataset_id'].tolist(),'mapped_KP_core':sorted(selected&kp_genes)}
    (root/'results/expression_integration_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
