from pathlib import Path
import json
import re
import pandas as pd
import networkx as nx

def kegg_records(root):
    records={}
    for path in sorted((root/'raw/public').glob('KEGG_compounds_*.txt')):
        for block in path.read_text(encoding='utf-8').split('///'):
            fields={}
            current=None
            for line in block.splitlines():
                key=line[:12].strip()
                value=line[12:].strip()
                if key:
                    current=key
                    fields.setdefault(key,[])
                if current:
                    fields[current].append(value)
            if 'ENTRY' in fields:
                identifier=fields['ENTRY'][0].split()[0]
                records[identifier]={key:' '.join(values) for key,values in fields.items()}
    return records

def main(root):
    records=kegg_records(root)
    rhea=pd.read_csv(root/'raw/public/Rhea_reactions.tsv',sep='\t').fillna('')
    directions=pd.read_csv(root/'raw/public/Rhea_directions.tsv',sep='\t').astype(str)
    master={}
    for values in directions.itertuples(index=False,name=None):
        for value in values:
            master[re.sub(r'\D','',value)]=re.sub(r'\D','',values[0])
    uniprot=pd.read_csv(root/'raw/public/UniProt_human_Rhea.tsv',sep='\t').fillna('')
    proteins={}
    protein_rows=[]
    for row in uniprot.itertuples(index=False,name=None):
        accession,gene,name,reactions=row
        if not gene:
            continue
        for reaction in re.findall(r'RHEA:(\d+)',reactions):
            reaction=master.get(reaction,reaction)
            for symbol in gene.split(';'):
                symbol=symbol.strip()
                if symbol:
                    proteins.setdefault(reaction,set()).add((symbol,accession,name))
                    protein_rows.append((reaction,symbol,accession,name))
    pd.DataFrame(sorted(set(protein_rows)),columns=['Rhea','gene_symbol','UniProt','protein_name']).to_csv(root/'data/human_reaction_proteins.csv',index=False)
    reaction_kegg={}
    reaction_chebi={}
    reaction_records={}
    for row in rhea.itertuples(index=False,name=None):
        reaction,equation,chebi,kegg,pubmed=row
        identifier=re.sub(r'\D','',reaction)
        reaction_records[identifier]=(equation,pubmed)
        if identifier not in proteins:
            continue
        for k in re.findall(r'R\d{5}',kegg):
            reaction_kegg.setdefault(k,set()).add(identifier)
        for c in re.findall(r'CHEBI:(\d+)',chebi):
            reaction_chebi.setdefault(c,set()).add(identifier)
    microspecies={}
    for item in json.loads((root/'raw/public/ChEBI_microspecies.json').read_text(encoding='utf-8'))['results']['bindings']:
        microspecies.setdefault(item['chebi']['value'].split('_')[-1],set()).add(item['major']['value'].split('_')[-1])
    candidates=pd.read_csv(root/'data/compound_candidates.csv',dtype={'ChEBI':str}).fillna('')
    compound_edges=[]
    compounds=[]
    for identifier,group in candidates.groupby('CompoundID'):
        kegg=group.KEGG.iloc[0]
        record=records.get(kegg,{})
        links={}
        for kr in re.findall(r'R\d{5}',record.get('REACTION','')):
            for rr in reaction_kegg.get(kr,set()):
                links.setdefault(rr,set()).add('KEGG_reaction:'+kr)
        chebi_match=re.search(r'ChEBI:\s*([0-9 ]+)',record.get('DBLINKS',''))
        chebi_ids=set(re.findall(r'\d+',chebi_match.group(1)) if chebi_match else [])
        chebi_ids.update(c for c in group.ChEBI if c)
        for chebi in chebi_ids:
            for rr in reaction_chebi.get(chebi,set()):
                links.setdefault(rr,set()).add('exact_ChEBI:'+chebi)
            for major in microspecies.get(chebi,set()):
                for rr in reaction_chebi.get(major,set()):
                    links.setdefault(rr,set()).add('ChEBI_microspecies:'+chebi+'>'+major)
        compounds.append({'CompoundID':identifier,'database_name':record.get('NAME',group.reference_name.iloc[0]),'database_formula':record.get('FORMULA',group.reference_formula.iloc[0]),'mapped_ChEBI':';'.join(sorted(chebi_ids)),'human_reactions':len(links)})
        for reaction,sources in links.items():
            equation,pubmed=reaction_records[reaction]
            for gene,accession,name in proteins[reaction]:
                compound_edges.append({'CompoundID':identifier,'Rhea':reaction,'gene_symbol':gene,'UniProt':accession,'protein_name':name,'mapping_route':';'.join(sorted(sources)),'equation':equation,'PubMed':pubmed})
    compounds=pd.DataFrame(compounds)
    compounds.to_csv(root/'data/compound_reference.csv',index=False)
    compound_edges=pd.DataFrame(compound_edges).drop_duplicates().sort_values(['CompoundID','Rhea','gene_symbol','UniProt']).reset_index(drop=True)
    candidate_checked=candidates.merge(compounds,on='CompoundID',how='left')
    candidate_checked['formula_match']=candidate_checked.formula.eq(candidate_checked.database_formula)
    candidate_checked.to_csv(root/'results/compound_identity_mapping.csv',index=False)
    valid=candidate_checked[candidate_checked.formula_match]
    edges=valid.merge(compound_edges,on='CompoundID',how='inner').sort_values(['FeaturesPeak','CompoundID','Rhea','gene_symbol','UniProt']).reset_index(drop=True)
    edges.to_csv(root/'results/metabolite_reaction_protein_edges.csv',index=False)
    annotations=pd.read_csv(root/'data/metabolite_annotations.csv')
    mapped=set(edges.FeaturesPeak)
    has_kegg=set(candidates.FeaturesPeak)
    valid_id=set(valid.FeaturesPeak)
    annotations['mapping_status']=annotations.FeaturesPeak.map(lambda x:'human_reaction_mapped' if x in mapped else 'no_curated_human_reaction' if x in valid_id else 'formula_conflict' if x in has_kegg else 'no_resolved_public_identifier')
    annotations['mapped_human_genes']=annotations.FeaturesPeak.map(edges.groupby('FeaturesPeak').gene_symbol.apply(lambda x:';'.join(sorted(set(x))))).fillna('')
    annotations.to_csv(root/'results/all_metabolite_mapping_inventory.csv',index=False)
    selected=edges[edges.selected]
    selected.to_csv(root/'results/differential_metabolite_reaction_protein_edges.csv',index=False)
    all_genes=sorted(set(edges.gene_symbol))
    selected_genes=sorted(set(selected.gene_symbol))
    pd.DataFrame({'gene_symbol':all_genes,'selected':pd.Series(all_genes).isin(selected_genes)}).to_csv(root/'data/protein_background.csv',index=False)
    (root/'data/Metascape_input_genes.txt').write_text('\n'.join(selected_genes)+'\n',encoding='utf-8')
    (root/'data/Metascape_background_genes.txt').write_text('\n'.join(all_genes)+'\n',encoding='utf-8')
    pd.DataFrame({'Input':pd.Series(selected_genes),'Background':pd.Series(all_genes)}).to_csv(root/'data/Metascape_input_background.csv',index=False)
    kp=pd.read_csv(root/'data/kp_gene_reference.csv')
    kp_genes=set(kp.loc[kp.is_kp_core,'gene_symbol'])
    graph=nx.Graph()
    for row in selected.itertuples():
        feature='F:'+row.physical_feature
        compound='C:'+row.CompoundID
        reaction='R:'+row.Rhea
        gene='G:'+row.gene_symbol
        graph.add_node(feature,node_type='LCMS_feature',label=row.physical_feature)
        graph.add_node(compound,node_type='compound',label=row.database_name.split(';')[0],KEGG=row.KEGG)
        graph.add_node(reaction,node_type='reaction',label='RHEA:'+row.Rhea)
        graph.add_node(gene,node_type='human_protein',label=row.gene_symbol,KP_core=row.gene_symbol in kp_genes)
        graph.add_edge(feature,compound,evidence='supplier_candidate_annotation')
        graph.add_edge(compound,reaction,evidence=row.mapping_route)
        graph.add_edge(reaction,gene,evidence='reviewed_UniProt_Rhea_annotation')
    nx.write_graphml(graph,root/'results/metabolite_reaction_protein_network.graphml')
    summary={'selected_annotations':int(annotations.selected_any_comparison.sum()),'selected_physical_features':annotations.loc[annotations.selected_any_comparison,'physical_feature'].nunique(),'selected_mapped_annotations':selected.FeaturesPeak.nunique(),'selected_mapped_features':selected.physical_feature.nunique(),'selected_mapped_compounds':selected.CompoundID.nunique(),'selected_human_reactions':selected.Rhea.nunique(),'selected_human_genes':len(selected_genes),'background_human_genes':len(all_genes),'KP_mapped_genes':sorted(kp_genes&set(selected_genes)),'network_nodes':len(graph),'network_edges':graph.number_of_edges()}
    (root/'results/network_mapping_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
