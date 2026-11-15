import argparse
import json
from pathlib import Path
import urllib.parse
import pandas as pd
from public_sources import fetch
from compound_resolution import resolve

def main(root,online):
    fetch(root,'KEGG_info_compound.txt','https://rest.kegg.jp/info/compound',online=online)
    identifiers=pd.read_csv(root/'data/KEGG_query_ids.csv').KEGG.tolist()
    for i in range(0,len(identifiers),10):
        chunk=identifiers[i:i+10]
        fetch(root,f'KEGG_compounds_{i//10:03d}.txt','https://rest.kegg.jp/get/' + '+'.join(chunk),online=online)
        print('KEGG compounds',min(i+10,len(identifiers)),flush=True)
    query='(organism_id:9606) AND (reviewed:true) AND (fragment:false) AND (cc_catalytic_activity:"rhea:*")'
    url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode({'query':query,'format':'tsv','fields':'accession,gene_primary,protein_name,rhea'})
    fetch(root,'UniProt_human_Rhea.tsv',url,online=online)
    url='https://www.rhea-db.org/rhea/?'+urllib.parse.urlencode({'query':'','columns':'rhea-id,equation,chebi-id,reaction-xref(KEGG),pubmed','format':'tsv'})
    fetch(root,'Rhea_reactions.tsv',url,online=online)
    fetch(root,'Rhea_directions.tsv','https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv',online=online)
    queries={
        'ChEBI_crossreferences.json':'PREFIX obo: <http://www.geneontology.org/formats/oboInOwl#> PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> SELECT ?chebi ?xref ?label ?formula WHERE { ?chebi obo:hasDbXref ?xref . FILTER(REGEX(STR(?xref),"HMDB|KEGG", "i")) OPTIONAL { ?chebi rdfs:label ?label } OPTIONAL { ?chebi <https://w3id.org/chemrof/generalized_empirical_formula> ?formula } }',
        'ChEBI_microspecies.json':'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> PREFIX owl: <http://www.w3.org/2002/07/owl#> SELECT DISTINCT ?chebi ?major WHERE { { ?chebi <https://w3id.org/chemrof/has_major_microspecies_at_pH7_3> ?major } UNION { ?chebi rdfs:subClassOf ?restriction . ?restriction owl:onProperty <https://w3id.org/chemrof/has_major_microspecies_at_pH7_3> ; owl:someValuesFrom ?major } }'
    }
    for name,query in queries.items():
        url='https://sparql.rhea-db.org/sparql?'+urllib.parse.urlencode({'query':query,'format':'json'})
        fetch(root,name,url,online=online)
    additional=sorted(set(resolve(root))-set(identifiers))
    for i in range(0,len(additional),10):
        fetch(root,f'KEGG_compounds_extra_{i//10:03d}.txt','https://rest.kegg.jp/get/' + '+'.join(additional[i:i+10]),online=online)
    for name,path in [('KEGG_compound_pathways.tsv','link/pathway/compound'),('KEGG_human_pathways.tsv','list/pathway/hsa')]:
        fetch(root,name,'https://rest.kegg.jp/'+path,online=online)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--online',action='store_true')
    args=parser.parse_args()
    main(Path(__file__).resolve().parents[1],args.online)
