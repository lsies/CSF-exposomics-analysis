from pathlib import Path
import argparse
import json
import pandas as pd
from public_sources import fetch

def main(root,online):
    frame=pd.read_csv(root/'data/protein_background.csv')
    genes=sorted(frame.loc[frame.selected,'gene_symbol'])
    base='https://version-12-0.string-db.org/api/json/'
    fetch(root,'STRING_version.json',base+'version',online=online)
    payload={'identifiers':'\r'.join(genes),'species':9606,'limit':1,'echo_query':1,'caller_identity':'MetaboliteReactionNetwork'}
    mapping=pd.DataFrame(json.loads(fetch(root,'STRING_mapping.json',base+'get_string_ids',payload,online=online)))
    mapping.to_csv(root/'results/STRING_identifier_mapping.csv',index=False)
    assert mapping.stringId.is_unique
    assert mapping.queryItem.is_unique
    payload={'identifiers':'\r'.join(sorted(mapping.stringId)),'species':9606,'required_score':400,'network_type':'functional','add_nodes':0,'caller_identity':'MetaboliteReactionNetwork'}
    network=pd.DataFrame(json.loads(fetch(root,'STRING_network.json',base+'network',payload,online=online)))
    aliases=dict(zip(mapping.stringId,mapping.queryItem))
    network['gene_A']=network.stringId_A.map(aliases)
    network['gene_B']=network.stringId_B.map(aliases)
    assert network.gene_A.notna().all() and network.gene_B.notna().all()
    network.to_csv(root/'results/STRING_PPI_edges.csv',index=False)
    print(json.dumps({'submitted_genes':len(genes),'STRING_mapped':len(mapping),'PPI_edges':len(network),'connected_genes':len(set(network.gene_A)|set(network.gene_B))}))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--online',action='store_true')
    args=parser.parse_args()
    main(Path(__file__).resolve().parents[1],args.online)
