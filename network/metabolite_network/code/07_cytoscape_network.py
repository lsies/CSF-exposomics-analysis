from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import pandas as pd

def request(path, data=None):
    payload=None if data is None else json.dumps(data).encode()
    req=urllib.request.Request('http://localhost:1234/v1/'+path,data=payload,headers={'Content-Type':'application/json','Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=180) as response:
        return json.loads(response.read())

def main(root, export):
    destination=root/'raw/cytoscape'
    destination.mkdir(exist_ok=True)
    if export:
        identity=json.loads((destination/'network_identity.json').read_text())
        nid=identity['networkSUID']
        table=request(f'networks/{nid}/tables/defaultnode')
        (destination/'node_table.json').write_text(json.dumps(table,indent=2),encoding='utf-8')
        pd.DataFrame(table['rows']).to_csv(destination/'node_table.csv',index=False)
        result=request('commands/session/save',{'file':str(root/'results/Step6_network.cys')})
        (destination/'session_save.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
        print('Exported node table and session')
        return
    edges=pd.read_csv(root/'results/STRING_PPI_edges.csv')
    nodes=sorted(set(edges.gene_A)|set(edges.gene_B))
    graph={'data':{'name':'Step6_STRING_connected_209'},'elements':{'nodes':[{'data':{'id':x,'name':x,'gene_symbol':x}} for x in nodes],'edges':[{'data':{'source':r.gene_A,'target':r.gene_B,'interaction':'STRING functional','score':r.score}} for r in edges.itertuples()]}}
    identity=request('networks?collection=Step6',graph)
    (destination/'input_hash.json').write_text(json.dumps({'STRING_PPI_edges.csv':hashlib.sha256((root/'results/STRING_PPI_edges.csv').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
    (destination/'network_identity.json').write_text(json.dumps(identity,indent=2),encoding='utf-8')
    nid=identity['networkSUID']
    request('commands/layout/force-directed',{'network':f'SUID:{nid}','isDeterministic':'true','numIterations':'500'})
    parameters={'degreeCutoff':'2','nodeScoreCutoff':'0.2','kCore':'2','maxDepthFromStart':'100','haircut':'true','fluff':'false','includeLoops':'false','scope':'NETWORK'}
    result=request('commands/mcode/cluster',parameters)
    (destination/'MCODE_parameters.json').write_text(json.dumps(parameters,indent=2),encoding='utf-8')
    (destination/'MCODE_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (destination/'software_version.json').write_text(json.dumps(request('version'),indent=2),encoding='utf-8')
    print(json.dumps(identity),json.dumps(result)[:1200])

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--export',action='store_true')
    args=parser.parse_args()
    main(Path(__file__).resolve().parents[1],args.export)
