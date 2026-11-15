from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import io
import re
import ssl
import urllib.request
import zipfile
import certifi
import pandas as pd

CATEGORIES={'GO_BP':'GO Biological Processes','GO_CC':'GO Cellular Components','GO_MF':'GO Molecular Functions','KEGG':'KEGG Pathway'}

def valid_payload(payload,filename,category):
    if filename=='all.zip':
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            if 'AnalysisReport.html' not in archive.namelist():
                return False
            table=pd.read_csv(archive.open('Enrichment_GO/GO_AllLists.csv'))
    else:
        table=pd.read_excel(io.BytesIO(payload),sheet_name='Enrichment')
    return set(table.Category.dropna())=={CATEGORIES[category]}

def main(root,category,session,online):
    folder=root/'raw/public'
    for filename in ['all.zip','metascape_result.xlsx']:
        path=folder/f'Metascape_{category}_{filename}'
        source=path.with_name(path.name+'.source.json')
        complete=path.exists()
        if complete:
            complete=valid_payload(path.read_bytes(),filename,category)
        if not complete:
            if not online:
                raise FileNotFoundError(path)
            url=f'https://metascape.org/gp_server/get_file?session_id={session}&file_name={filename}&rename=true'
            with urllib.request.urlopen(url,context=ssl.create_default_context(cafile=certifi.where()),timeout=180) as response:
                payload=response.read()
            assert payload.startswith(b'PK')
            assert valid_payload(payload,filename,category),'Report tables are not yet complete for the requested category'
            path.write_bytes(payload)
            source.write_text(json.dumps({'url':url,'session_id':session,'category':category,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(payload).hexdigest()},indent=2),encoding='utf-8')
        metadata=json.loads(source.read_text(encoding='utf-8'))
        assert hashlib.sha256(path.read_bytes()).hexdigest()==metadata['sha256']
    with zipfile.ZipFile(folder/f'Metascape_{category}_all.zip') as archive:
        report=re.sub(r'\s+',' ',re.sub('<[^>]+>',' ',archive.read('AnalysisReport.html').decode()))
        start=report.index('ontology sources')
        description=report[start:start+500]
        assert CATEGORIES[category] in description
        for other,name in CATEGORIES.items():
            if other!=category:
                assert name not in description
        (root/'reports'/f'Metascape_{category}_settings.txt').write_text(description,encoding='utf-8')
    print(category,description)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--category',choices=CATEGORIES,required=True)
    parser.add_argument('--session',required=True)
    parser.add_argument('--online',action='store_true')
    args=parser.parse_args()
    main(Path(__file__).resolve().parents[1],args.category,args.session,args.online)
