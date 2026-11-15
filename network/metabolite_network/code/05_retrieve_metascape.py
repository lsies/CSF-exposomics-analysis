from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import ssl
import urllib.request
import zipfile
import certifi

def main(root,online):
    session='tfcfqjle6'
    for filename in ['all.zip','metascape_result.xlsx']:
        path=root/'raw/public'/('Metascape_'+filename)
        sidecar=path.with_name(path.name+'.source.json')
        url=f'https://metascape.org/gp_server/get_file?session_id={session}&file_name={filename}&rename=true'
        if not path.exists():
            if not online:
                raise FileNotFoundError(path)
            with urllib.request.urlopen(url,context=ssl.create_default_context(cafile=certifi.where()),timeout=180) as response:
                data=response.read()
            assert data.startswith(b'PK')
            path.write_bytes(data)
            sidecar.write_text(json.dumps({'url':url,'session_id':session,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest()},indent=2),encoding='utf-8')
        metadata=json.loads(sidecar.read_text(encoding='utf-8'))
        assert hashlib.sha256(path.read_bytes()).hexdigest()==metadata['sha256']
    with zipfile.ZipFile(root/'raw/public/Metascape_all.zip') as zipped:
        assert 'AnalysisReport.html' in zipped.namelist()
    print('Original Metascape archive verified')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--online',action='store_true')
    args=parser.parse_args()
    main(Path(__file__).resolve().parents[1],args.online)
