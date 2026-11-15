from datetime import datetime,timezone
import hashlib
import json
import ssl
import time
import urllib.parse
import urllib.request
import certifi

def fetch(root,name,url,payload=None,json_body=False,online=False):
    path=root/'raw/public'/name
    record=path.with_name(path.name+'.source.json')
    if path.exists() and record.exists():
        meta=json.loads(record.read_text(encoding='utf-8'))
        if meta['url']==url and meta['payload']==payload:
            assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['sha256'],f'Changed public snapshot: {name}'
            return path.read_text(encoding='utf-8-sig')
    if not online:
        raise FileNotFoundError(f'Missing public snapshot: {name}; use --online')
    headers={'User-Agent':'MetaboliteReactionNetwork/1.0'}
    body=None
    if payload is not None:
        body=(json.dumps(payload) if json_body else urllib.parse.urlencode(payload)).encode('utf-8')
        headers['Content-Type']='application/json' if json_body else 'application/x-www-form-urlencoded'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,data=body,headers=headers),timeout=180,context=ssl.create_default_context(cafile=certifi.where())) as response:
                data=response.read()
            text=data.decode('utf-8-sig')
            if '<html' in text[:500].lower():
                raise ValueError('Unexpected HTML response')
            path.write_bytes(data)
            record.write_text(json.dumps({'url':url,'payload':payload,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest()},indent=2),encoding='utf-8')
            time.sleep(0.4)
            return text
        except Exception:
            if attempt==2:
                raise
            time.sleep(2)
