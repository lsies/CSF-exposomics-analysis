import json
import re
import pandas as pd

def resolve(root):
    records={}
    for row in json.loads((root/'raw/public/ChEBI_crossreferences.json').read_text(encoding='utf-8'))['results']['bindings']:
        chebi=row['chebi']['value'].split('_')[-1]
        value=records.setdefault(chebi,{'names':set(),'formulas':set(),'HMDB':set(),'KEGG':set()})
        value['names'].add(row.get('label',{}).get('value',''))
        value['formulas'].add(row.get('formula',{}).get('value',''))
        xref=row['xref']['value']
        if re.search(r'HMDB\d+',xref,re.I):
            hmdb=re.search(r'HMDB(\d+)',xref,re.I).group(1)
            value['HMDB'].add('HMDB'+hmdb.zfill(7))
        if re.search(r'C\d{5}',xref):
            value['KEGG'].add(re.search(r'C\d{5}',xref).group())
    lookup={}
    for chebi,record in records.items():
        for hmdb in record['HMDB']:
            lookup.setdefault(hmdb,set()).add(chebi)
    annotations=pd.read_csv(root/'data/metabolite_annotations.csv').fillna('')
    candidates=[]
    for row in annotations.itertuples():
        provided=set(re.findall(r'C\d{5}',str(row.kegg)))
        chebis=set()
        for identifier in re.findall(r'HMDB(\d+)',str(row.HMDB),re.I):
            chebis.update(lookup.get('HMDB'+identifier.zfill(7),set()))
        chebis={c for c in chebis if row.Formula in records[c]['formulas']}
        common={'FeaturesPeak':row.FeaturesPeak,'physical_feature':row.physical_feature,'metabolite_name':row.Metabolites,'formula':row.Formula,'selected':row.selected_any_comparison,'selected_q05':row.selected_q05,'supplier_KEGG':';'.join(sorted(provided))}
        if chebis:
            for chebi in sorted(chebis):
                record=records[chebi]
                identifiers=sorted(record['KEGG']) or ['']
                for identifier in identifiers:
                    candidates.append({**common,'CompoundID':identifier or 'CHEBI:'+chebi,'KEGG':identifier,'ChEBI':chebi,'reference_name':';'.join(sorted(record['names'])),'reference_formula':row.Formula,'identity_source':'HMDB_ChEBI_crossreference','supplier_id_discordant':bool(provided and identifier and identifier not in provided)})
        else:
            for identifier in sorted(provided):
                candidates.append({**common,'CompoundID':identifier,'KEGG':identifier,'ChEBI':'','reference_name':'','reference_formula':'','identity_source':'supplier_KEGG','supplier_id_discordant':False})
    frame=pd.DataFrame(candidates).drop_duplicates()
    frame.to_csv(root/'data/compound_candidates.csv',index=False)
    return sorted(set(frame.KEGG)-{''})
