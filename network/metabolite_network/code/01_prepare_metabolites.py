from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def main(root):
    frames=[]
    for comparison in ['LVSH','LVSM','MVSH']:
        frame=pd.read_csv(root/'raw/study'/f'{comparison}_diff_featurespeaks.txt',sep='\t')
        expected=frame.p_value.lt(0.05)&frame.VIP.ge(1)&(frame.FC.gt(1.2)|frame.FC.lt(1/1.2))
        assert expected.eq(frame.Significant.eq('yes')).all()
        assert np.allclose(frame.q_value,multipletests(frame.p_value,method='fdr_bh')[1],rtol=1e-10,atol=1e-12)
        frame['comparison']=comparison
        frame['physical_feature']=frame.Type+'_'+frame.MZ.map(lambda x:f'{x:.4f}')+'_'+frame.RT.map(lambda x:f'{x:.1f}')
        frame['selected']=expected
        frames.append(frame)
    full=pd.concat(frames,ignore_index=True)
    columns=['FeaturesPeak','physical_feature','Type','MZ','RT','Metabolites','ConfidenceScore','Formula','HMDB','kegg','comparison','FC','p_value','q_value','VIP','regulation','selected']
    full[columns].to_csv(root/'data/comparison_statistics.csv',index=False)
    annotations=full.drop_duplicates('FeaturesPeak')[columns[:10]].copy()
    selected_ids=set(full.loc[full.selected,'FeaturesPeak'])
    selected_features=set(full.loc[full.selected,'physical_feature'])
    annotations['selected_any_comparison']=annotations.FeaturesPeak.isin(selected_ids)
    annotations['physical_feature_selected']=annotations.physical_feature.isin(selected_features)
    annotations['selected_q05']=annotations.FeaturesPeak.isin(full.loc[full.selected & full.q_value.lt(.05),'FeaturesPeak'])
    assert annotations.selected_any_comparison.eq(annotations.selected_q05).all()
    annotations.to_csv(root/'data/metabolite_annotations.csv',index=False)
    edges=[]
    for row in annotations.itertuples():
        for identifier in sorted(set(re.findall(r'C\d{5}',str(row.kegg)))):
            edges.append({'FeaturesPeak':row.FeaturesPeak,'physical_feature':row.physical_feature,'metabolite_name':row.Metabolites,'formula':row.Formula,'KEGG':identifier,'selected':row.selected_any_comparison,'selected_q05':row.selected_q05})
    pd.DataFrame(edges).to_csv(root/'data/metabolite_KEGG_candidates.csv',index=False)
    pd.DataFrame({'KEGG':sorted({x['KEGG'] for x in edges})}).to_csv(root/'data/KEGG_query_ids.csv',index=False)
    rows=[]
    for comparison,group in full.groupby('comparison'):
        subset=group[group.selected]
        rows.append({'comparison':comparison,'selected_annotations':len(subset),'selected_physical_features':subset.physical_feature.nunique(),'selected_annotations_q05':int(subset.q_value.lt(.05).sum())})
    pd.DataFrame(rows).to_csv(root/'results/differential_selection_summary.csv',index=False)
    print(json.dumps({'annotations':len(annotations),'physical_features':annotations.physical_feature.nunique(),'selected_annotations':len(selected_ids),'selected_physical_features':len(selected_features),'KEGG_background':len({x['KEGG'] for x in edges})}))

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
