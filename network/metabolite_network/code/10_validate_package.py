from pathlib import Path
import ast
import hashlib
import io
import json
import tokenize
import pandas as pd
import networkx as nx
from PIL import Image

def main(root):
    checks=[]
    for path in [root/'run_analysis.py',*sorted((root/'code').glob('*.py'))]:
        source=path.read_text(encoding='utf-8')
        tree=ast.parse(source)
        assert not any(token.type==tokenize.COMMENT for token in tokenize.generate_tokens(io.StringIO(source).readline))
        assert not any(ast.get_docstring(node) for node in ast.walk(tree) if isinstance(node,(ast.Module,ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)))
        checks.append({'check':'clean_python_syntax','file':path.relative_to(root).as_posix(),'passed':True})
    for path in sorted((root/'raw/public').glob('*.source.json')):
        metadata=json.loads(path.read_text(encoding='utf-8'))
        original=path.with_name(path.name.removesuffix('.source.json'))
        assert hashlib.sha256(original.read_bytes()).hexdigest()==metadata['sha256'],str(original)
    checks.append({'check':'public_source_SHA256','file':'raw/public','passed':True})
    for path in sorted((root/'figures').rglob('*.tiff')):
        with Image.open(path) as image:
            assert image.info['dpi']==(600,600)
            assert image.width<=6300 and image.height<=8280
        checks.append({'check':'600_dpi_scaled_canvas','file':path.relative_to(root).as_posix(),'passed':True})
    typography=pd.read_csv(root/'data/figure_source/figure_typography.csv')
    assert len(typography)==6 and typography.font_scale.eq(1.8).all()
    checks.append({'check':'all_figure_text_enlarged_1_8_fold','file':'data/figure_source/figure_typography.csv','passed':True})
    for path in sorted((root/'figures').rglob('*.svg')):
        assert '<text' in path.read_text(encoding='utf-8')
    checks.append({'check':'editable_SVG_text','file':'figures','passed':True})
    coverage=pd.read_csv(root/'results/STRING_gene_coverage.csv')
    ppi=nx.read_graphml(root/'results/complete_PPI_network.graphml')
    original=pd.read_csv(root/'results/STRING_PPI_edges.csv')
    assert set(ppi)==set(coverage.loc[coverage.STRING_mapped,'gene_symbol'])
    assert {frozenset(e) for e in ppi.edges}=={frozenset((r.gene_A,r.gene_B)) for r in original.itertuples()}
    local=nx.read_graphml(root/'results/KYNU_local_biochemistry.graphml')
    assert set(local.neighbors('KYNU'))=={'RHEA:16813','RHEA:25143','hsa00380'}
    assert (root/'figures/Fig6_1A_PPI_full.svg').is_file()
    assert (root/'figures/Fig6_1B_KYNU_local.svg').is_file()
    checks.append({'check':'complete_PPI_and_separate_KYNU_sources','file':'results/complete_PPI_network.graphml','passed':True})
    evidence=pd.read_csv(root/'results/KP_expression_evidence_matrix.csv')
    assert evidence[evidence.gene_symbol.eq('KYNU')].status.eq('not_testable').all()
    checks.append({'check':'KYNU_expression_missingness','file':'results/KP_expression_evidence_matrix.csv','passed':True})
    enrichment=pd.read_csv(root/'results/Metascape_enrichment_all_returned_terms.csv')
    assert not enrichment.duplicated(['Category','GO']).any()
    assert set(enrichment.analysis_family)=={'GO_BP','GO_CC','GO_MF','KEGG'}
    assert enrichment.groupby('analysis_family').Category.nunique().eq(1).all()
    assert enrichment.q_value.ge(enrichment.p_value-1e-12).all()
    assert 'KYNU' in enrichment.loc[enrichment.GO.eq('hsa00380'),'Hits'].iloc[0].split('|')
    checks.append({'check':'category_specific_FDR_and_KYNU_membership','file':'results/Metascape_enrichment_all_returned_terms.csv','passed':True})
    annotations=pd.read_csv(root/'data/metabolite_annotations.csv')
    assert annotations.selected_any_comparison.eq(annotations.selected_q05).all()
    checks.append({'check':'differential_union_supported_by_q05','file':'data/metabolite_annotations.csv','passed':True})
    pd.DataFrame(checks).to_csv(root/'reports/validation_checks.csv',index=False)
    print(f'Passed {len(checks)} validation checks')

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
