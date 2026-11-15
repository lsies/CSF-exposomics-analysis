from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from figure_style import COLORS, save
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot_kynu(root, colors=COLORS, export=save):
    result = root / 'results'
    sources = root / 'data/figure_source'
    kp = pd.read_csv(result / 'KP_core_metabolite_edges.csv')
    kp = kp.loc[kp.gene_symbol.eq('KYNU')].drop_duplicates(['Rhea', 'FeaturesPeak'])
    assert set(kp.Rhea.astype(int)) == {16813, 25143}
    assert kp.metabolite_name.eq('Alanine').all()
    pathway_edges = pd.read_csv(result / 'KP_KEGG_annotation_edges.csv')
    assert pathway_edges.gene_symbol.eq('KYNU').any()
    local = nx.Graph()
    local.add_node('Alanine', node_type='metabolite_annotation', label='Alanine', FeaturesPeak=kp.FeaturesPeak.iloc[0], physical_feature=kp.physical_feature.iloc[0])
    local.add_node('KYNU', node_type='human_gene', label='KYNU', UniProt='Q16719')
    local.add_node('hsa00380', node_type='KEGG_pathway', label='Tryptophan metabolism', pathway_q_value=float(pathway_edges.pathway_q_value.iloc[0]))
    local_positions = {'Alanine': (.10, .60), 'RHEA:16813': (.38, .78), 'RHEA:25143': (.38, .42), 'KYNU': (.65, .60), 'hsa00380': (.90, .60)}
    for r in kp.itertuples():
        reaction = f'RHEA:{int(r.Rhea)}'
        local.add_node(reaction, node_type='biochemical_reaction', label=reaction, equation=r.equation)
        local.add_edge('Alanine', reaction, edge_type='reaction_participant', source_database='Rhea')
        local.add_edge(reaction, 'KYNU', edge_type='catalysis_annotation', source_database='Rhea/UniProt')
    local.add_edge('KYNU', 'hsa00380', edge_type='KEGG_annotation', source_database='KEGG')
    fig = plt.figure(figsize=(7, 1.8))
    ax = fig.add_axes([.025, .16, .95, .82])
    ax.set_xlim(.015, 1.005)
    ax.set_ylim(.34, 1.03)
    ax.axis('off')
    for u, v, d in local.edges(data=True):
        (x0, y0), (x1, y1) = local_positions[u], local_positions[v]
        ax.plot([x0, x1], [y0, y1], color=colors['blue'] if d['edge_type'] == 'KEGG_annotation' else colors['orange'], lw=2.16, ls='--' if d['edge_type'] == 'KEGG_annotation' else '-', zorder=1)
    for name, (x, y) in local_positions.items():
        color = colors['purple'] if name == 'Alanine' else colors['orange'] if name == 'KYNU' else colors['blue'] if name == 'hsa00380' else '#8B9CAA'
        marker = 'D' if name == 'hsa00380' else 's' if name.startswith('RHEA') else 'o'
        ax.scatter(x, y, s=300 if name == 'KYNU' else 120, color=color, marker=marker, edgecolors='none', zorder=3)
        label = 'Tryptophan\nmetabolism' if name == 'hsa00380' else name
        ax.text(x, y + .10, label, ha='center', va='bottom', fontsize=8, weight='bold' if name == 'KYNU' else 'normal', color=colors['orange'] if name == 'KYNU' else colors['ink'])
    handles = [Line2D([], [], marker='o', ls='', color=colors['purple'], markeredgewidth=0, label='Study annotation'), Line2D([], [], color=colors['orange'], lw=2.16, label='Biochemical connection'), Line2D([], [], color=colors['blue'], lw=2.16, ls='--', label='KEGG membership')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .012), ncol=3, fontsize=7, columnspacing=.75, handletextpad=.35, handlelength=1.2, borderpad=0, borderaxespad=0, labelspacing=.15)
    export(fig, root, 'Fig6_1B_KYNU_local')
    nx.set_node_attributes(local, {n: {'x': float(p[0]), 'y': float(p[1])} for n, p in local_positions.items()})
    nx.write_graphml(local, result / 'KYNU_local_biochemistry.graphml')
    pd.DataFrame([{'node': n, **d} for n, d in local.nodes(data=True)]).to_csv(sources / 'KYNU_local_nodes.csv', index=False)
    pd.DataFrame([{'source': u, 'target': v, **d} for u, v, d in local.edges(data=True)]).to_csv(sources / 'KYNU_local_edges.csv', index=False)
    kp.to_csv(sources / 'KYNU_reaction_equations.csv', index=False)
    full = nx.read_graphml(result / 'complete_PPI_network.graphml')
    integrated = nx.compose(full, local)
    for row in pathway_edges.itertuples():
        integrated.add_edge(row.gene_symbol, 'hsa00380', edge_type='KEGG_annotation', source_database='KEGG')
    nx.write_graphml(integrated, result / 'integrated_KP_network.graphml')
    pd.DataFrame([{'node': n, **d} for n, d in integrated.nodes(data=True)]).to_csv(sources / 'integrated_network_nodes.csv', index=False)
    pd.DataFrame([{'source': u, 'target': v, **d} for u, v, d in integrated.edges(data=True)]).to_csv(sources / 'integrated_network_edges.csv', index=False)
    pd.DataFrame([{'node': n, 'x': d['x'], 'y': d['y'], 'layout_panel': d.get('layout_panel', 'local_biochemistry')} for n, d in integrated.nodes(data=True)]).to_csv(sources / 'integrated_network_positions.csv', index=False)


if __name__ == '__main__':
    plot_kynu(Path(__file__).resolve().parents[1])
