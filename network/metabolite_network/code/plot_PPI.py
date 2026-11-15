from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from figure_style import COLORS, save
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from ppi_labels import place_nearby_labels, separate_nodes


def plot_ppi(root, nodes, graph, colors, save):
    result = root / 'results'
    sources = root / 'data/figure_source'
    coverage = pd.read_csv(result / 'STRING_gene_coverage.csv')
    isolates = sorted(coverage.loc[coverage.STRING_mapped & ~coverage.connected_network, 'gene_symbol'])
    positions = nx.kamada_kawai_layout(graph, weight=None)
    angle = np.arctan2(positions['KYNU'][1], positions['KYNU'][0])
    rotation = np.pi / 4 - angle
    matrix = np.array([[np.cos(rotation), -np.sin(rotation)], [np.sin(rotation), np.cos(rotation)]])
    positions = {g: matrix @ p for g, p in positions.items()}
    positions = {g: p / max(np.linalg.norm(v) for v in positions.values()) for g, p in positions.items()}
    positions = {g: p * max(np.linalg.norm(p), 1e-8) ** (-.30) for g, p in positions.items()}
    enrichment = pd.read_csv(result / 'Metascape_enrichment_all_returned_terms.csv')
    trp_rows = []
    for row in enrichment.loc[enrichment.GO.isin(['hsa00380', 'GO:0015827'])].itertuples():
        for gene in row.Hits.split('|'):
            if gene in graph:
                trp_rows.append({'gene_symbol': gene, 'term_id': row.GO, 'description': row.Description, 'category': row.Category, 'term_q_value': row.q_value})
    pd.DataFrame(trp_rows).to_csv(sources / 'PPI_tryptophan_annotations.csv', index=False)
    trp_genes = set(r['gene_symbol'] for r in trp_rows)
    ordered = sorted(graph, key=lambda g: (g == 'KYNU', g in trp_genes, bool(nodes.loc[g, 'consensus_hub']), bool(nodes.loc[g, 'exposure_responsive'])))
    fills = {g: colors['orange'] if g == 'KYNU' else colors['blue'] if g in trp_genes else colors['pink'] if nodes.loc[g, 'consensus_hub'] else colors['teal'] if nodes.loc[g, 'exposure_responsive'] else '#AEBBC6' for g in ordered}
    sizes = {g: 170 if g == 'KYNU' else 105 if nodes.loc[g, 'consensus_hub'] else 85 if g in trp_genes else 75 if nodes.loc[g, 'exposure_responsive'] else 15 + 1.5 * nodes.loc[g, 'Degree'] for g in ordered}
    sizes = {g: value * 1.4 for g, value in sizes.items()}
    focus = [g for g in ordered if g in trp_genes or nodes.loc[g, 'consensus_hub'] or nodes.loc[g, 'exposure_responsive']]
    fig = plt.figure(figsize=(7, 6.1))
    ax = fig.add_axes([.012, .15, .976, .843])
    coordinates = np.array(list(positions.values()))
    ax.set_xlim(coordinates[:, 0].min() - .065, coordinates[:, 0].max() + .065)
    ax.set_ylim(coordinates[:, 1].min() - .055, coordinates[:, 1].max() + .055)
    ax.set_aspect('auto')
    ax.axis('off')
    positions = separate_nodes(ax, positions, sizes)
    nx.draw_networkx_edges(graph, positions, ax=ax, width=.42, alpha=.27, edge_color='#8B9CAA')
    nx.draw_networkx_nodes(graph, positions, ax=ax, nodelist=ordered, node_color=[fills[g] for g in ordered], node_size=[sizes[g] for g in ordered], edgecolors='none')
    label_records = place_nearby_labels(ax, graph, nodes, positions, focus, fills, 8 * 1.7)
    pd.DataFrame(label_records).to_csv(sources / 'PPI_standalone_label_positions.csv', index=False)
    isolated_positions = {gene: (.035 + (i % 3) * .345, .75 - (i // 3) * .50) for i, gene in enumerate(isolates)}
    handles = [Line2D([], [], marker='o', ls='', color=c, markersize=6 * np.sqrt(1.4) * 2.7, markeredgewidth=0, label=label) for label, c in [('Other mapped genes', '#AEBBC6'), ('Consensus hubs', colors['pink']), ('KYNU', colors['orange']), ('Exposure-responsive', colors['teal']), ('TRP metabolism / transport', colors['blue'])]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .008), ncol=2, fontsize=7 * 1.7, columnspacing=.8, handletextpad=.45, handlelength=1.6, borderpad=0, borderaxespad=0, labelspacing=.25)
    full = graph.copy()
    full.add_nodes_from(isolates)
    nx.set_edge_attributes(full, 'STRING_association', 'edge_type')
    for gene in full:
        connected = gene in graph
        pos = positions[gene] if connected else isolated_positions[gene]
        full.nodes[gene].update({'node_type': 'human_gene', 'label': gene, 'x': float(pos[0]), 'y': float(pos[1]), 'connected_network': connected, 'displayed_in_figure': connected, 'tryptophan_related': gene in trp_genes, 'layout_panel': 'connected' if connected else 'isolates', 'fill_color': fills.get(gene, '#AEBBC6'), 'node_size': float(sizes.get(gene, 28 * 1.4))})
    nx.write_graphml(full, result / 'complete_PPI_network.graphml')
    pd.DataFrame([{'gene_symbol': g, **d} for g, d in full.nodes(data=True)]).to_csv(sources / 'full_PPI_nodes.csv', index=False)
    pd.DataFrame([{'gene_symbol': g, 'x': float(p[0]), 'y': float(p[1])} for g, p in positions.items()]).to_csv(sources / 'STRING_positions.csv', index=False)
    pd.DataFrame([{'source': u, 'target': v, **d} for u, v, d in full.edges(data=True)]).to_csv(sources / 'full_PPI_edges.csv', index=False)
    save(fig, root, 'Fig6_1A_PPI_full')

def main(root):
    nodes = pd.read_csv(root / 'results/network_node_statistics.csv').set_index('gene_symbol')
    edges = pd.read_csv(root / 'results/STRING_PPI_edges.csv')
    graph = nx.from_pandas_edgelist(edges, 'gene_A', 'gene_B', edge_attr='score')
    plot_ppi(root, nodes, graph, COLORS, save)


if __name__ == '__main__':
    main(Path(__file__).resolve().parents[1])
