import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as pe


def place_nearby_labels(ax, graph, nodes, positions, focus, fills, fontsize):
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    inverse = ax.transData.inverted()
    origin = inverse.transform((0, 0))
    boxes = []
    leaders = []
    records = []
    priority = sorted(focus, key=lambda g: (g != 'KYNU', -sum(np.linalg.norm(positions[g] - positions[h]) < .25 for h in focus), g))
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    for gene in priority:
        weight = 'bold' if gene == 'KYNU' or nodes.loc[gene, 'consensus_hub'] else 'normal'
        prop = FontProperties(family='Arial', size=fontsize * 1.8 / 1.5, weight=weight)
        width, height, descent = renderer.get_text_width_height_descent(gene, prop, ismath=False)
        dx, dy = inverse.transform((width + 8, height + descent + 8)) - origin
        x, y = positions[gene]
        candidates = []
        for radius in [.07, .10, .14, .19, .25, .32, .40, .50, .65, .80]:
            for angle in np.linspace(0, 2 * np.pi, 32, endpoint=False):
                center = np.array([x + radius * np.cos(angle), y + radius * np.sin(angle)])
                cx, cy = center
                box = (cx - dx / 2, cy - dy / 2, cx + dx / 2, cy + dy / 2)
                if box[0] < xlim[0] + .015 or box[2] > xlim[1] - .015 or box[1] < ylim[0] + .015 or box[3] > ylim[1] - .015:
                    continue
                overlap = sum(max(0, min(box[2], q[2]) - max(box[0], q[0])) * max(0, min(box[3], q[3]) - max(box[1], q[1])) for q in boxes)
                covered_focus = sum(box[0] - .025 < positions[h][0] < box[2] + .025 and box[1] - .025 < positions[h][1] < box[3] + .025 for h in focus)
                covered_other = sum(box[0] - .015 < positions[h][0] < box[2] + .015 and box[1] - .015 < positions[h][1] < box[3] + .015 for h in graph if h not in focus)
                crossings = sum(segments_cross(np.array([x, y]), center, a, b) for a, b in leaders)
                leader_text = sum(segment_hits_box(np.array([x, y]), center, q) for q in boxes)
                text_leader = sum(segment_hits_box(a, b, box) for a, b in leaders)
                penalty = overlap * 1e12 + (leader_text + text_leader) * 1e10 + covered_focus * 1e7 + covered_other * 1e6 + crossings * 1e8 + radius * 75
                candidates.append((penalty, radius, box, center))
        _, radius, box, center = min(candidates, key=lambda c: c[0])
        collisions = sum(box[0] < q[2] and box[2] > q[0] and box[1] < q[3] and box[3] > q[1] for q in boxes)
        assert collisions == 0, gene
        boxes.append(box)
        leaders.append((np.array([x, y]), center))
        text = ax.annotate(gene, xy=(x, y), xytext=center, ha='center', va='center', fontsize=fontsize, weight=weight, color=fills[gene], zorder=10, arrowprops={'arrowstyle': '-', 'lw': .55, 'color': fills[gene], 'shrinkA': 3, 'shrinkB': 4})
        text.set_path_effects([pe.withStroke(linewidth=2, foreground='white')])
        records.append({'gene_symbol': gene, 'node_x': x, 'node_y': y, 'label_x': center[0], 'label_y': center[1], 'leader_length_layout_units': radius, 'text_overlap_count': collisions})
    for index, record in enumerate(records):
        record['leader_text_overlap_count'] = sum(segment_hits_box(*leaders[index], box) for j, box in enumerate(boxes) if j != index)
        record['leader_crossing_count'] = sum(segments_cross(*leaders[index], *line) for j, line in enumerate(leaders) if j != index)
    assert not any(r['leader_text_overlap_count'] or r['leader_crossing_count'] for r in records)
    return records


def separate_nodes(ax, positions, sizes):
    names = list(positions)
    coordinates = ax.transData.transform(np.array([positions[g] for g in names]))
    radii = np.array([np.sqrt(sizes[g]) * ax.figure.dpi / 144 / 1.5 for g in names])
    threshold = radii[:, None] + radii[None, :] + 1.5
    np.fill_diagonal(threshold, 0)
    for iteration in range(400):
        difference = coordinates[:, None, :] - coordinates[None, :, :]
        distance = np.linalg.norm(difference, axis=2)
        overlap = np.maximum(threshold - distance, 0)
        if overlap.max() < .01:
            break
        coordinates += (difference / np.maximum(distance[:, :, None], 1e-9) * overlap[:, :, None] * .52).sum(axis=1)
    assert overlap.max() < .01
    return dict(zip(names, ax.transData.inverted().transform(coordinates)))


def segment_hits_box(a, b, box):
    x0, y0, x1, y1 = box
    corners = [np.array([x0, y0]), np.array([x1, y0]), np.array([x1, y1]), np.array([x0, y1])]
    return any(segments_cross(a, b, corners[i], corners[(i + 1) % 4]) for i in range(4)) or (x0 < a[0] < x1 and y0 < a[1] < y1)


def segments_cross(a, b, c, d):
    def side(p, q, r):
        u, v = q - p, r - p
        return u[0] * v[1] - u[1] * v[0]
    return side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0

