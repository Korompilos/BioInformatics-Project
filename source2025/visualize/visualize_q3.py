import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

from config import (
    HMM_PROFILE_BEFORE, HMM_PROFILE_AFTER,
    VISUALIZATION_DIR,
    ensure_dirs
)

DEFAULT_JSON = str(HMM_PROFILE_AFTER)
DEFAULT_OUTDIR = str(VISUALIZATION_DIR)

edge_threshold = 0.001
top_k_outgoing = None
label_top_k = 3
forward_only = True
max_jump = None
prob_decimals = 2

layer_gap = 7.0
h_pad = 220.0
node_size = 400
type_x_offset = 0.45
label_offset = 0.45
label_min_dist = 0.38


def state_num(s):
    if s in ('Start', 'Begin'):
        return -1
    if s == 'End':
        return None
    d = [c for c in s if c.isdigit()]
    return int("".join(d)) if d else 0


def node_type(s):
    if s in ('Start', 'Begin'):
        return 'Start'
    if s == 'End':
        return 'End'
    p = s[:1]
    return p if p in ('M', 'I', 'D') else 'Other'


def sort_key(s):
    t = node_type(s)
    r = {'Start': -1, 'M': 0, 'I': 1, 'D': 2, 'End': 9999, 'Other': 5000}.get(t, 5000)
    n = state_num(s)
    if n is None:
        n = 10 ** 9
    return (n, r, s)


def build_all_states(emissions, transitions):
    nodes = set()
    if isinstance(emissions, dict):
        nodes.update(emissions.keys())
    if isinstance(transitions, dict):
        nodes.update(transitions.keys())
        for d in transitions.values():
            if isinstance(d, dict):
                nodes.update(d.keys())
    return sorted(list(nodes), key=sort_key)


def build_index_map(states):
    nums = [state_num(s) for s in states if state_num(s) is not None]
    mx = max(nums) if nums else 0
    idx = {}
    for s in states:
        if s in ('Start', 'Begin'):
            idx[s] = -1
        elif s == 'End':
            idx[s] = mx + 1
        else:
            idx[s] = state_num(s) if state_num(s) is not None else 0
    return idx


def layer_of(s):
    return {'Start': 0, 'M': 1, 'I': 2, 'D': 3, 'End': 4, 'Other': 5}.get(node_type(s), 5)


def plot_heatmap(mat, xlabels, ylabels, title, outpath, yaxis='From state'):
    if mat is None or mat.size == 0:
        return False
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(111)
    im = ax.imshow(mat, aspect='auto', interpolation='nearest')
    ax.set_title(title)
    ax.set_xlabel('Symbols')
    ax.set_ylabel(yaxis)
    nx = len(xlabels)
    ny = mat.shape[0]
    xt = list(range(0, nx, max(1, nx // 25)))
    yt = list(range(0, ny, max(1, ny // 18)))
    ax.set_xticks(xt)
    ax.set_xticklabels([xlabels[i] for i in xt], rotation=90)
    ax.set_yticks(yt)
    ax.set_yticklabels([ylabels[i] for i in yt])
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Probability', rotation=90, va='center')
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout(pad=1.2)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return True


def draw_graph(states, transitions, index_map, out_png):
    if not states or not isinstance(transitions, dict):
        return False
    import matplotlib.patheffects as pe
    xs = [index_map[s] for s in states]
    xmin, xmax = min(xs), max(xs)
    span = max(1, xmax - xmin + 2)
    scale_x = h_pad / span
    tx = {'Start': 0.0, 'M': -type_x_offset, 'I': 0.0, 'D': type_x_offset, 'End': 0.0, 'Other': 0.0}
    counts = {}
    for s in states:
        k = (index_map[s], layer_of(s))
        counts[k] = counts.get(k, 0) + 1
    offsets = {k: list(np.linspace(-0.18, 0.18, max(1, c))) for k, c in counts.items()}
    used = {}
    pos = {}
    for s in states:
        k = (index_map[s], layer_of(s))
        i = used.get(k, 0)
        used[k] = i + 1
        base_x = (k[0] - xmin + 1) * scale_x
        x = base_x + tx.get(node_type(s), 0.0)
        y = k[1] * layer_gap + offsets[k][i]
        pos[s] = (x, y)
    fig = plt.figure(figsize=(40, 15))
    ax = fig.add_subplot(111)
    types = [node_type(s) for s in states]
    lv = {'Start': 0, 'M': 1, 'I': 2, 'D': 3, 'End': 4, 'Other': 5}
    cols = [lv.get(t, 5) for t in types]
    ax.scatter([pos[s][0] for s in states], [pos[s][1] for s in states],
               s=node_size, c=cols, cmap='viridis', zorder=3)
    for s in states:
        ax.text(pos[s][0], pos[s][1], s, ha='center', va='center',
                fontsize=12, color='white', weight='bold', zorder=4)
    for i, name in enumerate(['Start', 'M', 'I', 'D', 'End', 'Other']):
        yy = i * layer_gap
        ax.hlines(yy, 0, (xmax - xmin + 2) * scale_x + 1, linestyles='dotted',
                  linewidth=0.5, alpha=0.22, zorder=1)
        ax.text(0, yy + 0.35, name, fontsize=12, ha='left', va='bottom',
                alpha=0.5, zorder=4)
    edges = []
    for src in states:
        row = transitions.get(src, {})
        if not isinstance(row, dict):
            continue
        items = []
        for dst, p in row.items():
            if not isinstance(p, (int, float)) or dst not in pos:
                continue
            if forward_only:
                if index_map[dst] < index_map[src]:
                    continue
                if max_jump is not None and (index_map[dst] - index_map[src]) > max_jump and dst != 'End':
                    continue
            items.append((dst, float(p)))
        items.sort(key=lambda x: x[1], reverse=True)
        keep = [(src, d, p) for d, p in items if p >= edge_threshold]
        if top_k_outgoing and len(keep) > top_k_outgoing:
            keep = keep[:top_k_outgoing]
        edges.extend(keep)
    if not edges:
        plt.close(fig)
        return False
    probs = [p for (_, _, p) in edges]
    pmin, pmax = min(probs), max(probs)

    def w(p):
        if pmax == pmin:
            return 2.2
        return 1.0 + 4.8 * (p - pmin) / (pmax - pmin)

    def a(p):
        return 0.35 + 0.55 * (p - pmin) / (pmax - pmin + 1e-12)

    labels = []
    labels_drawn = {}
    for src, dst, p in edges:
        x1, y1 = pos[src]
        x2, y2 = pos[dst]
        rad = 0.11 if abs(y2 - y1) < 0.3 else 0.0
        if src == dst:
            rad = 0.5
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                     connectionstyle=f"arc3,rad={rad}",
                                     arrowstyle='-|>', mutation_scale=13, lw=w(p),
                                     alpha=a(p), color='black', zorder=2))
        labels_drawn[src] = labels_drawn.get(src, 0)
        if labels_drawn[src] < label_top_k:
            mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            dx, dy = (x2 - x1), (y2 - y1)
            nx, ny = (-dy, dx)
            nrm = np.hypot(nx, ny) + 1e-12
            nx, ny = nx / nrm, ny / nrm
            px, py = mx + nx * label_offset, my + ny * label_offset
            text = ax.text(px, py, f'{p:.{prob_decimals}f}', fontsize=12, alpha=0.98, zorder=5,
                           path_effects=[pe.withStroke(linewidth=2.0, foreground='white')])
            labels.append({
                "artist": text,
                "mx": mx, "my": my,
                "nx": nx, "ny": ny,
                "px": px, "py": py,
                "src": src, "dst": dst, "p": p
            })
            labels_drawn[src] += 1
    if labels:
        fig.canvas.draw()
        inv = ax.transData.inverted()
        rend = fig.canvas.get_renderer()

        def data_size(ta):
            bb = ta.get_window_extent(renderer=rend)
            (x0, y0) = inv.transform((bb.x0, bb.y0))
            (x1, y1) = inv.transform((bb.x1, bb.y1))
            return abs(x1 - x0), abs(y1 - y0)

        bin_w = scale_x * 0.9
        bins = {}
        for i, L in enumerate(labels):
            key = int(L["px"] / bin_w)
            bins.setdefault(key, []).append(i)
        for key, idxs in bins.items():
            idxs.sort(key=lambda i: (labels[i]["py"], -labels[i]["p"]))
            y_prev = None
            for j, i in enumerate(idxs):
                ta = labels[i]["artist"]
                w, h = data_size(ta)
                min_gap = max(label_min_dist * 0.9, h * 1.15)
                x, y = labels[i]["px"], labels[i]["py"]
                if y_prev is not None and y < y_prev + min_gap:
                    y = y_prev + min_gap
                y_prev = y
                ta.set_position((x, y))
        for key, idxs in bins.items():
            for kpos, i in enumerate(idxs):
                ta = labels[i]["artist"]
                x, y = ta.get_position()
                ta.set_position((x + (kpos % 2) * 0.05, y))
    ax.set_xlim(-0.5, (xmax - xmin + 2) * scale_x + 0.5)
    ax.set_ylim(-layer_gap * 0.6, layer_gap * 4.6)
    ax.axis('off')
    ax.set_title('Transition Graph')
    fig.tight_layout(pad=1.0)
    fig.savefig(out_png, dpi=320, bbox_inches='tight')
    plt.close(fig)
    return True


def run(json_path: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    emissions = data.get('emissions')
    transitions = data.get('transitions')
    states = []
    symbols = []
    if isinstance(emissions, dict) and emissions:
        states = sorted(list(emissions.keys()), key=sort_key)
        symset = set()
        for st in states:
            v = emissions.get(st, {})
            if isinstance(v, dict):
                symset.update(v.keys())
        symbols = sorted(list(symset))
        rows = []
        for st in states:
            v = emissions.get(st, {})
            row = []
            for sym in symbols:
                x = v.get(sym, np.nan)
                row.append(float(x) if isinstance(x, (int, float)) else np.nan)
            rows.append(row)
        mat_em = np.array(rows, dtype=float) if rows else None
        if mat_em is not None and mat_em.size:
            plot_heatmap(mat_em, symbols, states, 'Emission Probabilities',
                         outdir / 'emissions_heatmap.png', yaxis='State')
    if isinstance(transitions, dict) and transitions:
        all_states = build_all_states(emissions, transitions)
        index_map = build_index_map(all_states)
        draw_graph(all_states, transitions, index_map, outdir / 'transitions_graph.png')
        src_states = sorted(list(transitions.keys()), key=sort_key)
        dst_set = set()
        for s in src_states:
            t = transitions.get(s, {})
            if isinstance(t, dict):
                dst_set.update(t.keys())
        dst_states = sorted(list(dst_set), key=sort_key)
        rows = []
        for s in src_states:
            t = transitions.get(s, {})
            row = []
            for d in dst_states:
                x = t.get(d, np.nan)
                row.append(float(x) if isinstance(x, (int, float)) else np.nan)
            rows.append(row)
        mat_tr = np.array(rows, dtype=float) if rows else None
        if mat_tr is not None and mat_tr.size:
            plot_heatmap(mat_tr, dst_states, src_states, 'Transition Probabilities',
                         outdir / 'transitions_heatmap.png', yaxis='From state')


def _tag_from_json_name(p: Path) -> str:
    name = p.name.lower()
    if "before" in name:
        return "before"
    if "after" in name:
        return "after"
    return "custom"


def parse_args():
    p = argparse.ArgumentParser(description="Visualize a Profile HMM (emissions, transitions, graph).")
    p.add_argument(
        "json_path",
        nargs="?",
        default=DEFAULT_JSON,
        help="Path to HMM profile JSON. You can also pass 'before' or 'after'."
    )
    p.add_argument(
        "outdir",
        nargs="?",
        default=DEFAULT_OUTDIR,
        help="Output directory for generated figures (default: config VIZ_DIR)."
    )
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs()
    jp = str(args.json_path).strip().lower()
    if jp == "before":
        json_path = Path(HMM_PROFILE_BEFORE)
    elif jp == "after":
        json_path = Path(HMM_PROFILE_AFTER)
    else:
        json_path = Path(args.json_path)
    user_supplied_default = str(args.outdir) == DEFAULT_OUTDIR
    outdir = Path(args.outdir)
    if user_supplied_default:
        tag = _tag_from_json_name(json_path)
        outdir = Path(DEFAULT_OUTDIR) / tag
    run(str(json_path), outdir)
