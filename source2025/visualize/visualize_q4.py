import argparse
import json
import math
import random
import sys
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt

from config import RESULTS_DIR, VIZ_Q4, ensure_dirs, Q4_RESULTS_PATH


def resolve_results(user_arg: str) -> Path:
    p = Path(user_arg)
    if p.is_absolute() and p.is_file():
        return p
    cand1 = RESULTS_DIR / p.name
    if cand1.is_file():
        return cand1
    cand2 = Path.cwd() / p
    if cand2.is_file():
        return cand2
    tried = [str(cand1), str(cand2), str(p)]
    print("Could not find q4 results JSON. Tried:\n" + "\n".join(f" - {t}" for t in tried), file=sys.stderr)
    sys.exit(1)


def load_group(items):
    out = []
    for x in items:
        sc = x.get("score")
        if isinstance(sc, (int, float)):
            out.append((x.get("id", "?"), float(sc)))
    return out


def jittered_y(base_y, n, spread=0.15, seed=1337):
    rng = random.Random(seed)
    return [base_y + rng.uniform(-spread, spread) for _ in range(n)]


def plot_points(ax, xs, ys, label):
    ax.scatter(xs, ys, s=38, alpha=0.8, label=label, edgecolor="black")


def annotate_top(ax, pts, top_k, align="left"):
    if not pts or top_k <= 0:
        return
    sorted_by = sorted(pts, key=lambda t: t[1])
    worst = sorted_by[:min(top_k, len(sorted_by))]
    best = sorted_by[-min(top_k, len(sorted_by)):]
    for pid, sc, x, y in worst + best:
        ax.annotate(f"{pid} ({sc:.1f})", (x, y),
                    xytext=(5 if align == "left" else -5, 5),
                    textcoords="offset points", fontsize=8)


def overlay_hist(c_scores, r_scores, out_png):
    if not c_scores and not r_scores:
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    all_vals = (c_scores or []) + (r_scores or [])
    bins = max(8, min(40, int(math.sqrt(len(all_vals)) * 2)))
    ax.hist(c_scores, bins=bins, density=True, alpha=0.55, label="Dataset C", edgecolor="black")
    ax.hist(r_scores, bins=bins, density=True, alpha=0.55, label="Random", edgecolor="black")
    if c_scores:
        mc = mean(c_scores)
        ax.axvline(mc, linestyle="--", linewidth=2, label=f"Average Dataset C = {mc:.2f}")
    if r_scores:
        mr = mean(r_scores)
        ax.axvline(mr, linestyle="-.", linewidth=2, label=f"Average Random = {mr:.2f}")
    ax.set_title("Q4 — Dataset C vs Random (score distribution)")
    ax.set_xlabel("Alignment score")
    ax.set_ylabel("Density")
    ax.grid(True, linewidth=0.35, alpha=0.4)
    ax.legend()
    fig.tight_layout(pad=1.0)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main():
    ensure_dirs()
    VIZ_Q4.mkdir(parents=True, exist_ok=True)

    p = argparse.ArgumentParser(description="Q4: compare Dataset C vs Random (overlap + per-sequence).")
    p.add_argument("results_json", nargs="?", default=str(Q4_RESULTS_PATH),
                   help="Path to q4_results.json (default: results/q4_results.json)")
    p.add_argument("--outdir", default=str(VIZ_Q4),
                   help="Output directory (default: results/visualization/Q4)")
    p.add_argument("--label-top", type=int, default=0,
                   help="Annotate top N best and worst points per group")
    args = p.parse_args()

    in_json = Path(args.results_json)
    if in_json == Q4_RESULTS_PATH or not in_json.is_file():
        in_json = resolve_results(str(in_json))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with in_json.open("r", encoding="utf-8") as f:
        data = json.load(f)

    c_group = load_group(data.get("datasetC", []))
    r_group = load_group(data.get("random", []))
    cx = [sc for (_id, sc) in c_group]
    rx = [sc for (_id, sc) in r_group]

    overlay_hist(cx, rx, outdir / "overlap_hist.png")

    fig, ax = plt.subplots(figsize=(11, 6))
    cy = jittered_y(1.0, len(cx), spread=0.12, seed=112)
    ry = jittered_y(0.0, len(rx), spread=0.12, seed=221)
    plot_points(ax, cx, cy, "Dataset C")
    plot_points(ax, rx, ry, "Random")
    if cx:
        mc = mean(cx)
        ax.axvline(mc, ymin=0.5, ymax=1.0, linestyle="--", linewidth=2, label=f"Average Dataset C = {mc:.2f}")
    if rx:
        mr = mean(rx)
        ax.axvline(mr, ymin=0.0, ymax=0.5, linestyle="-.", linewidth=2, label=f"Average Random = {mr:.2f}")
    if args.label_top > 0:
        c_pts = [(pid, sc, sc, y) for (pid, sc), y in zip(c_group, cy)]
        r_pts = [(pid, sc, sc, y) for (pid, sc), y in zip(r_group, ry)]
        annotate_top(ax, c_pts, args.label_top, align="left")
        annotate_top(ax, r_pts, args.label_top, align="right")
    ax.set_title("Q4 — Where each sequence landed (score per sequence)")
    ax.set_xlabel("Alignment score")
    ax.set_yticks([0.0, 1.0], labels=["random", "datasetC"])
    ax.grid(True, axis="x", linewidth=0.35, alpha=0.4)
    ax.legend(loc="best")
    fig.tight_layout(pad=1.0)
    fig.savefig(outdir / "where_each_landed.png", dpi=300)
    plt.close(fig)

    print(f"Read results from: {in_json}")
    print(f"Saved: {outdir}")
    if cx and rx:
        print(f"Average Dataset C = {mean(cx):.2f} | Average Random = {mean(rx):.2f} | Δ = {(mean(cx) - mean(rx)):.2f}")


if __name__ == "__main__":
    main()
