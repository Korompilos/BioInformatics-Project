import argparse
import json
import os
import random
import sys
from statistics import mean, median

ALPHABET = ["A", "C", "G", "T"]
ALPHA = 2
nrand = 20
seed = 42


def read_sequences(path):
    with open(path, "r", encoding="utf-8") as f:
        return [l.strip() for l in f if l.strip()]


def pairwise_align_with_ops(s1, s2, alpha=2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    back = [[None] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = -alpha * i
        back[i][0] = ("D", i - 1, 0)
    for j in range(1, n + 1):
        dp[0][j] = -alpha * j
        back[0][j] = ("I", 0, j - 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            sc_match = dp[i - 1][j - 1] + (1 if s1[i - 1] == s2[j - 1] else -alpha / 2)
            sc_del = dp[i - 1][j] - alpha
            sc_ins = dp[i][j - 1] - alpha
            best = max(sc_match, sc_del, sc_ins)
            dp[i][j] = best
            if best == sc_match:
                back[i][j] = ("M" if s1[i - 1] == s2[j - 1] else "S", i - 1, j - 1)
            elif best == sc_del:
                back[i][j] = ("D", i - 1, j)
            else:
                back[i][j] = ("I", i, j - 1)
    ops = []
    i, j = m, n
    while i > 0 or j > 0:
        op, pi, pj = back[i][j]
        ops.append(op)
        i, j = pi, pj
    ops.reverse()
    a1, a2 = [], []
    i = j = 0
    for op in ops:
        if op in ("M", "S"):
            a1.append(s1[i]);
            a2.append(s2[j]);
            i += 1;
            j += 1
        elif op == "D":
            a1.append(s1[i]);
            a2.append("-");
            i += 1
        else:
            a1.append("-");
            a2.append(s2[j]);
            j += 1
    return "".join(a1), "".join(a2), dp[m][n], "".join(ops)


def consensus_from_hmm_emissions(emissions):
    ms = [k for k in emissions.keys() if k.startswith("M") and k[1:].isdigit()]
    ms_sorted = sorted(ms, key=lambda x: int(x[1:]))
    cons = []
    for st in ms_sorted:
        e = emissions.get(st, {})
        best = max(ALPHABET, key=lambda a: e.get(a, 0.0))
        cons.append(best)
    return "".join(cons)


def random_sequences(n, min_len, max_len):
    out = []
    for _ in range(n):
        L = random.randint(min_len, max_len)
        out.append("".join(random.choice(ALPHABET) for _ in range(L)))
    return out


def summarize(scores):
    return {
        "n": len(scores),
        "min": min(scores) if scores else None,
        "max": max(scores) if scores else None,
        "mean": mean(scores) if scores else None,
        "median": median(scores) if scores else None,
    }


def resolve_existing(here, user_path, extra_dirs):
    tried = []
    p = user_path
    if not os.path.isabs(p):
        p = os.path.normpath(os.path.join(here, p))
    tried.append(p)
    if os.path.isfile(p):
        return p, tried
    base = os.path.basename(user_path)
    for d in extra_dirs:
        cand = os.path.normpath(os.path.join(here, d, base))
        tried.append(cand)
        if os.path.isfile(cand):
            return cand, tried
    return None, tried


def main():
    parser = argparse.ArgumentParser(description="Q4: align datasetC and random sequences to HMM consensus.")
    parser.add_argument("hmm_json", nargs="?", default="../results/hmm_profile_after.json",
                        help="Path to HMM JSON (default: ../results/hmm_profile_after.json)")
    parser.add_argument("--datasetC", default=None,
                        help="Path to datasetC.txt (default: ../auxiliary2025/datasetC.txt)")
    parser.add_argument("--out", default="q4_results.json",
                        help="Output JSON file (filename or full path). Default drops under ../results/")
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))

    hmm_path, tried = resolve_existing(
        here,
        args.hmm_json,
        extra_dirs=["..\\results", "../results", "..\\auxiliary2025", "../auxiliary2025"]
    )
    if hmm_path is None:
        msg = ["Could not find the HMM JSON. Tried:"] + [f" - {t}" for t in tried]
        print("\n".join(msg), file=sys.stderr)
        sys.exit(1)

    if args.datasetC:
        c_path, triedC = resolve_existing(
            here,
            args.datasetC,
            extra_dirs=["..\\auxiliary2025", "../auxiliary2025"]
        )
        if c_path is None:
            msg = ["Could not find datasetC. Tried:"] + [f" - {t}" for t in triedC]
            print("\n".join(msg), file=sys.stderr)
            sys.exit(1)
    else:
        c_path = os.path.normpath(os.path.join(here, "..", "auxiliary2025", "datasetC.txt"))

    results_dir = os.path.normpath(os.path.join(here, "..", "results"))
    os.makedirs(results_dir, exist_ok=True)
    out_path = args.out
    if not os.path.isabs(out_path) and os.path.dirname(out_path) == "":
        out_path = os.path.join(results_dir, out_path)

    with open(hmm_path, "r", encoding="utf-8") as f:
        profile = json.load(f)
    emissions = profile["emissions"]

    consensus = consensus_from_hmm_emissions(emissions)
    datasetC = read_sequences(c_path)
    c_lengths = [len(s) for s in datasetC]
    random.seed(seed)
    rnd = random_sequences(nrand, min(c_lengths), max(c_lengths)) if c_lengths else []

    results = {"consensus": consensus, "datasetC": [], "random": [], "summary": {}, "comparison": {}, "comment": ""}

    c_scores = []
    for idx, seq in enumerate(datasetC, 1):
        a1, a2, score, ops = pairwise_align_with_ops(consensus, seq, ALPHA)
        c_scores.append(score)
        results["datasetC"].append({
            "id": f"C{idx}", "sequence": seq, "alignment1": a1, "alignment2": a2, "score": score, "path": ops
        })

    r_scores = []
    for idx, seq in enumerate(rnd, 1):
        _, _, score, ops = pairwise_align_with_ops(consensus, seq, ALPHA)
        r_scores.append(score)
        results["random"].append({"id": f"R{idx}", "sequence": seq, "score": score, "path": ops})

    c_sum = summarize(c_scores);
    r_sum = summarize(r_scores)
    results["summary"]["datasetC"] = c_sum
    results["summary"]["random"] = r_sum

    c_mean = c_sum["mean"] or 0.0
    r_mean = r_sum["mean"] or 0.0
    c_median = c_sum["median"] or 0.0
    r_median = r_sum["median"] or 0.0
    results["comparison"] = {
        "mean_diff": c_mean - r_mean,
        "median_diff": c_median - r_median,
        "mean_ratio": (c_mean / r_mean) if r_mean != 0 else None
    }
    results["comment"] = (
        "datasetC aligns better than random (higher mean score)." if (c_scores and r_scores and c_mean > r_mean)
        else "No stronger mean-score advantage observed for datasetC vs random."
    )

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print("=== Summary ===")
    print("DatasetC:", c_sum)
    print("Random:", r_sum)
    print("Comparison:", results["comparison"])
    print("Comment:", results["comment"])


if __name__ == "__main__":
    main()
