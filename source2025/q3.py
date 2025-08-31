import json
from collections import defaultdict, Counter

from config import DATASET_A_PATH, DATASET_B_PATH, JSON_DIR, Q2_ALIGNMENT_PATH, ensure_dirs

ALPHABET = ["A", "C", "G", "T"]
ALPHA = 2
GAP_THRESHOLD = 0.5
LAPLACE = 1e-3


def read_sequences(filename):
    with open(filename, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def pairwise_align(s1, s2, alpha=2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    back = [[(0, 0)] * (n + 1) for _ in range(m + 1)]
    a = alpha
    for i in range(1, m + 1):
        dp[i][0] = -a * i
        back[i][0] = (i - 1, 0)
    for j in range(1, n + 1):
        dp[0][j] = -a * j
        back[0][j] = (0, j - 1)
    half = -a / 2
    for i in range(1, m + 1):
        si = s1[i - 1]
        for j in range(1, n + 1):
            sc_match = dp[i - 1][j - 1] + (1 if si == s2[j - 1] else half)
            sc_del = dp[i - 1][j] - a
            sc_ins = dp[i][j - 1] - a
            if sc_match >= sc_del and sc_match >= sc_ins:
                dp[i][j] = sc_match
                back[i][j] = (i - 1, j - 1)
            elif sc_del >= sc_ins:
                dp[i][j] = sc_del
                back[i][j] = (i - 1, j)
            else:
                dp[i][j] = sc_ins
                back[i][j] = (i, j - 1)
    a1, a2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        pi, pj = back[i][j]
        if pi == i - 1 and pj == j - 1:
            a1.append(s1[i - 1]);
            a2.append(s2[j - 1])
        elif pi == i - 1:
            a1.append(s1[i - 1]);
            a2.append("-")
        else:
            a1.append("-");
            a2.append(s2[j - 1])
        i, j = pi, pj
    a1.reverse();
    a2.reverse()
    return "".join(a1), "".join(a2)


class ProfileHMM:
    def __init__(self):
        self.match_cols = []
        self.M_states = []
        self.I_states = []
        self.D_states = []
        self.emissions = defaultdict(lambda: Counter())
        self.transitions = defaultdict(lambda: Counter())

    def normalize_probs(self):
        trans_probs = {}
        for st, cnts in self.transitions.items():
            total = float(sum(cnts.values()))
            if total:
                trans_probs[st] = {k: v / total for k, v in cnts.items()}
        emit_probs = {}
        for st, cnts in self.emissions.items():
            total = float(sum(cnts[a] for a in ALPHABET))
            if total:
                emit_probs[st] = {a: cnts[a] / total for a in ALPHABET}
        return emit_probs, trans_probs


def build_profile_hmm_from_alignment(aln, gap_threshold=GAP_THRESHOLD):
    hmm = ProfileHMM()
    nseq = len(aln)
    L = len(aln[0])
    match_cols = []
    for c in range(L):
        if sum(1 for s in aln if s[c] == "-") / nseq <= gap_threshold:
            match_cols.append(c)
    hmm.match_cols = match_cols
    K = len(match_cols)
    hmm.M_states = [f"M{i}" for i in range(1, K + 1)]
    hmm.D_states = [f"D{i}" for i in range(1, K + 1)]
    hmm.I_states = [f"I{i}" for i in range(0, K + 1)]
    for s in hmm.M_states + hmm.I_states:
        e = hmm.emissions[s]
        lap = LAPLACE
        e["A"] += lap;
        e["C"] += lap;
        e["G"] += lap;
        e["T"] += lap
    legal_next = defaultdict(list)
    for i in range(0, K + 1):
        ii = f"I{i}"
        if i == 0:
            legal_next["Start"] += [ii, "D1" if K >= 1 else "End", "M1" if K >= 1 else "End"]
        if i < K:
            Mi = f"M{i + 1}"
            Di = f"D{i + 1}"
            legal_next[Mi] += [f"M{i + 2}" if i + 1 < K else "End", f"I{i + 1}", f"D{i + 2}" if i + 1 < K else "End"]
            legal_next[Di] += [f"D{i + 2}" if i + 1 < K else "End", f"M{i + 2}" if i + 1 < K else "End"]
        legal_next[ii] += [ii] + ([f"M{i + 1}"] if i < K else ["End"])
    for st, nxts in legal_next.items():
        trow = hmm.transitions[st]
        lap = LAPLACE
        for nx in nxts:
            trow[nx] += lap
    mc = match_cols
    for row in aln:
        path = ["Start"]
        if K > 0:
            lead = row[:mc[0]]
            for ch in lead:
                if ch != "-":
                    path.append("I0");
                    hmm.emissions["I0"][ch] += 1
        else:
            for ch in row:
                if ch != "-":
                    path.append("I0");
                    hmm.emissions["I0"][ch] += 1
        for mi, col in enumerate(mc, start=1):
            ch = row[col]
            if ch == "-":
                path.append(f"D{mi}")
            else:
                path.append(f"M{mi}");
                hmm.emissions[f"M{mi}"][ch] += 1
            next_col = mc[mi] if mi < K else None
            ins_state = f"I{mi}"
            seg = row[col + 1:next_col] if next_col is not None else row[col + 1:]
            for mid_ch in seg:
                if mid_ch != "-":
                    path.append(ins_state);
                    hmm.emissions[ins_state][mid_ch] += 1
        path.append("End")
        tr = hmm.transitions
        for a, b in zip(path, path[1:]):
            tr[a][b] += 1
    return hmm


def consensus_from_profile(hmm, emit_probs):
    out = []
    for i in range(1, len(hmm.M_states) + 1):
        st = f"M{i}"
        if st in emit_probs:
            ep = emit_probs[st]
            out.append(max(ALPHABET, key=ep.get))
        else:
            out.append("A")
    return "".join(out)


def update_hmm_with_datasetB(hmm, emit_probs, datasetB):
    consensus = consensus_from_profile(hmm, emit_probs)
    K = len(hmm.M_states)
    em = hmm.emissions
    tr = hmm.transitions
    for seq in datasetB:
        a_cons, a_seq = pairwise_align(consensus, seq, ALPHA)
        path = ["Start"]
        i_match = 0
        for cc, ss in zip(a_cons, a_seq):
            if cc != "-" and ss != "-":
                i_match += 1
                st = f"M{i_match}"
                em[st][ss] += 1
                path.append(st)
            elif cc != "-" and ss == "-":
                i_match += 1
                path.append(f"D{i_match}")
            elif cc == "-" and ss != "-":
                ins_idx = i_match if i_match <= K else K
                st = f"I{ins_idx}"
                em[st][ss] += 1
                path.append(st)
        path.append("End")
        for a, b in zip(path, path[1:]):
            tr[a][b] += 1


def pretty_print_transitions(trans_probs):
    print("\n# === Transition Probabilities ===")
    for st in sorted(trans_probs.keys(), key=lambda s: (s != "Start", s != "End", s)):
        row = {k: round(v, 3) for k, v in sorted(trans_probs[st].items())}
        print(f"{st:>6} -> {row}")


def pretty_print_emissions(emit_probs):
    print("\n# === Emission Probabilities (A,C,G,T) ===")
    for st in sorted(emit_probs.keys(), key=lambda s: (s[0] != "M", s[0] != "I", s)):
        ep = emit_probs[st]
        row = [round(ep.get(a, 0.0), 3) for a in ALPHABET]
        print(f"{st:>6} : {{'A': {row[0]}, 'C': {row[1]}, 'G': {row[2]}, 'T': {row[3]}}}")


if __name__ == "__main__":
    ensure_dirs()
    align_path = Q2_ALIGNMENT_PATH if Q2_ALIGNMENT_PATH.exists() else DATASET_A_PATH
    alignment = read_sequences(str(align_path))
    L = len(alignment[0])
    assert all(len(s) == L for s in alignment)
    hmm = build_profile_hmm_from_alignment(alignment, GAP_THRESHOLD)
    emit_probs, trans_probs = hmm.normalize_probs()
    pretty_print_transitions(trans_probs)
    pretty_print_emissions(emit_probs)
    out_before = JSON_DIR / "hmm_profile_before.json"
    with open(out_before, "w", encoding="utf-8") as f:
        json.dump({"emissions": emit_probs, "transitions": trans_probs}, f, indent=2)
    datasetB = read_sequences(str(DATASET_B_PATH))
    update_hmm_with_datasetB(hmm, emit_probs, datasetB)
    emit_probs, trans_probs = hmm.normalize_probs()
    pretty_print_transitions(trans_probs)
    pretty_print_emissions(emit_probs)
    out_after = JSON_DIR / "hmm_profile_after.json"
    with open(out_after, "w", encoding="utf-8") as f:
        json.dump({"emissions": emit_probs, "transitions": trans_probs}, f, indent=2)
