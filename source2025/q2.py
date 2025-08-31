from pathlib import Path

A = 2
MS = 1
MM = -1
GP = -2

try:
    from config import DATASET_A_PATH, Q2_ALIGNMENT_PATH, ensure_dirs
except Exception:
    ROOT = Path(__file__).resolve().parents[1]
    AUX_DIR = ROOT / "auxiliary2025"
    DATASET_A_PATH = AUX_DIR / "datasetA.txt"
    Q2_ALIGNMENT_PATH = AUX_DIR / "q2_alignment.txt"


    def ensure_dirs():
        AUX_DIR.mkdir(parents=True, exist_ok=True)


def read_sequences(p: Path):
    with p.open('r', encoding='utf-8') as f:
        return [l.strip().upper() for l in f if l.strip()]


def nw(a, b):
    m, n = len(a), len(b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    bk = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i * GP
        bk[i][0] = 1
    for j in range(1, n + 1):
        dp[0][j] = j * GP
        bk[0][j] = 2
    for i in range(1, m + 1):
        ai = a[i - 1]
        for j in range(1, n + 1):
            bj = b[j - 1]
            d = dp[i - 1][j - 1] + (MS if ai == bj else MM)
            u = dp[i - 1][j] + GP
            l = dp[i][j - 1] + GP
            if d >= u and d >= l:
                dp[i][j] = d;
                bk[i][j] = 0
            elif u >= l:
                dp[i][j] = u;
                bk[i][j] = 1
            else:
                dp[i][j] = l;
                bk[i][j] = 2
    i, j = m, n
    ra, rb = [], []
    while i > 0 or j > 0:
        t = bk[i][j]
        if t == 0:
            ra.append(a[i - 1]);
            rb.append(b[j - 1]);
            i -= 1;
            j -= 1
        elif t == 1:
            ra.append(a[i - 1]);
            rb.append('-');
            i -= 1
        else:
            ra.append('-');
            rb.append(b[j - 1]);
            j -= 1
    return ''.join(reversed(ra)), ''.join(reversed(rb)), dp[m][n]


def col_counts(block, col):
    c = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
    for s in block:
        c[s[col]] += 1
    return c


def profile_align(p1, p2):
    L1 = len(p1[0]);
    L2 = len(p2[0])
    cnt1 = [col_counts(p1, j) for j in range(L1)]
    cnt2 = [col_counts(p2, j) for j in range(L2)]

    def sp(c1, c2):
        s = 0
        for x in 'ACGT':
            for y in 'ACGT':
                s += (MS if x == y else MM) * c1[x] * c2[y]
        return s

    def gap_col(c):
        return (c['A'] + c['C'] + c['G'] + c['T']) * GP

    dp = [[0] * (L2 + 1) for _ in range(L1 + 1)]
    bk = [[0] * (L2 + 1) for _ in range(L1 + 1)]
    for i in range(1, L1 + 1):
        dp[i][0] = dp[i - 1][0] + gap_col(cnt1[i - 1]);
        bk[i][0] = 1
    for j in range(1, L2 + 1):
        dp[0][j] = dp[0][j - 1] + gap_col(cnt2[j - 1]);
        bk[0][j] = 2
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            d = dp[i - 1][j - 1] + sp(cnt1[i - 1], cnt2[j - 1])
            u = dp[i - 1][j] + gap_col(cnt1[i - 1])
            l = dp[i][j - 1] + gap_col(cnt2[j - 1])
            if d >= u and d >= l:
                dp[i][j] = d;
                bk[i][j] = 0
            elif u >= l:
                dp[i][j] = u;
                bk[i][j] = 1
            else:
                dp[i][j] = l;
                bk[i][j] = 2
    i, j = L1, L2
    map1, map2 = [], []
    while i > 0 or j > 0:
        t = bk[i][j]
        if t == 0:
            map1.append(1);
            map2.append(1);
            i -= 1;
            j -= 1
        elif t == 1:
            map1.append(1);
            map2.append(0);
            i -= 1
        else:
            map1.append(0);
            map2.append(1);
            j -= 1
    map1.reverse();
    map2.reverse()
    out = []
    for row in p1:
        r = [];
        k = 0
        for t in map1:
            r.append(row[k] if t else '-')
            if t: k += 1
        out.append(''.join(r))
    for row in p2:
        r = [];
        k = 0
        for t in map2:
            r.append(row[k] if t else '-')
            if t: k += 1
        out.append(''.join(r))
    return out


def guide_tree_blocks(seqs):
    n = len(seqs)
    scores = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            _, _, sc = nw(seqs[i], seqs[j])
            scores[i][j] = sc;
            scores[j][i] = sc
    clusters = {i: ({'block': [seqs[i]], 'members': [i]}) for i in range(n)}
    active = set(range(n))
    next_id = n

    def sim(mem_a, mem_b):
        return max(scores[i][j] for i in mem_a for j in mem_b)

    while len(active) > 1:
        a, b = max(((x, y) for x in active for y in active if x < y),
                   key=lambda ab: sim(clusters[ab[0]]['members'], clusters[ab[1]]['members']))
        merged = profile_align(clusters[a]['block'], clusters[b]['block'])
        new_id = next_id;
        next_id += 1
        clusters[new_id] = {'block': merged, 'members': clusters[a]['members'] + clusters[b]['members']}
        active.remove(a);
        active.remove(b);
        active.add(new_id)
        del clusters[a];
        del clusters[b]
    return clusters[next(iter(active))]['block']


def main():
    ensure_dirs()
    seqs = read_sequences(Path(DATASET_A_PATH))
    msa = guide_tree_blocks(seqs)
    out_path = Path(Q2_ALIGNMENT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', encoding='utf-8') as f:
        for r in msa:
            f.write(r + '\n')
    for r in msa:
        print(r)
    print(f"[Q2] Wrote: {out_path}")


if __name__ == '__main__':
    main()
