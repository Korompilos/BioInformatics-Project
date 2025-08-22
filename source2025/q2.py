def read_sequences(filename):
    with open(filename, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]

alpha = 2  # gap penalty

def pairwise_align(s1, s2, alpha=2):
    m, n = len(s1), len(s2)
    dp = [[0]*(n+1) for _ in range(m+1)]
    back = [[None]*(n+1) for _ in range(m+1)]

    for i in range(1, m+1):
        dp[i][0] = -alpha * i
        back[i][0] = (i-1, 0)
    for j in range(1, n+1):
        dp[0][j] = -alpha * j
        back[0][j] = (0, j-1)

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + (1 if s1[i-1] == s2[j-1] else -alpha/2)
            delete = dp[i-1][j] - alpha
            insert = dp[i][j-1] - alpha
            max_score = max(match, delete, insert)
            dp[i][j] = max_score
            if max_score == match:
                back[i][j] = (i-1, j-1)
            elif max_score == delete:
                back[i][j] = (i-1, j)
            else:
                back[i][j] = (i, j-1)

    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        prev = back[i][j]
        if prev is None:
            break
        if prev == (i-1, j-1):
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i, j = i-1, j-1
        elif prev == (i-1, j):
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            j -= 1
    return align1, align2

def multiple_alignment(seqs, alpha=2):
    center = seqs[0]
    alignments = [center]
    for seq in seqs[1:]:
        aligned_center, aligned_seq = pairwise_align(alignments[0], seq, alpha)
        new_alignments = []
        for s in alignments:
            new_s = ""
            k = 0
            for c in aligned_center:
                if c == "-":
                    new_s += "-"
                else:
                    new_s += s[k]
                    k += 1
            new_alignments.append(new_s)
        new_alignments.append(aligned_seq)
        alignments = new_alignments
        alignments[0] = aligned_center
    return alignments

if __name__ == "__main__":
    sequences = read_sequences("auxiliary2025\\datasetA.txt")
    result = multiple_alignment(sequences, alpha)
    for row in result:
        print(row)