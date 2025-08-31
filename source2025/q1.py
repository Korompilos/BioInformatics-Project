import os
import random

from config import DATASET_A_PATH, DATASET_B_PATH, DATASET_C_PATH, ensure_dirs

ALPHABET = ["A", "C", "G", "T"]
PATTERNS = ["ATTAGA", "ACGCATTT", "AGGACTCAA", "ATTTCAGT"]


def random_symbols(min_n, max_n):
    return "".join(random.choices(ALPHABET, k=random.randint(min_n, max_n)))


def mutate_pattern(pattern):
    p = list(pattern)
    for pos in random.sample(range(len(p)), k=random.randint(0, 2)):
        if random.choice([True, False]):
            choices = [c for c in ALPHABET if c != p[pos]]
            p[pos] = random.choice(choices)
        else:
            p[pos] = ""
    return "".join(p)


def synthesize_string():
    s = random_symbols(1, 3)
    for pat in PATTERNS:
        s += mutate_pattern(pat)
    s += random_symbols(1, 2)
    return s


def write_list(lines, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for seq in lines:
            f.write(seq + "\n")
    print(f"Saved {len(lines)} sequences -> {path}")


def main():
    ensure_dirs()
    all_strings = [synthesize_string() for _ in range(100)]
    random.shuffle(all_strings)
    datasetA = all_strings[:10]
    datasetB = all_strings[10:80]
    datasetC = all_strings[80:]
    write_list(datasetA, DATASET_A_PATH)
    write_list(datasetB, DATASET_B_PATH)
    write_list(datasetC, DATASET_C_PATH)


if __name__ == "__main__":
    main()
