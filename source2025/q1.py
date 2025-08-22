import random

ALPHABET = ['A', 'C', 'G', 'T']
PATTERNS = [
    "ATTAGA",
    "ACGCATTT",
    "AGGACTCAA",
    "ATTTCAGT"
]

def random_symbols(min_n, max_n):
    n = random.randint(min_n, max_n)
    return ''.join(random.choices(ALPHABET, k=n))

def mutate_pattern(pattern):
    pattern = list(pattern)
    positions = random.sample(range(len(pattern)), k=random.randint(0, 2))
    for pos in positions:
        if random.choice([True, False]):
            choices = [c for c in ALPHABET if c != pattern[pos]]
            pattern[pos] = random.choice(choices)
        else:
            pattern[pos] = ''
    return ''.join(pattern)

def synthesize_string():
    s = random_symbols(1, 3)
    for pat in PATTERNS:
        mutated = mutate_pattern(pat)
        s += mutated
    s += random_symbols(1, 2)
    return s

def main():
    all_strings = [synthesize_string() for _ in range(100)]
    random.shuffle(all_strings)

    datasetA = set(all_strings[:10])
    datasetB = set(all_strings[10:80])
    datasetC = set(all_strings[80:])

    print("=== datasetA ({} sequences) ===".format(len(datasetA)))
    for seq in datasetA:
        print(seq)
    print("\n")

    print("=== datasetB ({} sequences) ===".format(len(datasetB)))
    for seq in datasetB:
        print(seq)
    print("\n")

    print("=== datasetC ({} sequences) ===".format(len(datasetC)))
    for seq in datasetC:
        print(seq)
    print("\n")

    # Αποθήκευση datasetA σε αρχείο (για ιι)
    with open("datasetA.txt", "w") as f:
        for seq in datasetA:
            f.write(seq + "\n")

     # Αποθήκευση datasetB σε αρχείο (για ιιι)
    with open("datasetB.txt", "w") as f:
        for seq in datasetB:
            f.write(seq + "\n")

    # Αποθήκευση datasetC σε αρχείο (για ιν)
    with open("datasetC.txt", "w") as f:
        for seq in datasetC:
            f.write(seq + "\n")        


if __name__ == "__main__":
    main()