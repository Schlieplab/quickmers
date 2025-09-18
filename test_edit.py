import quickmers
import random
import time
import numpy as np
import logging
import sys
# -----------------------------
# Logging setup
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    stream=sys.stdout,
    format='### INFO - %(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# -----------------------------
# Helper functions
# -----------------------------
def random_kmer(k, alphabet="ACGT"):
    """Generate a random kmer string of length k."""
    return "".join(random.choice(alphabet) for _ in range(k))


C2B = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
B2C = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}

def string2bitstring(s):
    return int('0B'+''.join([C2B[c] for c in s]),2)

def naive_edit_distance(str1: str, str2: str) -> int:
    m, n = len(str1), len(str2)

    # Create a DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i  # deletions
    for j in range(n + 1):
        dp[0][j] = j  # insertions

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                cost = 0
            else:
                cost = 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,      # deletion
                dp[i][j - 1] + 1,      # insertion
                dp[i - 1][j - 1] + cost # substitution
            )

    return dp[m][n]

# -----------------------------
# Benchmark function
# -----------------------------
def benchmark(n, k):
    # Generate random kmers
    kmers_str = [random_kmer(k) for _ in range(n)]
    target_kmer_str = random_kmer(k)

    # Convert kmers to integers
    target_kmer_int = string2bitstring(target_kmer_str)

    # -----------------------------
    # Naive Python benchmark
    # -----------------------------
    start = time.time()
    distances_naive_py = [naive_edit_distance(target_kmer_str, km) for km in kmers_str]
    t_py = time.time() - start
    logging.info(f"Naive Python time: {t_py:.6f} seconds")

    # -----------------------------
    # quickmers C extension benchmark
    # -----------------------------
    start = time.time()
    distances_quickmers = [quickmers.levenshtein(target_kmer_str, km) for km in kmers_str]
    t_np = time.time() - start
    logging.info(f"quickmers levenshtein myers time: {t_np:.6f} seconds")

    # -----------------------------
    # C extension benchmarks
    # -----------------------------
    start = time.time()
    distances_quickmers_list = quickmers.levenshtein_list(target_kmer_str, kmers_str)
    t_np = time.time() - start
    logging.info(f"quickmers levenshtein myers list time: {t_np:.6f} seconds")
    
    # -----------------------------
    # Result check
    # -----------------------------
    if np.all(distances_quickmers == distances_naive_py):
        logging.info("Correct results for naive python calculation ✅")
    else:
        logging.warning("Mismatch for encoded 32bit array results ❌")

    if np.all(distances_quickmers_list == distances_naive_py):
        logging.info("Correct results for naive python calculation ✅")
    else:
        logging.warning("Mismatch for encoded 32bit array results ❌")

    # For extra logging in case of mismatches
    # logging.info(target_kmer_str)
    # logging.info(kmers_str)
    # logging.info(distances_turboshtein)
    # logging.info(distances_naive_py)
    # logging.info(distances_quickmers)
    # logging.info(distances_c)
    # logging.info(distances_np)
    # logging.info(distances_py)

# -----------------------------
# Run benchmark
# -----------------------------
if __name__ == "__main__":
    benchmark(n=1_000, k=32)
    # print("HERE!!" , levenshtein("GACAGGCC", "TCGTCTCA"))
    # print("HERE!!" , quickmers.levenshtein("GGCGCCGGACCGACTGGGGTCAATTGGGCC", "TTCAATGTTATGGTTCGACGTAACGGGCAA"))
    # print("HERE!!" , quickmers.levenshtein("TTCAATGTTATGGTTCGACGTAACGGGCAA", "GGCGCCGGACCGACTGGGGTCAATTGGGCC"))

    # print("HERE!!" , naive_edit_distance("GACAGGCC", "TCGTCTCA"))
    # print("HERE!!" , quickmers.levenshtein("GACAGGCC", "TCGTCTCA"))
