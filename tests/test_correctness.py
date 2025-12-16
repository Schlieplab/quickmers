import quickmers
import random
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
    return "".join(random.choice(alphabet) for _ in range(k))

def hamming_naive(s1: str, s2: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein_naive(s1: str, s2: str) -> int:
    m, n = len(s1), len(s2)
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1): dp[i][0] = i
    for j in range(n+1): dp[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            cost = 0 if s1[i-1]==s2[j-1] else 1
            dp[i][j] = min(dp[i-1][j]+1, dp[i][j-1]+1, dp[i-1][j-1]+cost)
    return dp[m][n]

# -----------------------------
# Test parameters
# -----------------------------
N = 1000
K = 16

kmers = [random_kmer(K) for _ in range(N)]
target = random_kmer(K)

# -----------------------------
# Test 1: hamming_distance
# -----------------------------
def test_hamming_distance_correctness():
    logging.info("Test 1: hamming_distance (single string vs string)")

    for kmer in kmers[:10]:  # check first 10 for quick logging
        dist_py = hamming_naive(target, kmer)
        dist_c = quickmers.hamming_distance(target, kmer)
        assert dist_py == dist_c, f"Hamming mismatch: {dist_py} vs {dist_c}"

    logging.info("✅ hamming_distance test passed")

# -----------------------------
# Test 2: hamming_distance_array
# -----------------------------
def test_hamming_distance_array_correctness():
    logging.info("Test 2: hamming_distance_array (single string vs list)")

    distances_py = [hamming_naive(target, kmer) for kmer in kmers]
    distances_c = quickmers.hamming_distance_array(target, kmers)
    assert np.all(distances_py == distances_c), "Hamming array mismatch"
    logging.info("✅ hamming_distance_array test passed")

# -----------------------------
# Test 3: levenshtein (scalar)
# -----------------------------
def test_levenshtein_correctness():
    logging.info("Test 3: levenshtein (single string vs string)")

    for kmer in kmers[:10]:  # check first 10 for quick logging
        dist_py = levenshtein_naive(target, kmer)
        dist_c = quickmers.levenshtein_distance(target, kmer)
        assert dist_py == dist_c, f"Levenshtein mismatch: {dist_py} vs {dist_c}"

    logging.info("✅ levenshtein scalar test passed")

# -----------------------------
# Test 4: levenshtein_array
# -----------------------------
def test_levenshtein_array_correctness():
    logging.info("Test 4: levenshtein_array (single string vs list)")

    distances_py = [levenshtein_naive(target, kmer) for kmer in kmers]
    distances_c = quickmers.levenshtein_distance_array(target, kmers)
    assert np.all(distances_py == distances_c), "Levenshtein array mismatch"

    logging.info("✅ levenshtein_array test passed")
