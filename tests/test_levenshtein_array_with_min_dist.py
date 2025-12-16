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
N = 100
K = 8

kmers = [random_kmer(K) for _ in range(N)]
target = random_kmer(K)

# -----------------------------
# Test 1: levenshtein_distance_array_with_min_dist - all distances above min_distance
# -----------------------------
def test_levenshtein_array_min_dist_true():
    logging.info("Test 1: levenshtein_distance_array_with_min_dist - all distances ≥ min_distance")

    min_dist = 0  # guarantees no early exit
    flag, distances = quickmers.levenshtein_distance_array_with_min_dist(target, kmers, min_dist)
    distances_expected = [levenshtein_naive(target, kmer) for kmer in kmers]

    assert flag is True, "Flag should be True when all distances ≥ min_distance"
    assert distances == distances_expected, "Distances do not match naive calculation"

    logging.info("✅ Test 1 passed")

# -----------------------------
# Test 2: levenshtein_distance_array_with_min_dist - early exit
# -----------------------------
def test_levenshtein_array_min_dist_false():
    logging.info("Test 2: levenshtein_distance_array_with_min_dist - early exit")

    # Construct kmers so first kmer is below min_distance
    target_local = "AAAA"
    kmers_local = ["AAAC", "CCCC", "GGGG", "TTTT"]  # first distance = 1
    min_dist = 2  # triggers early exit

    flag, distances = quickmers.levenshtein_distance_array_with_min_dist(target_local, kmers_local, min_dist)
    distances_expected = [levenshtein_naive(target_local, kmer) for kmer in kmers_local]

    assert flag is False, "Flag should be False due to early exit"
    assert distances[0] == distances_expected[0], "First distance should match naive"
    assert all(d is None for d in distances[1:]), "All distances after early exit should be None"

    logging.info("✅ Test 2 passed")

# -----------------------------
# Test 3: levenshtein_distance_array_with_min_dist - empty input
# -----------------------------
def test_levenshtein_array_min_dist_empty():
    logging.info("Test 3: levenshtein_distance_array_with_min_dist - empty input list")

    flag, distances = quickmers.levenshtein_distance_array_with_min_dist(target, [], 1)
    assert flag is True, "Flag should be True for empty input"
    assert distances == [], "Distances should be empty list"

    logging.info("✅ Test 3 passed")
