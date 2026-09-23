import quickmers
import logging
import sys
import time
from collections import Counter
from itertools import product

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
# Constants
# -----------------------------
ALPHABET = ["A", "C", "G", "T"]

# -----------------------------
# Naive Levenshtein
# -----------------------------
def levenshtein(a: str, b: str) -> int:
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost
            )

    return dp[n][m]


# -----------------------------
# Brute-force ball (ground truth)
# -----------------------------
def brute_force_edit_ball(kmer: str, radius: int) -> set:
    k = len(kmer)
    ball = set()

    for tup in product(ALPHABET, repeat=k):
        s = "".join(tup)
        if levenshtein(kmer, s) <= radius:
            ball.add(s)

    return ball


# -----------------------------
# Test parameters
# -----------------------------
# You can tune these; keep k small when radius is nontrivial.
K = 6
RADIUS = 2
TEST_KMER = "ACGTAC"  # len=6; adjust if you change K

# -----------------------------
# Test 1: fixed_length_levenshtein_ball (list version)
# -----------------------------
def test_fixed_length_levenshtein_ball_correctness():
    logging.info("Test 1: fixed_length_levenshtein_ball (list version)")

    start_py = time.time()
    ball_py = brute_force_edit_ball(TEST_KMER, RADIUS)
    t_py = time.time() - start_py
    logging.info(f"Python brute-force ball size: {len(ball_py)}, time: {t_py:.4f}s")

    start_c = time.time()
    ball_c = set(quickmers.fixed_length_levenshtein_ball(TEST_KMER, RADIUS, 0))
    t_c = time.time() - start_c
    logging.info(f"C list ball size: {len(ball_c)}, time: {t_c:.6f}s")

    missing = ball_py - ball_c
    extra = ball_c - ball_py

    assert len(missing) == 0 and len(extra) == 0, (
        f"Ball mismatch for list version:\n"
        f"  missing: {len(missing)}\n"
        f"  extra:   {len(extra)}\n"
        f"  sample missing: {list(missing)[:5]}\n"
        f"  sample extra:   {list(extra)[:5]}"
    )

    logging.info("✅ fixed_length_levenshtein_ball (list) test passed")


# -----------------------------
# Test 2: fixed_length_levenshtein_ball_iterator
# -----------------------------
def test_fixed_length_levenshtein_ball_iterator_correctness():
    logging.info("Test 2: fixed_length_levenshtein_ball_iterator")

    start_py = time.time()
    ball_py = brute_force_edit_ball(TEST_KMER, RADIUS)
    t_py = time.time() - start_py
    logging.info(f"Python brute-force ball size: {len(ball_py)}, time: {t_py:.4f}s")

    start_c_iter = time.time()
    it = quickmers.fixed_length_levenshtein_ball_iterator(TEST_KMER, RADIUS, 0)
    ball_c_iter = set()
    for s in it:
        ball_c_iter.add(s)
    t_c_iter = time.time() - start_c_iter
    logging.info(f"C iterator ball size: {len(ball_c_iter)}, time: {t_c_iter:.6f}s")

    missing = ball_py - ball_c_iter
    extra = ball_c_iter - ball_py

    assert len(missing) == 0 and len(extra) == 0, (
        f"Ball mismatch for iterator version:\n"
        f"  missing: {len(missing)}\n"
        f"  extra:   {len(extra)}\n"
        f"  sample missing: {list(missing)[:5]}\n"
        f"  sample extra:   {list(extra)[:5]}"
    )

    logging.info("✅ fixed_length_levenshtein_ball_iterator test passed")


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    logging.info(f"Testing fixed-length Levenshtein ball for kmer={TEST_KMER}, radius={RADIUS}")

    test_fixed_length_levenshtein_ball_correctness()
    test_fixed_length_levenshtein_ball_iterator_correctness()

    logging.info("All fixed-length Levenshtein ball correctness tests passed.")