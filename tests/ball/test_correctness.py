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


def brute_force_variable_edit_ball(kmer: str, radius: int) -> set:
    """All strings of any length (k - radius .. k + radius) within radius."""
    k = len(kmer)
    ball = set()

    for length in range(max(0, k - radius), k + radius + 1):
        for tup in product(ALPHABET, repeat=length):
            s = "".join(tup)
            if levenshtein(kmer, s) <= radius:
                ball.add(s)

    return ball


BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode_kmer(kmer: str) -> int:
    val = 0
    for c in kmer:
        val = (val << 2) | BASE_TO_BITS[c]
    return val


def assert_same_ball(name, ball_py, result):
    """Compare a C result (list) to the reference set, including duplicates."""
    ball_c = set(result)
    missing = ball_py - ball_c
    extra = ball_c - ball_py

    assert len(missing) == 0 and len(extra) == 0, (
        f"Ball mismatch for {name}:\n"
        f"  missing: {len(missing)}\n"
        f"  extra:   {len(extra)}\n"
        f"  sample missing: {list(missing)[:5]}\n"
        f"  sample extra:   {list(extra)[:5]}"
    )
    assert len(result) == len(ball_c), (
        f"{name} returned {len(result) - len(ball_c)} duplicates"
    )


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
# Test 3: fixed-length ball for small radii and edge cases
# -----------------------------
def test_fixed_length_levenshtein_ball_small_radii():
    logging.info("Test 3: fixed_length_levenshtein_ball small radii")

    for kmer, radius in [("ACG", 0), ("ACG", 1), ("ACGT", 1), ("ACGTA", 3), ("A", 1)]:
        ball_py = brute_force_edit_ball(kmer, radius)
        assert_same_ball(f"fixed list kmer={kmer} r={radius}", ball_py,
                         quickmers.fixed_length_levenshtein_ball(kmer, radius, 0))
        assert_same_ball(f"fixed iterator kmer={kmer} r={radius}", ball_py,
                         list(quickmers.fixed_length_levenshtein_ball_iterator(kmer, radius, 0)))

    logging.info("✅ fixed-length small radii test passed")


# -----------------------------
# Test 4: levenshtein_ball (variable length, list version)
# -----------------------------
def test_levenshtein_ball_correctness():
    logging.info("Test 4: levenshtein_ball (variable length, list version)")

    start_py = time.time()
    ball_py = brute_force_variable_edit_ball(TEST_KMER, RADIUS)
    t_py = time.time() - start_py
    logging.info(f"Python brute-force ball size: {len(ball_py)}, time: {t_py:.4f}s")

    start_c = time.time()
    result = quickmers.levenshtein_ball(TEST_KMER, RADIUS)
    t_c = time.time() - start_c
    logging.info(f"C list ball size: {len(result)}, time: {t_c:.6f}s")

    assert_same_ball("levenshtein_ball (list)", ball_py, result)

    logging.info("✅ levenshtein_ball (list) test passed")


# -----------------------------
# Test 5: levenshtein_ball_iterator (variable length)
# -----------------------------
def test_levenshtein_ball_iterator_correctness():
    logging.info("Test 5: levenshtein_ball_iterator (variable length)")

    start_py = time.time()
    ball_py = brute_force_variable_edit_ball(TEST_KMER, RADIUS)
    t_py = time.time() - start_py
    logging.info(f"Python brute-force ball size: {len(ball_py)}, time: {t_py:.4f}s")

    start_c_iter = time.time()
    result = [s for s in quickmers.levenshtein_ball_iterator(TEST_KMER, RADIUS)]
    t_c_iter = time.time() - start_c_iter
    logging.info(f"C iterator ball size: {len(result)}, time: {t_c_iter:.6f}s")

    assert_same_ball("levenshtein_ball_iterator", ball_py, result)

    logging.info("✅ levenshtein_ball_iterator test passed")


# -----------------------------
# Test 6: variable-length ball for small radii and edge cases
# -----------------------------
def test_levenshtein_ball_small_radii():
    logging.info("Test 6: levenshtein_ball small radii and edge cases")

    # ("A", 1) and ("", 2) include the empty string
    for kmer, radius in [("ACG", 0), ("ACG", 1), ("ACGTA", 3), ("A", 1), ("", 2)]:
        ball_py = brute_force_variable_edit_ball(kmer, radius)
        assert_same_ball(f"list kmer={kmer!r} r={radius}", ball_py,
                         quickmers.levenshtein_ball(kmer, radius))
        assert_same_ball(f"iterator kmer={kmer!r} r={radius}", ball_py,
                         list(quickmers.levenshtein_ball_iterator(kmer, radius)))

    logging.info("✅ levenshtein_ball small radii test passed")


# -----------------------------
# Test 7: variable-length ball binary output
# -----------------------------
def test_levenshtein_ball_binary():
    logging.info("Test 7: levenshtein_ball binary output (encoded, length)")

    ball_py = brute_force_variable_edit_ball(TEST_KMER, RADIUS)
    bin_py = {(encode_kmer(s), len(s)) for s in ball_py}

    assert_same_ball("levenshtein_ball binary (list)", bin_py,
                     quickmers.levenshtein_ball(TEST_KMER, RADIUS, 1))
    assert_same_ball("levenshtein_ball binary (iterator)", bin_py,
                     list(quickmers.levenshtein_ball_iterator(TEST_KMER, RADIUS, 1)))

    logging.info("✅ levenshtein_ball binary test passed")


# -----------------------------
# Test 8: invalid arguments
# -----------------------------
def test_levenshtein_ball_invalid_args():
    logging.info("Test 8: invalid arguments raise ValueError")

    calls = [
        lambda: quickmers.levenshtein_ball("A" * 30, 2),          # k + radius > 31
        lambda: quickmers.levenshtein_ball_iterator("A" * 30, 2),
        lambda: quickmers.levenshtein_ball("ACG", -1),
        lambda: quickmers.fixed_length_levenshtein_ball("A" * 33, 1, 0),  # k > 32
        lambda: quickmers.fixed_length_levenshtein_ball_iterator("ACG", -1),
    ]
    for call in calls:
        try:
            call()
        except ValueError:
            continue
        raise AssertionError("expected ValueError")

    logging.info("✅ invalid arguments test passed")


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    logging.info(f"Testing Levenshtein balls for kmer={TEST_KMER}, radius={RADIUS}")

    test_fixed_length_levenshtein_ball_correctness()
    test_fixed_length_levenshtein_ball_iterator_correctness()
    test_fixed_length_levenshtein_ball_small_radii()
    test_levenshtein_ball_correctness()
    test_levenshtein_ball_iterator_correctness()
    test_levenshtein_ball_small_radii()
    test_levenshtein_ball_binary()
    test_levenshtein_ball_invalid_args()

    logging.info("All Levenshtein ball correctness tests passed.")