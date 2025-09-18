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

def popcount32(x):
    """Vectorized population count for uint32 array"""
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return ((x * 0x01010101) >> 24) & 0xFF

def hamming_distance_naive(a, b):
    """Compute Hamming distance between two 2-bit encoded kmers."""
    mismatches = sum(ch1 != ch2 for ch1, ch2 in zip(a, b))
    return mismatches

def hamming_distance_array_numpy(kmer, kmers_array):
    """Vectorized Hamming distance using NumPy, bitwise operations."""
    b = 16 * 2  # assuming 16-mers
    all_lo = (4 ** (b // 2) - 1) // 3
    st_xor = np.bitwise_xor(kmer, kmers_array)
    mask = np.bitwise_and(np.bitwise_or(np.right_shift(st_xor, 1), st_xor), all_lo)
    hamming_dist = popcount32(mask)
    return hamming_dist

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
    kmers_int32 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint32)
    distances_py = [hamming_distance_naive(target_kmer_str, km) for km in kmers_str]
    t_py = time.time() - start
    logging.info(f"Naive Python time: {t_py:.6f} seconds")

    # -----------------------------
    # NumPy vectorized benchmark
    # -----------------------------
    start = time.time()
    kmers_int32 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint32)
    distances_np = hamming_distance_array_numpy(target_kmer_int, kmers_int32)
    t_np = time.time() - start
    logging.info(f"Vectorized NumPy time: {t_np:.6f} seconds")

    # -----------------------------
    # C extension benchmarks
    # -----------------------------

    # Encoded in python functions
    ## Array version 32bit
    start = time.time()
    kmers_int32 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint32)
    distances_c_array_32bit = quickmers.hamming_distance_encoded_array_32bit(target_kmer_int, kmers_int32)
    t_c_array = time.time() - start
    logging.info(f"C library array 32bit time: {t_c_array:.6f} seconds")

    ## Array version 64bit
    start = time.time()
    kmers_int64 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint64)
    distances_c_array_64bit = quickmers.hamming_distance_encoded_array_64bit(target_kmer_int, kmers_int64)
    t_c_array = time.time() - start
    logging.info(f"C library array 64bit time: {t_c_array:.6f} seconds")

    ## One by one version 32bit
    distances_c_single_32bit = []
    start = time.time()
    kmers_int32 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint32)
    for kmer in kmers_int32:
        distances_c_single_32bit.append(quickmers.hamming_distance_encoded_32bit(target_kmer_int, kmer))
    t_c_single = time.time() - start
    logging.info(f"C library 32bit one by one time: {t_c_single:.6f} seconds")

    ## One by one version 64bit
    distances_c_single_64bit = []
    start = time.time()
    kmers_int64 = np.array([string2bitstring(km) for km in kmers_str], dtype=np.uint64)
    for kmer in kmers_int64:
        distances_c_single_64bit.append(quickmers.hamming_distance_encoded_64bit(int(target_kmer_int), int(kmer)))
    t_c_single = time.time() - start
    logging.info(f"C library 64bit one by one time: {t_c_single:.6f} seconds")

    # Encoded in C functions directly from strings
    ## Array version 32bit
    start = time.time()
    distances_c_array_32bit_strings = quickmers.hamming_distance_array_32bit(target_kmer_str, kmers_str)
    t_c_array = time.time() - start
    logging.info(f"C library array 32bit for strings time: {t_c_array:.6f} seconds")

    ## Array version 64bit
    start = time.time()
    distances_c_array_64bit_strings = quickmers.hamming_distance_array_64bit(target_kmer_str, kmers_str)
    t_c_array = time.time() - start
    logging.info(f"C library array 64bit for strings time: {t_c_array:.6f} seconds")

    ## One by one version 32bit
    distances_c_single_32bit_string = []
    start = time.time()
    for kmer in kmers_str:
        distances_c_single_32bit_string.append(quickmers.hamming_distance_32bit(target_kmer_str, kmer))
    t_c_single = time.time() - start
    logging.info(f"C library 32bit one by one for strings time: {t_c_single:.6f} seconds")

    ## One by one version 64bit
    distances_c_single_64bit_string = []
    start = time.time()
    for kmer in kmers_str:
        distances_c_single_64bit_string.append(quickmers.hamming_distance_64bit(target_kmer_str, kmer))
    t_c_single = time.time() - start
    logging.info(f"C library 64bit one by one for strings time: {t_c_single:.6f} seconds")

    # -----------------------------
    # Result check
    # -----------------------------
    if np.all(distances_py == distances_c_array_32bit):
        logging.info("Correct results for encoded 32bit array results ✅")
    else:
        logging.warning("Mismatch for encoded 32bit array results ❌")

    if np.all(distances_py == distances_np):
        logging.info("Correct results for NumPy results ✅")
    else:
        logging.warning("Mismatch for NumPy results ❌")

    if np.all(distances_py == distances_c_single_32bit):
        logging.info("Correct results for encoded 32bit one by one results ✅")
    else:
        logging.warning("Mismatch for encoded 32bit one by one results ❌")

    if np.all(distances_py == distances_c_single_64bit):
        logging.info("Correct results for encoded 64bit one by one results ✅")
    else:
        logging.warning("Mismatch for encoded 64bit one by one results ❌")

    if np.all(distances_py == distances_c_array_64bit):
        logging.info("Correct results for encoded 64bit array results ✅")
    else:
        logging.warning("Mismatch for encoded 64bit array results ❌")

    if np.all(distances_py == distances_c_array_32bit_strings):
        logging.info("Correct results for 32bit array results ✅")
    else:
        logging.warning("Mismatch for C 32bit array results ❌")

    if np.all(distances_py == distances_c_array_64bit_strings):
        logging.info("Correct results for 64bit array results ✅")
    else:
        logging.warning("Mismatch for C 64bit array results ❌")

    if np.all(distances_py == distances_c_single_32bit_string):
        logging.info("Correct results for 32bit one by one results ✅")
    else:
        logging.warning("Mismatch for 32bit one by one results ❌")

    if np.all(distances_py == distances_c_single_64bit_string):
        logging.info("Correct results for 64bit one by one results ✅")
    else:
        logging.warning("Mismatch for 64bit one by one results ❌")

    # For extra logging in case of mismatches
    # logging.info(target_kmer_str)
    # logging.info(kmers_str)
    # logging.info(distances_c_array_32bit)
    # logging.info(distances_c_array_32bit_strings)
    # logging.info(distances_c)
    # logging.info(distances_np)
    # logging.info(distances_py)

# -----------------------------
# Run benchmark
# -----------------------------
if __name__ == "__main__":
    benchmark(n=1_000_000, k=16)
