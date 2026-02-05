#!/usr/bin/env python3
import random
import time
import quickmers

K_LENGTHS = [5, 10, 15, 20, 25, 30]
N_PAIRS = 1_000_000

FUNCS = [
    ("hamming_distance", quickmers.hamming_distance),
    ("hamming_distance_array", quickmers.hamming_distance_array),
    ("levenshtein_distance", quickmers.levenshtein_distance),
    ("levenshtein_distance_array", quickmers.levenshtein_distance_array),
]

def random_kmer(k):
    return "".join(random.choice("ACGT") for _ in range(k))

def benchmark_function(func, *args):
    """Return time in seconds for processing all pairs."""
    start = time.perf_counter()
    func(*args)
    return time.perf_counter() - start

def benchmark_single(func, k):
    """Benchmark scalar version: call f(a, b) N times."""
    a = random_kmer(k)
    b = random_kmer(k)
    start = time.perf_counter()
    for _ in range(N_PAIRS):
        func(a, b)
    end = time.perf_counter()
    elapsed = end - start
    return N_PAIRS / elapsed  # calculations per second

def benchmark_array(func, k):
    """Benchmark array version: call f(query, list_of_kmers) once."""
    query = random_kmer(k)
    kmers = [random_kmer(k) for _ in range(N_PAIRS)]
    elapsed = benchmark_function(func, query, kmers)
    return N_PAIRS / elapsed  # calculations per second

def test_performance():
    results = {name: [] for name, _ in FUNCS}

    print("Running benchmarks... (this will take time)")

    for k in K_LENGTHS:
        for name, fn in FUNCS:
            if "array" in name:
                cps = benchmark_array(fn, k)
            else:
                cps = benchmark_single(fn, k)
            results[name].append(cps)

    # ---- Print Markdown table ----
    print("\n### Benchmark Results")
    print(f"Calculations for N_PAIRS={N_PAIRS}...")

    header = "| Function | " + " | ".join(f"k={k}" for k in K_LENGTHS) + " |"
    sep = "|" + "----------|" * (len(K_LENGTHS) + 1)

    print(header)
    print(sep)

    for name in results:
        row = [f"{int(t):,}" for t in results[name]]
        print("| `{}` | ".format(name) + " | ".join(row) + " |")
