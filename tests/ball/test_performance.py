#!/usr/bin/env python3
import time
import multiprocessing as mp
import quickmers

# -----------------------------
# Benchmark parameters
# -----------------------------
K_LENGTHS = [8, 12, 16, 20, 24, 28]
RADII = [2, 3, 4, 5]
N_REPEAT = 5  # how many times to repeat each measurement and take median
MAX_WORKERS = None  # None => use mp.cpu_count(); set to e.g. 8 if you want to limit

# -----------------------------
# Helpers
# -----------------------------
def deterministic_kmer(k: int) -> str:
    """Generate kmer as repeated 'ACGT', truncated to length k."""
    base = "ACGT"
    repeats = (k // len(base)) + 1
    return (base * repeats)[:k]


def time_fixed_ball(kmer: str, radius: int) -> float:
    """Time the list version: fixed_length_levenshtein_ball."""
    start = time.perf_counter()
    _ = quickmers.fixed_length_levenshtein_ball(kmer, radius, 0)
    return time.perf_counter() - start


def time_fixed_ball_iterator(kmer: str, radius: int) -> float:
    """Time the iterator version: fixed_length_levenshtein_ball_iterator (just iterate)."""
    start = time.perf_counter()
    it = quickmers.fixed_length_levenshtein_ball_iterator(kmer, radius, 0)
    count = 0
    for _ in it:
        count += 1  # just consume
    return time.perf_counter() - start


def benchmark(func, kmer: str, radius: int) -> float:
    """Run func multiple times and return median time."""
    times = []
    for _ in range(N_REPEAT):
        t = func(kmer, radius)
        times.append(t)
    times.sort()
    return times[len(times) // 2]  # median


# -----------------------------
# Worker for multiprocessing
# -----------------------------
def run_benchmark_task(args):
    """
    args: (func_name, k, radius)
    returns: (func_name, k, radius, median_time)
    """
    func_name, k, radius = args
    kmer = deterministic_kmer(k)

    if func_name == "fixed_length_levenshtein_ball":
        t = benchmark(time_fixed_ball, kmer, radius)
    elif func_name == "fixed_length_levenshtein_ball_iterator":
        t = benchmark(time_fixed_ball_iterator, kmer, radius)
    else:
        raise ValueError(f"Unknown function name: {func_name}")

    return func_name, k, radius, t


# -----------------------------
# Performance test
# -----------------------------
def test_performance_fixed_length_balls():
    func_names = [
        "fixed_length_levenshtein_ball",
        "fixed_length_levenshtein_ball_iterator",
    ]

    # Build all tasks: one per (func, k, radius)
    tasks = []
    for fn in func_names:
        for k in K_LENGTHS:
            for r in RADII:
                tasks.append((fn, k, r))

    print(f"\nRunning performance benchmarks with {len(tasks)} tasks... (this will take a while)")

    with mp.Pool(processes=MAX_WORKERS) as pool:
        results_raw = pool.map(run_benchmark_task, tasks)

    # Organize results: results[func_name][k][radius] = median time
    results = {fn: {} for fn in func_names}
    for fn, k, r, t in results_raw:
        results[fn].setdefault(k, {})[r] = t

    # ---- Print Markdown tables ----
    print("\n### Benchmark Results: fixed_length_levenshtein_ball (list version)")
    print("Times are median over N_REPEAT={} runs.\n".format(N_REPEAT))

    header = "| k \\ radius | " + " | ".join(f"r={r}" for r in RADII) + " |"
    sep = "|" + "------------|" * (len(RADII) + 1)

    print(header)
    print(sep)

    for k in K_LENGTHS:
        row_times = [f"{results['fixed_length_levenshtein_ball'][k][r]:.6f}" for r in RADII]
        print(f"| {k} | " + " | ".join(row_times) + " |")

    print("\n### Benchmark Results: fixed_length_levenshtein_ball_iterator")
    print("Times are median over N_REPEAT={} runs.\n".format(N_REPEAT))

    print(header)
    print(sep)

    for k in K_LENGTHS:
        row_times = [f"{results['fixed_length_levenshtein_ball_iterator'][k][r]:.6f}" for r in RADII]
        print(f"| {k} | " + " | ".join(row_times) + " |")


if __name__ == "__main__":
    # On Windows / some environments, this guard is important:
    mp.set_start_method("spawn", force=True)
    test_performance_fixed_length_balls()