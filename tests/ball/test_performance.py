#!/usr/bin/env python3
import time
import multiprocessing as mp
import quickmers

# -----------------------------
# Benchmark parameters
# -----------------------------
K_LENGTHS = [8, 12, 16, 20, 24, 28]
RADII = [2, 3, 4, 5]
# levenshtein_ball requires k + r <= 31 (such cells are shown as n/a), and
# r=5 needs ~100+ GB per task (the ball holds hundreds of millions of k-mers),
# so its tables only go up to MAX_VARIABLE_RADIUS.
MAX_VARIABLE_LEN = 31
MAX_VARIABLE_RADIUS = 4
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


def time_variable_ball(kmer: str, radius: int) -> float:
    """Time the list version: levenshtein_ball."""
    start = time.perf_counter()
    _ = quickmers.levenshtein_ball(kmer, radius, 0)
    return time.perf_counter() - start


def time_variable_ball_iterator(kmer: str, radius: int) -> float:
    """Time the iterator version: levenshtein_ball_iterator (just iterate)."""
    start = time.perf_counter()
    it = quickmers.levenshtein_ball_iterator(kmer, radius, 0)
    count = 0
    for _ in it:
        count += 1  # just consume
    return time.perf_counter() - start


# func_name -> (timing function, table title)
BENCHMARKS = {
    "fixed_length_levenshtein_ball": (time_fixed_ball, "fixed_length_levenshtein_ball (list version)"),
    "fixed_length_levenshtein_ball_iterator": (time_fixed_ball_iterator, "fixed_length_levenshtein_ball_iterator"),
    "levenshtein_ball": (time_variable_ball, "levenshtein_ball (list version)"),
    "levenshtein_ball_iterator": (time_variable_ball_iterator, "levenshtein_ball_iterator"),
}


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

    if func_name not in BENCHMARKS:
        raise ValueError(f"Unknown function name: {func_name}")
    t = benchmark(BENCHMARKS[func_name][0], kmer, radius)

    return func_name, k, radius, t


def run_benchmarks(func_names, radii, is_supported=lambda k, r: True):
    """
    Run all (func, k, radius) tasks in parallel and print one Markdown table
    per function. Unsupported (k, radius) cells are printed as n/a.
    """
    # Build all tasks: one per (func, k, radius)
    tasks = []
    for fn in func_names:
        for k in K_LENGTHS:
            for r in radii:
                if is_supported(k, r):
                    tasks.append((fn, k, r))

    print(f"\nRunning performance benchmarks with {len(tasks)} tasks... (this will take a while)")

    # Longest tasks (large k and r) first so they don't end up running last
    tasks.sort(key=lambda t: (t[2], t[1]), reverse=True)

    with mp.Pool(processes=MAX_WORKERS) as pool:
        results_raw = pool.map(run_benchmark_task, tasks, chunksize=1)

    # Organize results: results[func_name][k][radius] = median time
    results = {fn: {} for fn in func_names}
    for fn, k, r, t in results_raw:
        results[fn].setdefault(k, {})[r] = t

    # ---- Print Markdown tables ----
    header = "| k \\ radius | " + " | ".join(f"r={r}" for r in radii) + " |"
    sep = "|" + "------------|" * (len(radii) + 1)

    for fn in func_names:
        print(f"\n### Benchmark Results: {BENCHMARKS[fn][1]}")
        print("Times are median over N_REPEAT={} runs.\n".format(N_REPEAT))

        print(header)
        print(sep)

        for k in K_LENGTHS:
            row = results[fn].get(k, {})
            row_times = [f"{row[r]:.6f}" if r in row else "n/a" for r in radii]
            print(f"| {k} | " + " | ".join(row_times) + " |")


# -----------------------------
# Performance tests
# -----------------------------
def test_performance_fixed_length_balls():
    run_benchmarks(
        ["fixed_length_levenshtein_ball", "fixed_length_levenshtein_ball_iterator"],
        RADII,
    )


def test_performance_variable_length_balls():
    run_benchmarks(
        ["levenshtein_ball", "levenshtein_ball_iterator"],
        [r for r in RADII if r <= MAX_VARIABLE_RADIUS],
        is_supported=lambda k, r: k + r <= MAX_VARIABLE_LEN and r <= MAX_VARIABLE_RADIUS,
    )


if __name__ == "__main__":
    # On Windows / some environments, this guard is important:
    mp.set_start_method("spawn", force=True)
    test_performance_fixed_length_balls()
    test_performance_variable_length_balls()
