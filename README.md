# QuickMers

[![PyPI Version](https://img.shields.io/pypi/v/quickmers)](https://pypi.org/project/quickmers)
[![License](https://img.shields.io/badge/License-LGPL--3.0--or--later-blue)](LICENSE)


## Overview

**QuickMers** is a high-performance library for computing distances between k-mers (short DNA sequences), implemented in C with Python bindings. It provides fast functions for:

- **Hamming distance** – number of mismatches between equal-length sequences.
- **Levenshtein (edit) distance** – minimum number of insertions, deletions, or substitutions needed to transform one sequence into another.
- **Fixed-length Levenshtein balls** – all k-mers at edit distance ≤ r from a query, restricted to the same length k.
- **Levenshtein balls** – all sequences of any length at edit distance ≤ r from a query.

The library leverages **bitwise operations for Hamming distance** and the **Myers bit-parallel algorithm** for edit distance. On systems with **AVX2 instructions**, computations are further accelerated.


## Installation

**Planned:** Quickmers will be available on [PyPI](https://pypi.org/) in the future.  
Install using pip:
```bash
python -m pip install quickmers
```

or manually:

```bash
python setup.py install
```

## Usage

### 1. Hamming Distance

#### Signature:
```python
hamming_distance(sequence1: str, sequence2: str) -> int
```

Compute the Hamming distance between two sequences of equal length.

```python
d = hamming_distance("ACGT", "TCGA")
print(d)  # Output: 2
```


### 2. Hamming Distance (array)

#### Signature
```python
hamming_distance_array(query: str, targets: Union[List[str], np.ndarray]) -> np.ndarray
```

Compute Hamming distances between one query sequence and a list of target sequences.

```python
query = "ACGT"
targets_list = ["ACGT", "TCGA", "CGTA"]
targets_array = np.array(targets_list, dtype=object)

distances1 = hamming_distance_array(query, targets_list)
distances2 = hamming_distance_array(query, targets_array)

print(distances1) # Output: array([0, 2, 4])
print(distances2) # Output: array([0, 2, 4])
```


### 3. Levenshtein Distance

#### Signature
```python
levenshtein_distance(sequence1: str, sequence2: str) -> int
```

Compute the edit distance between two sequences.

```python
d = levenshtein_distance("ACGT", "CGTA")
print(d)  # Output: 2
```

### 4. Levenshtein Distance (array)

#### Signature
```python
levenshtein_distance_array(query: str, targets: Union[List[str], np.ndarray]) -> np.ndarray
```

Compute edit distances between a query sequence and a list of sequences.

```python
query = "ACGT"
targets = ["ACGT", "TCGA", "CGTA"]
distances = levenshtein_distance_array(query, targets)
print(distances)  # Output: array([0, 2, 2])
```

### 5. Levenshtein Distance with Minimum Distance

#### Signature
```python
levenshtein_distance_array_with_min_dist(query: str, targets: Union[List[str], np.ndarray], min_distance: int) -> Tuple[bool, List[Optional[int]]]
```
Compute the Levenshtein (edit) distance between a query k-mer and a list of target k-mers, with **early exit** if a minimum distance is violated.  

This function is useful when you want to skip unnecessary calculations once a target is "too close" to the query, which can save significant computation time for long lists.

```python
query = "ACGT"
targets = ["AGGG", "TCGA", "CGTA"]
early_exit, distances = levenshtein_distance_array_with_min_dist(query, targets, 1)
print(early_exit) # Output: False
print(distances)  # Output: array([2, 2, 2])
```


### 6. Fixed-length Levenshtein ball

A **Levenshtein ball** of radius `r` around a k-mer `x` is the set of all k-mers whose edit distance to `x` is at most `r`.  
In many k-mer applications we only care about sequences of the **same length** as `x`. The **fixed-length Levenshtein ball** is therefore:

```text
B_r(x) = { y in {A, C, G, T}^k | edit_distance(x, y) <= r }
```

All output k-mers have length exactly `k`, the same as the input k-mer. This is different from a standard Levenshtein neighborhood, which can include shorter or longer strings due to insertions and deletions.

QuickMers provides two interfaces:

- A function that returns the entire ball as a Python list.
- An iterator that yields members of the ball one by one, useful for large radii or memory-constrained settings.

Both functions assume:
- Alphabet: `A`, `C`, `G`, `T`
- Input k-mer length `k` ≤ 32 (due to the underlying bit-parallel implementation).

#### Signature (list version)
```python
fixed_length_levenshtein_ball(kmer: str, radius: int, return_binary: int = 0) -> List[str]
```

- `kmer`: input k-mer (e.g. `"ACGTACGT"`).
- `radius`: maximum Levenshtein distance.
- `return_binary`: if `0` (default), returns k-mers as strings; if `1`, returns internal 64-bit integer encodings.

Example:

```python
kmer = "ACGT"
radius = 2

ball = fixed_length_levenshtein_ball(kmer, radius, 0)
print(len(ball))        # Output: 85
print(ball[:5])         # first 5 k-mers in the ball, output: ['CAGT', 'CGGT', 'CTGT', 'GAGT', 'GGGT']
```

All k-mers in `ball` have length `len(kmer)` and satisfy `levenshtein_distance(kmer, y) <= radius`.


### 7. Fixed-length Levenshtein ball (iterator)

#### Signature
```python
fixed_length_levenshtein_ball_iterator(kmer: str, radius: int, return_binary: int = 0) -> Iterator[str]
```

Returns a Python iterator that yields all k-mers in the fixed-length Levenshtein ball of `kmer` with given `radius`.

- `kmer`: input k-mer.
- `radius`: maximum Levenshtein distance.
- `return_binary`: if `0` (default), yields strings; if `1`, yields 64-bit integer encodings.

Example:

```python
kmer = "ACGT"
radius = 2

it = fixed_length_levenshtein_ball_iterator(kmer, radius, 0)
for y in it:
    # process each k-mer in the ball
    pass
```

This is equivalent to:

```python
ball = fixed_length_levenshtein_ball(kmer, radius, 0)
for y in ball:
    pass
```

but avoids materializing the entire list in memory at once as the size can grow quite fast.


### 8. Levenshtein ball

The **Levenshtein ball** drops the fixed-length restriction of the previous sections. It contains every sequence, of any length, whose edit distance to `x` is at most `r`.

Since insertions and deletions no longer have to cancel each other out, the output contains sequences of every length from `k - r` to `k + r`. The ball is therefore much larger than the fixed-length ball for the same `k` and `r`, and it is a superset of it.

QuickMers again provides two interfaces:

- A function that returns the entire ball as a Python list.
- An iterator that yields members of the ball one by one.

Both functions assume:
- Alphabet: `A`, `C`, `G`, `T`
- `k + r` ≤ 31, i.e. the longest possible output sequence has at most 31 bases (see [Levenshtein ball limitations](#levenshtein-ball-limitations)). A `ValueError` is raised otherwise.

#### Signature (list version)
```python
levenshtein_ball(kmer: str, radius: int, return_binary: int = 0) -> List[str]
```

- `kmer`: input k-mer (e.g. `"ACGTACGT"`).
- `radius`: maximum Levenshtein distance.
- `return_binary`: if `0` (default), returns sequences as strings; if `1`, returns `(encoded, length)` tuples, where `encoded` is the 2-bit encoding of the sequence (`A=00`, `C=01`, `G=10`, `T=11`) and `length` its number of bases. The length is needed because the encoding alone cannot distinguish e.g. `AC` from `AAC`.

Example:

```python
kmer = "ACGT"
radius = 2

ball = levenshtein_ball(kmer, radius)
print(len(ball))        # Output: 431 (the fixed-length ball has 85)
print(ball[:5])         # first 5 sequences in the ball, output: ['CAGT', 'CGGT', 'CTGT', 'GAGT', 'GGGT']

ball = levenshtein_ball(kmer, radius, 1)
print(ball[:3])         # Output: [(75, 4), (107, 4), (123, 4)]
```

The 431 sequences have lengths 2 to 6 (2: 6, 3: 26, 4: 85, 5: 160, 6: 154), and all satisfy `levenshtein_distance(kmer, y) <= radius`.


### 9. Levenshtein ball (iterator)

#### Signature
```python
levenshtein_ball_iterator(kmer: str, radius: int, return_binary: int = 0) -> Iterator[str]
```

Returns a Python iterator that yields all sequences in the Levenshtein ball of `kmer` with given `radius`.

- `kmer`: input k-mer.
- `radius`: maximum Levenshtein distance.
- `return_binary`: if `0` (default), yields strings; if `1`, yields `(encoded, length)` tuples as described above.

Example:

```python
kmer = "ACGT"
radius = 2

it = levenshtein_ball_iterator(kmer, radius)
for y in it:
    # process each sequence in the ball
    pass
```

The iterator avoids building the Python list of all sequences, so it needs less memory than `levenshtein_ball`. It still has to remember every sequence it has already yielded to avoid duplicates, so its memory use still grows with the size of the ball.


## Implementation Details
- **Hamming distance:** computed using bitwise operations on 64-bit integers for speed. The sequence length has to be 32 or lower.
- **Levenshtein distance:** implemented with **Myers bit-parallel algorithm** for rapid edit distance computation. The sequence length has to be 64 or lower.
- **Array functions:** can take Python lists or NumPy object arrays as input.
- **AVX2 acceleration:** if your CPU supports AVX2, internal loops are vectorized for faster computation.
- **Memory management:** all functions use efficient pre-allocated arrays to minimize Python overhead.
- **Fixed-length Levenshtein balls:** generated using a combinatorial enumeration of operation strings (matches, mismatches, insertions, deletions) that preserve output length `k`. Only k-mers of length exactly `k` are produced; no shorter or longer sequences appear. Supported for `k` ≤ 32.
- **Levenshtein balls:** same enumeration, but without requiring the number of insertions and deletions to be equal. Results are deduplicated in a hash set, using 64-bit integers that carry an extra bit marking the sequence length. Supported for `k + r` ≤ 31.


## Performance

Benchmarking was performed on randomly generated k-mers of varying lengths.  
Values in the first table represent the number of string pairs distances calculated per second.  
Values in the Levenshtein ball tables represent median runtime in seconds over 5 runs.

### Distance function benchmarks

| Function | k=5 | k=10 | k=15 | k=20 | k=25 | k=30 |
|----------|----------|----------|----------|----------|----------|----------|
| `hamming_distance` | 9,505,940 | 9,192,009 | 8,866,906 | 8,555,942 | 8,249,809 | 7,946,207 |
| `hamming_distance_array` | 49,038,568 | 47,807,278 | 44,236,469 | 38,861,603 | 35,788,375 | 33,238,156 |
| `levenshtein_distance` | 7,499,375 | 6,407,775 | 5,797,060 | 5,131,730 | 4,688,698 | 4,267,584 |
| `levenshtein_distance_array` | 26,354,634 | 20,227,311 | 16,396,735 | 13,399,246 | 11,772,501 | 10,360,021 |


### Fixed-length Levenshtein ball benchmarks

The following tables show median runtimes (over 5 runs) for generating fixed-length Levenshtein balls of radius `r` around a k-mer of length `k`. All times are in seconds.
Note that, unlike Hamming balls, fixed-length Levenshtein ball sizes depend on the k-mer’s sequence: different k-mers of the same length can yield different ball sizes for the same radius.

#### `fixed_length_levenshtein_ball` (list version)

| k \ radius | r=2 | r=3 | r=4 | r=5 |
|------------|------|------|------|-------|
| 8  | 0.000120 | 0.001434 | 0.007199 | 0.043119 |
| 12 | 0.000440 | 0.007060 | 0.096053 | 0.985557 |
| 16 | 0.000918 | 0.013795 | 0.491295 | 6.161375 |
| 20 | 0.001075 | 0.034931 | 1.348028 | 29.785755 |
| 24 | 0.001494 | 0.093763 | 3.741409 | 76.865599 |
| 28 | 0.003478 | 0.173104 | 7.473514 | 207.178461 |

#### `fixed_length_levenshtein_ball_iterator` (iterator version)

| k \ radius | r=2 | r=3 | r=4 | r=5 |
|------------|------|------|------|-------|
| 8  | 0.000276 | 0.002309 | 0.019280 | 0.085768 |
| 12 | 0.000640 | 0.009446 | 0.119095 | 1.479240 |
| 16 | 0.001180 | 0.016352 | 0.536328 | 7.978017 |
| 20 | 0.001077 | 0.071314 | 1.441640 | 30.129923 |
| 24 | 0.002852 | 0.076538 | 3.723360 | 84.883196 |
| 28 | 0.003997 | 0.167211 | 7.499478 | 209.908604 |

### Levenshtein ball benchmarks

Median runtimes (over 5 runs) for generating Levenshtein balls (all lengths) of radius `r` around a k-mer of length `k`. All times are in seconds. `n/a` marks combinations that are not supported (`k + r` > 31). `r = 5` is not benchmarked, see the limitations below.

#### `levenshtein_ball` (list version)

| k \ radius | r=2 | r=3 | r=4 |
|------------|------|------|------|
| 8  | 0.000450 | 0.013566 | 0.193951 |
| 12 | 0.001666 | 0.071871 | 1.160621 |
| 16 | 0.003254 | 0.204384 | 4.177346 |
| 20 | 0.005811 | 0.328056 | 10.125348 |
| 24 | 0.009441 | 0.564906 | 25.353289 |
| 28 | 0.014549 | 1.107160 | n/a |

#### `levenshtein_ball_iterator` (iterator version)

| k \ radius | r=2 | r=3 | r=4 |
|------------|------|------|------|
| 8  | 0.000629 | 0.019127 | 0.254369 |
| 12 | 0.001352 | 0.084010 | 1.373834 |
| 16 | 0.004341 | 0.140489 | 4.661887 |
| 20 | 0.006934 | 0.318951 | 12.333837 |
| 24 | 0.011000 | 0.601258 | 24.105267 |
| 28 | 0.015928 | 1.031021 | n/a |

#### Levenshtein ball limitations

> **Memory and runtime grow very quickly with the radius.** Only use the Levenshtein ball functions for small radii; `r = 5` and above is impractical for most k.

- **Ball size and memory.** The ball grows exponentially in `r`. For `k = 24` and `r = 4` it already contains about 27 million sequences, and generating it with `levenshtein_ball` needs about 4.5 GB of memory. Increasing the radius to 5 grows the ball roughly 20–30×, which would need an estimated 100 GB or more. Every sequence of the ball is kept in memory for deduplication, in both the list and the iterator version, so the iterator reduces but does not remove this cost.
- **Runtime.** For each combination of edit operations, every insertion can be any of 4 bases and every substitution any of 3. The fixed-length ball has at most `r / 2` insertions, while the Levenshtein ball can have up to `r`. At `r = 4` this makes it about 8× slower than the fixed-length ball (see the tables above). At `r = 5` the difference grows to about 11×, which for `k = 24` means an estimated 14 minutes per call.
- **Why `k = 28`, `r = 4` is `n/a`.** Sequences are stored as 64-bit integers with 2 bits per base, plus one extra bit that marks the sequence length (so that e.g. `AC` and `AAC` are distinguishable). This allows sequences of at most 31 bases. The longest sequence in a ball is `k + r` bases long (`r` insertions), so `k + r` must be at most 31. For `k = 28`, `r = 4` this is 32, and the functions raise a `ValueError`. The fixed-length functions do not need the extra bit and support `k` ≤ 32.

> Results may vary depending on CPU architecture, compiler optimizations, and whether AVX2 instructions are available.


## References

[1] Gene Myers. 1999. *A fast bit-vector algorithm for approximate string matching based on dynamic programming*. J. ACM 46, 3 (May 1999), 395-415. [https://doi.org/10.1145/316542.316550](https://doi.org/10.1145/316542.316550)

## License

This software is licensed under **LGPL-3.0-or-later**.

Copyright 2025 Alexander Schliep