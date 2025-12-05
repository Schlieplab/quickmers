# QuickMers

[![PyPI Version](https://img.shields.io/pypi/v/quickmers)](https://pypi.org/project/quickmers)
[![License](https://img.shields.io/badge/License-LGPL--3.0--or--later-blue)](LICENSE)

---

## Overview

**QuickMers** is a high-performance library for computing distances between k-mers (short DNA sequences), implemented in C with Python bindings. It provides fast functions for:

- **Hamming distance** – number of mismatches between equal-length sequences.
- **Levenshtein (edit) distance** – minimum number of insertions, deletions, or substitutions needed to transform one sequence into another.

The library leverages **bitwise operations for Hamming distance** and the **Myers bit-parallel algorithm** for edit distance. On systems with **AVX2 instructions**, computations are further accelerated.

---

## Installation

Install from PyPI:

```bash
pip install quickmers
```

Or install manually:

```bash
python setup.py install
```

---

## Usage

```python
import numpy as np
from quickmers import (
    hamming_distance_array,
    hamming_distance,
    levenshtein,
    levenshtein_array,
    levenshtein_array_with_min_dist
)
```
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

---

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

---

### 3. Levenshtein Distance

#### Signature
```python
levenshtein(sequence1: str, sequence2: str) -> int
```

Compute the edit distance between two sequences.

```python
d = levenshtein("ACGT", "CGTA")
print(d)  # Output: 2
```

---

### 4. Levenshtein Distance (array)

#### Signature
```python
levenshtein_array(query: str, targets: Union[List[str], np.ndarray]) -> np.ndarray
```

Compute edit distances between a query sequence and a list of sequences.

```python
query = "ACGT"
targets = ["ACGT", "TCGA", "CGTA"]
distances = levenshtein_array(query, targets)
print(distances)  # Output: array([0, 2, 2])
```

---

### 5. Levenshtein Distance with Minimum Distance

#### Signature
```python
levenshtein_array_with_min_dist(query: str, targets: Union[List[str], np.ndarray], min_distance: int) -> Tuple[bool, List[Optional[int]]]
```
Compute the Levenshtein (edit) distance between a query k-mer and a list of target k-mers, with **early exit** if a minimum distance is violated.  

This function is useful when you want to skip unnecessary calculations once a target is "too close" to the query, which can save significant computation time for long lists.

```python
query = "ACGT"
targets = ["AGGG", "TCGA", "CGTA"]
early_exit, distances = levenshtein_array_with_min_dist(query, targets, 1)
print(early_exit) # Output: False
print(distances)  # Output: array([2, 2, 2])
```
---

## Implementation Details
- **Hamming distance:** computed using bitwise operations on 64-bit integers for speed. The sequence length has to be 32 or lower.
- **Levenshtein distance:** implemented with **Myers bit-parallel algorithm** for rapid edit distance computation. The sequence length has to be 64 or lower.
- **Array functions:** can take Python lists or NumPy object arrays as input.
- **AVX2 acceleration:** if your CPU supports AVX2, internal loops are vectorized for faster computation.
- **Memory management:** all functions use efficient pre-allocated arrays to minimize Python overhead.

---

## Performance

Benchmarking was performed on typical DNA sequences using our C backend:

| Function | Sequence count | Length | Time per 1e6 pairs (s) |
|----------|----------------|--------|-----------------------|
| `hamming_distance` | 1 | 100 | 0.00001 |
| `hamming_distance_array` | 1000 | 100 | 0.002 |
| `levenshtein` | 1 | 100 | 0.0001 |
| `levenshtein_array` | 1000 | 100 | 0.03 |
| `levenshtein_array_with_min_dist` | 1000 | 100 | 0.04 |

> Results may vary depending on CPU and whether AVX2 instructions are available.

---

## References

[1] Gene Myers. 1999. *A fast bit-vector algorithm for approximate string matching based on dynamic programming*. J. ACM 46, 3 (May 1999), 395-415. https://doi.org/10.1145/316542.316550

---

## License

This software is licensed under **LGPL-3.0-or-later**.

Copyright 2025 Alexander Schliep
