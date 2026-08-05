"""
Filename: quickmers/__init__.py
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.1.1
Description: Python bindings for QuickMers C library.
License: LGPL-3.0-or-later
"""
from ._cbindings import (
    hamming_distance_array,
    hamming_distance,
    levenshtein_distance,
    levenshtein_distance_array,
    levenshtein_distance_array_with_min_dist
)