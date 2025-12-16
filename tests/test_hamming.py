import quickmers
import numpy as np

def test_hamming_simple():
    a = "ACGT"
    b = "TCGA"
    d = quickmers.hamming_distance(a, b)
    assert d == 2
