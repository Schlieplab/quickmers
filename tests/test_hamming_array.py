import quickmers
import numpy as np

def test_hamming_array():
    query = "ACGT"
    kmers = ["ACGT", "TGCA", "AAAA"]
    result = quickmers.hamming_distance_array(query, kmers)
        
    assert result.tolist() == [0, 4, 3]

def test_hamming_array_numpy():
    import numpy as np
    query = "ACGT"
    kmers = np.array(["ACGT", "TGCA", "AAAA"], dtype=object)
    result = quickmers.hamming_distance_array(query, kmers)

    assert result.tolist() == [0, 4, 3]