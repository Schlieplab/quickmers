import quickmers
import numpy as np

def test_levenshtein_distance_array():
    query = "GATTACA"
    kmers = ["GATTACA", "GACTATA", "GATTTCA"]
    result = quickmers.levenshtein_distance_array(query, kmers)

    assert result.tolist() == [0, 2, 1]

def test_levenshtein_distance_array_numpy():
    import numpy as np
    query = "GATTACA"
    kmers = np.array(["GATTACA", "GACTATA", "GATTTCA"], dtype=object)
    result = quickmers.levenshtein_distance_array(query, kmers)

    assert result.tolist() == [0, 2, 1]