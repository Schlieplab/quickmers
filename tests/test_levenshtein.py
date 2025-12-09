import quickmers

def test_levenshtein_basic():
    assert quickmers.levenshtein_distance("GATTACA", "GACTATA") == 2