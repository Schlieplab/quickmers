import quickmers

def test_levenshtein_basic():
    assert quickmers.levenshtein("GATTACA", "GACTATA") == 2