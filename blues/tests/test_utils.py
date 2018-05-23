import os

from ..utils import get_directory


def test_get_directory():

    # Test bad input
    try:
        _ = get_directory("fake_dir")
    except(OSError):
        pass
    else:
        raise Exception("Fake directory name was accepted.")

    try:
        _ = get_directory(1234)
    except(TypeError):
        pass
    else:
        raise Exception("Numeric value accepted.")

    assert get_directory("tests") == os.path.dirname(__file__),\
        "Path does not match."

    return
