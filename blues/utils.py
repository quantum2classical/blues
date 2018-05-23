import os
from os.path import abspath, dirname, join


def get_directory(directory):
    """
    Retrieves absolute path to specified directory.

    Parameters
    ----------
    directory : str
        Name of directory

    Returns
    -------
    path_to_dir : str
        Absolute path to desired directory.
    """
    if not isinstance(directory, str):
        raise TypeError("Input must be a string.")
    else:
        pass

    path_to_dir = abspath(join(dirname(__file__), directory))

    if not os.path.exists(path_to_dir):
        raise OSError(f"{directory} does not exist.")
    else:
        pass

    return path_to_dir
