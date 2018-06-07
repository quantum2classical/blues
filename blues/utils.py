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
        raise OSError("'{}' directory does not exist.".format(directory))
    else:
        pass

    return path_to_dir


def file_path(filename, directory):
    """Retrieves absolute path to specified file in a directory.

    Parameters
    ----------
    filename : str
    directory : str

    Returns
    -------
    full_path : str
    """

    full_dirname = get_directory(directory)
    full_path = join(full_dirname, filename)

    if not os.path.exists(full_path):
        raise OSError("{} does not exist.".format(full_path))
    else:
        pass

    return full_path
