import gzip
import bz2

def adaptiveOpen(path, mode = "r"):
    """
    
    Arguments:
    - `path`:
    - `mode`:
    """

    if path.endswith('.gz'):
        return gzip.open(path, mode)
    elif path.endswith('.bz2'):
        return bz2.BZ2File(path, mode)
    return open(path, mode)

def FileType(mode):
    """
    
    Arguments:
    - `mode`:
    """

    return lambda path: adaptiveOpen(path, mode)
