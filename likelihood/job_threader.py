from multiprocessing import Pool
import atexit

__author__ = 'davidmurphy'

"""
Usage:
  - Init with init_threads();
  - Calc likelihood in parallel with runParallel().
"""


def init_threads(fixp):
    """create a list of Pool objects to process each chunk of data"""
    global _pools
    _pools = [Pool(1, initializer=_setfixp, initargs=(fp,)) for fp in fixp]
    # schedule cleanup of these resources when program is finished
    atexit.register(_cleanup)


def _setfixp(fp):
    """set the fixed parameters for one pool"""
    global _fixp
    _fixp = fp


def _fpart(func, var):
    """evaluate the function for the subset of data in _fixp"""
    args = tuple([var] + [fp for fp in _fixp])
    # evaluate args via tuple expansion
    return func(*args)


def run_threads(func, var):
    """evaluate the function for var at each subset of data in _fixp"""
    res = [p.apply_async(_fpart, args=(func, var)) for p in _pools]
    return [r.get() for r in res]


def _cleanup():
    """order processes to terminate and wait until they do"""
    for p in _pools:
        p.close()
        p.join()

