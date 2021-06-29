import atexit
from multiprocessing import Pipe, Process

__author__ = 'davidmurphy'


"""
Usage:
  - Init with init_threads(nproc);
  - Set fixed parameters with set_fixed_params(fixed_params)
  - Calc likelihood in parallel with runParallel(var).
"""


# GLOBAL MODULE PARAMETERS
_state = 'null_state'
_nproc = 0
_func = lambda x: x
_p, _q = [], []


def set_state(new_state='null_state'):
    """reset state to null condition"""
    global _state
    _state = new_state


def init_threads(nproc, func):
    """create list of Pipes to process each chunk of data"""
    # the current state of the threader should be 0 at initialization
    global _state
    assert _state == 'null_state'

    # initialize the other global variables
    global _nproc, _func, _p, _q
    _nproc = nproc
    _func = func
    _p, _q, = [], []

    # create processes by opening pipes
    for _ in xrange(_nproc):
        parent_conn, child_conn = Pipe()
        p = Process(target=_worker, args=(child_conn,))
        p.start()
        _p.append(p)
        _q.append(parent_conn)

    # reset _state to indicate that initialization is complete
    set_state('initialized')

    # schedule cleanup of these resources when program is finished
    atexit.register(_cleanup)


def set_fixed_params(fixed_params):
    """set the fixed parameters for each process"""
    global _state
    assert _state == 'initialized'
    assert len(fixed_params) == _nproc

    # have each parent send indexed fixed param to its child
    for (i, fp) in enumerate(fixed_params):
        _q[i].send(('set_fixed_params', fp))

    # reset global _state to indicate that fixed params are now loaded
    set_state('fixed_params_loaded')


def run_threads(inp):
    """evaluate the function for inp at each subset of data"""
    global _state
    assert _state == 'fixed_params_loaded'

    # have each parent send the calculate command to its child
    for q in _q:
        q.send(('calculate', inp))
    calculations = []
    for q in _q:
        check_inp, val = q.recv()
        try:
            assert all(a == b for (a, b) in zip(inp, check_inp))
        except AssertionError:
            msg = 'inp={} check_inp={}'.format(str(inp), str(check_inp))
            raise AssertionError(msg)
        calculations.append(val)
    return calculations


def _worker(conn):
    """child function that processes data and commands from parent"""
    fixp = None
    while True:
        cmd, inp = conn.recv()
        if cmd == 'set_fixed_params':
            fixp = inp
        elif cmd == 'calculate':
            args = tuple([inp] + [fp for fp in fixp])
            res = _func(*args)
            conn.send((inp, res))
        elif cmd == 'terminate':
            conn.close()
            break


def _cleanup():
    """terminate all processes and close all pipes"""
    for q in _q:
        q.send(('terminate', None))
    for p in _p:
        p.join()
