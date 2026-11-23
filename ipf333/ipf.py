"""Iterative proportional fitting (IPF, raking) for 3x3x3 tables in arbitrary precision, with its raking factors."""
import itertools
import mpmath as mp

CELLS = list(itertools.product(range(3), repeat=3))


def _margin(q, ax):
    out = {}
    for c in CELLS:
        key = tuple(x for a, x in enumerate(c) if a != ax)
        out[key] = out.get(key, 0) + q[c]
    return out


def ipf(n, start=None, dps=60, tol_digits=None, max_sweeps=100000):
    """Fit the no-three-way interaction model to the table n (27 counts, cell ijk at index 9i+3j+k, 0-based) by IPF
    from the start table (default: uniform). Each sweep matches n_ij+, then n_i+k, then n_+jk. Returns (p, (a, b, c), sweeps): the fitted
    table and the raking factors, the cumulative products of the adjustments applied to each margin, so that
    p_ijk = t_ijk a_ij b_ik c_jk exactly, t the start table. All are dicts keyed by 0-based indices, with mpf values."""
    tol = mp.mpf(10) ** -(tol_digits or dps - 5)
    with mp.workdps(dps):
        N = dict(zip(CELLS, n)); tgt = [_margin(N, ax) for ax in range(3)]
        p = {c: mp.mpf(1) if start is None else mp.mpf(start[9 * c[0] + 3 * c[1] + c[2]]) for c in CELLS}
        F = [{k: mp.mpf(1) for k in tgt[ax]} for ax in range(3)]   # F[2] = a_ij, F[1] = b_ik, F[0] = c_jk
        for sweep in range(1, max_sweeps + 1):
            for ax in (2, 1, 0):
                m = _margin(p, ax); r = {k: tgt[ax][k] / m[k] for k in m}
                for k in r: F[ax][k] *= r[k]
                for c in CELLS:
                    p[c] *= r[tuple(x for a, x in enumerate(c) if a != ax)]
            err = max(abs(_margin(p, ax)[k] - tgt[ax][k]) for ax in range(3) for k in tgt[ax])
            if err < tol:
                return p, (F[2], F[1], F[0]), sweep
    raise RuntimeError("IPF did not converge")


def factorization_residual(p, factors, start=None):
    """max |t_ijk a_ij b_ik c_jk / p_ijk - 1|."""
    a, b, c = factors
    t = lambda i, j, k: 1 if start is None else start[9 * i + 3 * j + k]
    return max(abs(t(i, j, k) * a[(i, j)] * b[(i, k)] * c[(j, k)] / p[(i, j, k)] - 1) for (i, j, k) in CELLS)
