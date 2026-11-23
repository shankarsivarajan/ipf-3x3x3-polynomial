"""The polynomial behind the IPF limit of a 3x3x3 table.

Given a 3x3x3 table of counts, iterative proportional fitting (raking) from the uniform table converges to the unique
positive table with the same two-way margins of the form a_ij b_ik c_jk. Each fitted cell is a root of a polynomial
of degree 511 with integer coefficients: its values at all complex solutions of the scaling equations. This script
computes the IPF limit, its adjustment factors, and that polynomial for a chosen cell.

Usage:
  python ipf_polynomial.py TABLE [--cell ijk] [--out FILE] [--real] [--msolve PATH]

TABLE is 27 nonnegative integers in the order 111, 112, 113, 121, ..., 333 (the last index varies fastest),
separated by commas or whitespace, or a path to a file holding them. Example:
  python ipf_polynomial.py 3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8 --cell 111 --out E.txt
"""
import argparse, sys, time
from fractions import Fraction
from math import lcm
from pathlib import Path
import mpmath as mp
from ipf333 import ipf, factorization_residual, eliminant

mp.mp.dps = 60


def read27(s, kind):
    txt = Path(s).read_text() if Path(s).is_file() else s
    vals = [Fraction(x) for x in txt.replace(",", " ").split()]
    if len(vals) != 27:
        sys.exit(f"{kind}: expected 27 entries, got {len(vals)}")
    return vals


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("table")
ap.add_argument("--cell", default="111", help="cell whose polynomial is computed (default 111)")
ap.add_argument("--out", help="write the coefficients here, one per line, constant term first")
ap.add_argument("--real", action="store_true", help="also list the real solutions and check that exactly one is positive")
ap.add_argument("--msolve", help="path to the msolve binary (default: msolve on PATH)")
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()

n = read27(a.table, "table")
if any(x.denominator != 1 or x < 0 for x in n):
    sys.exit("table: entries must be nonnegative integers")
n = [int(x) for x in n]
start = None
cell = tuple(int(ch) - 1 for ch in a.cell)
if len(cell) != 3 or any(not 0 <= c <= 2 for c in cell):
    sys.exit("--cell: three digits from 1 to 3, e.g. 213")
ci = 9 * cell[0] + 3 * cell[1] + cell[2]
show = lambda v: mp.nstr(v, 12)

# 1. IPF limit and adjustment factors
p, (A, B, C), sweeps = ipf(n, start=start)
print(f"IPF converged in {sweeps} sweeps (margins matched to 55 digits).")
print("Fitted table, slices i = 1, 2, 3 (rows j, columns k):")
for i in range(3):
    for j in range(3):
        print("   " + "  ".join(f"{show(p[(i, j, k)]):>16}" for k in range(3)))
    print()
print("Adjustment factors, p_ijk = " + ("t_ijk " if start else "") + "a_ij b_ik c_jk (cumulative products of the IPF steps"
      " matching n_ij+, n_i+k, n_+jk):")
for name, D in (("a", A), ("b", B), ("c", C)):
    for r in range(3):
        print(f"  {name if r == 0 else ' '}  " + "  ".join(f"{show(D[(r, s_)]):>16}" for s_ in range(3)))
print(f"  max |{'t ' if start else ''}a b c / p - 1| = {mp.nstr(factorization_residual(p, (A, B, C), start), 3)}\n")

# 2. the polynomial of the chosen cell
x = p[cell]
t0 = time.time()
co = eliminant(n, seed=start, msolve=a.msolve, threads=a.threads, cell=cell)
d = len(co) - 1
print(f"Polynomial of cell {a.cell}: degree {d}, integer coefficients of {min(len(str(abs(c))) for c in co)} to "
      f"{max(len(str(abs(c))) for c in co)} digits ({time.time() - t0:.0f} s).")
if a.out:
    Path(a.out).write_text("\n".join(map(str, co)) + "\n"); print(f"Coefficients written to {a.out} (constant term first).")
with mp.workdps(max(100, 2 * max(len(str(abs(c))) for c in co))):
    X = mp.mpf(x); E = mp.polyval([mp.mpf(c) for c in co[::-1]], X)
    dE = mp.polyval([mp.mpf(k * c) for k, c in list(enumerate(co))[1:][::-1]], X)
    print(f"The fitted value p_{a.cell} = {mp.nstr(x, 30)} is a root: Newton correction |E/E'| = {mp.nstr(abs(E / dE), 3)}.")

# 3. optionally, the real solutions
if a.real:
    from ipf333.realsol import real_solutions, sign_summary
    from ipf333.eliminant import move_cell
    nn = move_cell(n, cell) if cell != (0, 0, 0) else n; ss = None if start is None else (move_cell(start, cell) if cell != (0, 0, 0) else start)
    sols = real_solutions(nn, msolve=a.msolve, threads=a.threads, seed=ss); nr, npos = sign_summary(sols)
    print(f"Real solutions of the scaling equations: {nr}; positive tables among them: {npos} (the IPF limit).")
