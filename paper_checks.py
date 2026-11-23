"""Checks of the paper's statements (uniform start table only).

  python paper_checks.py appendix [--table T] [--coeffs FILE] [--real] [--json OUT]
      Appendix A: the IPF limit and adjustment factors, the polynomial E(x; n) of cell 111, and the checks reported
      there (root, end coefficients as products of the 103 boundary forms, forced factors of every coefficient,
      real roots, real solutions). Without --coeffs, E is computed with msolve (about a minute).

  python paper_checks.py line [--t T ...]
      The full eliminant on the line n_0 + t d of Section 3 (Remark "full line"): data/line holds the Newton
      interpolants in t, modulo two primes, of r_k = c_k F_511 / (c_511 F_k) for every k, from 759 exact eliminants
      (0 <= t <= 760). For each T (default 1000, far outside the interpolation range) the script computes E(x; n_0 + T d)
      with msolve and checks all 512 coefficients against the interpolants.

  python paper_checks.py crossings [--max M] [--type T]
      The isolated crossings of the line with the 325 boundary hyperplanes (data/line/crossings.json). At each
      crossing t, E(x; n_0 + t d) from the interpolants is divided, modulo both primes, by the eliminant G computed with
      msolve at the crossing: the quotient must be c x^e on a constant-side form (e > 0), a constant with
      deg G = 511 - |e| on a lead-side form (e < 0; 495 for the total N, where 16 branches escape), and of degree
      |e| of the type on a neutral form (e = 0)."""
import argparse, hashlib, json, time
from fractions import Fraction
from functools import reduce
from math import gcd, prod
from pathlib import Path
import mpmath as mp
import flint
from ipf333 import ipf, factorization_residual, eliminant, form_values
from ipf333.forced import forced
mp.mp.dps = 60

BASE = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3, 8]
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("check", choices=["appendix", "line", "crossings"])
ap.add_argument("--max", type=int, help="crossings: check only the first M")
ap.add_argument("--type", help="crossings: check only crossings of this type, e.g. N or B233")
ap.add_argument("--t", type=int, nargs="+", default=[1000], help="line: values of t to check")
ap.add_argument("--table", default=",".join(map(str, BASE))); ap.add_argument("--coeffs"); ap.add_argument("--msolve")
ap.add_argument("--json"); ap.add_argument("--write-coeffs")
ap.add_argument("--real", action="store_true", help="also compute all real solutions of the scaling equations (msolve)")
a = ap.parse_args(); n = [int(x) for x in a.table.split(",")]; assert len(n) == 27

if a.check in ("line", "crossings"):
    import sys
    import numpy as np
    from ipf333.forms import form_values
    from ipf333.forced import order
    D = Path(__file__).resolve().parent / "data" / "line"; L = json.loads((D / "line.json").read_text())
    NW = {q: np.load(D / f"newton_{q}.npy") for q in L["primes"]}; nodes = L["nodes"]

    def r_at(q, t):  # r_k(t) mod q, t an integer or a Fraction
        tq = t.numerator * pow(t.denominator, q - 2, q) % q if isinstance(t, Fraction) else t % q
        A_ = NW[q]; r = [int(x) for x in A_[-1]]
        for j in range(len(nodes) - 2, -1, -1):
            f = (tq - nodes[j]) % q; r = [(x * f + int(y)) % q for x, y in zip(r, A_[j])]
        return r

    def E_line(q, t):  # coefficients of E(x; n_0 + t d) mod q, scaled so that c_511 = F_511
        m = [Fraction(b) + t * v for b, v in zip(L["base"], L["v"])]
        fv = [(c, e, v.numerator * pow(v.denominator, q - 2, q) % q) for c, e, v in form_values(m)]
        r = r_at(q, t)
        return [r[k] * prod(pow(v, order(c, e, k), q) for c, e, v in fv) % q for k in range(512)]

    allok = True
    if a.check == "line":
        for t in a.t:
            m = [b + t * v for b, v in zip(L["base"], L["v"])]
            t0 = time.time(); co = eliminant(m, msolve=a.msolve, threads=16); sec = time.time() - t0
            Fk = [forced(m, k) for k in range(512)]
            for q in L["primes"]:
                r = r_at(q, t)
                pred = [r[k] * co[511] % q * Fk[k] % q for k in range(512)]
                got = [co[k] % q * Fk[511] % q for k in range(512)]
                bad = sum(x != y for x, y in zip(pred, got)); allok &= bad == 0
                print(f"t = {t}, p = {q}: {512 - bad}/512 coefficients match the interpolants (msolve {sec:.0f} s)")
    else:
        X = json.loads((D / "crossings.json").read_text())[: a.max]
        if a.type: X = [c for c in X if c["type"] == a.type]
        for c in X:
            t = Fraction(*c["t"]); e = c["e"]
            m = [b * t.denominator + t.numerator * v for b, v in zip(L["base"], L["v"])]  # integer multiple of n_0 + t d
            assert sum(l * x for l, x in zip(c["lam"], m)) == 0
            t0 = time.time(); G = eliminant(m, msolve=a.msolve, threads=16); sec = time.time() - t0
            res = []
            for q in L["primes"]:
                El = flint.nmod_poly(E_line(q, t), q)
                Gq = flint.nmod_poly([g * pow(t.denominator, i, q) % q for i, g in enumerate(G)], q)  # G(den x)
                Qt, Rm = divmod(El, Gq); qd = Qt.degree()
                if e > 0: ok = Rm == 0 and qd == e and all(Qt[i] == 0 for i in range(e))
                elif e < 0: ok = Rm == 0 and qd == 0 and len(G) - 1 == (495 if c["type"] == "N" else 511 + e)  # on N = 0 only 16 branches escape
                else: ok = Rm == 0 and qd == 511 - (len(G) - 1)
                res.append(ok); allok &= ok
            print(f"{c['type']:>5} e = {e:>4}  t = {str(t):>7}: deg G = {len(G) - 1}, quotient degree {qd}, "
                  f"{'as predicted' if all(res) else 'NOT as predicted'} mod both primes (msolve {sec:.0f} s)")
    sys.exit(0 if allok else 1)
R = {"table": n}
lab = lambda i, j, k: f"{i+1}{j+1}{k+1}"

# 1. IPF limit and raking factors
p, (A, B, C), sweeps = ipf(n, dps=60); resid = factorization_residual(p, (A, B, C))
R["ipf_sweeps"] = sweeps
R["fitted"] = {lab(*c): mp.nstr(v, 25) for c, v in p.items()}
R["raking"] = {"a": {f"{i+1}{j+1}": mp.nstr(v, 25) for (i, j), v in A.items()},
               "b": {f"{i+1}{k+1}": mp.nstr(v, 25) for (i, k), v in B.items()},
               "c": {f"{j+1}{k+1}": mp.nstr(v, 25) for (j, k), v in C.items()}, "max_rel_resid": mp.nstr(resid, 3)}
x = p[(0, 0, 0)]

# 2. the eliminant
if a.coeffs:
    co = [int(l) for l in Path(a.coeffs).read_text().split()]
else:
    t0 = time.time(); co = eliminant(n, msolve=a.msolve, threads=16); R["msolve_seconds"] = round(time.time() - t0)
txt = "\n".join(map(str, co)) + "\n"
if a.write_coeffs: Path(a.write_coeffs).write_text(txt)
R["degree"] = len(co) - 1; R["sha256"] = hashlib.sha256(txt.encode()).hexdigest()
dig = [len(str(abs(c))) for c in co]; R["digits"] = {"min": min(dig), "max": max(dig), "total": sum(dig)}
R["primitive"] = reduce(gcd, co) == 1 and co[-1] > 0

# 3. the IPF limit is a root: Newton distance |E/E'| at 60 digits
with mp.workdps(1500):
    X = mp.mpf(x); E = mp.polyval([mp.mpf(c) for c in co[::-1]], X)
    dE = mp.polyval([mp.mpf(k * c) for k, c in list(enumerate(co))[1:][::-1]], X)
    R["newton_distance"] = mp.nstr(abs(E / dE), 3); R["p111"] = mp.nstr(x, 40)

# 4. end coefficients: product of the 511 roots = prod_{e>0} l^e / prod_{e<0} l^|e|, exactly
vals = form_values(n)
num = prod(v ** e for _, e, v in vals if e > 0); den = prod(v ** -e for _, e, v in vals if e < 0)
R["root_product_identity"] = Fraction(-co[0], co[-1]) == Fraction(num, den)  # (-1)^511 c_0 / c_511
D = Fraction(den, co[-1]); R["content_D"] = str(D)                            # E_true = D * E_primitive
R["c0_equals_minus_D_prod"] = Fraction(co[0]) * D == -num
R["forms"] = [{"class": c, "exponent": e, "value": v} for c, e, v in vals]

# 5. forced factors of every coefficient: D c_k divisible by F_k(n) for all k
ok = 0
for k in range(len(co)):
    q = Fraction(co[k]) * D / forced(n, k)
    ok += q.denominator == 1
R["forced_divisibility_all_k"] = f"{ok}/{len(co)}"
R["Q1_digits"] = len(str(abs(int(Fraction(co[1]) * D / forced(n, 1)))))

# 6. real roots (certified isolation with Arb) and the root reached by IPF
t0 = time.time(); roots = flint.fmpz_poly(co).complex_roots()
real = sorted(float(r.real.mid()) for r, m in roots if r.imag == 0)
R["real_roots"] = len(real); R["positive_real_roots"] = sum(r > 0 for r in real)
R["real_root_values"] = [round(r, 6) for r in real]; R["roots_seconds"] = round(time.time() - t0)
R["ipf_root_index_among_positive"] = sorted(r for r in real if r > 0).index(min((r for r in real if r > 0), key=lambda r: abs(r - float(x))))
# 7. real solutions of the full system: exactly one is a positive table (the IPF limit)
if a.real:
    from ipf333.realsol import real_solutions, sign_summary
    sols = real_solutions(n, msolve=a.msolve, threads=16); nr, npos = sign_summary(sols)
    pos = [s_ for s_ in sols if all(lo > 0 for lo, hi in s_)]
    R["real_solutions"] = nr; R["positive_tables"] = npos
    lo, hi = pos[0][-1] if len(pos) == 1 else (1, 0)
    mid = (mp.mpf(lo.numerator) / lo.denominator + mp.mpf(hi.numerator) / hi.denominator) / 2
    R["positive_table_p111_minus_ipf"] = mp.nstr(abs(mid - x), 3)   # IPF value is accurate to about 55 digits
print(json.dumps({k: v for k, v in R.items() if k not in ("forms", "fitted", "raking", "real_root_values")}, indent=1))
if a.json: Path(a.json).write_text(json.dumps(R, indent=1, default=str))
