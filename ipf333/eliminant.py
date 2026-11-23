"""Exact eliminant E(x; n) of the fitted cell p_111 over all complex critical points, via msolve.
Requires the msolve binary (https://msolve.lip6.fr); pass its path or put it on PATH."""
import ast, itertools, re, shutil, subprocess, tempfile
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from pathlib import Path


def _rat_recon(a, m):
    r0, r1 = m, a % m
    t0, t1 = 0, 1
    lim = 1 << ((m // 2).bit_length() // 2)
    while r1 >= lim:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    if t1 < 0:
        r1, t1 = -r1, -t1
    return Fraction(r1, t1)


def _charpoly_mod_p(P_co, D_co, N_co, c_mul, q):
    import flint
    d = len(P_co) - 1
    P = flint.nmod_poly(P_co, q)
    D = flint.nmod_poly(D_co, q)
    N = flint.nmod_poly(N_co, q)
    g, inv_cD, _ = (c_mul * D).xgcd(P)
    if g != 1:
        raise ZeroDivisionError
    alpha = (-N * inv_cD) % P
    lc_inv = pow(int(P[d]), q - 2, q)
    a = [(int(P[i]) * lc_inv) % q for i in range(d)]
    tr_t = [d % q]
    for k in range(1, d):
        val = (k * a[d - k] + sum(a[d - k + i] * tr_t[i] for i in range(1, k))) % q
        tr_t.append((-val) % q)
    psums = [d % q]
    cur = flint.nmod_poly([1], q)
    for _ in range(1, d + 1):
        cur = (cur * alpha) % P
        psums.append(sum(int(cur[j]) * tr_t[j] for j in range(min(len(cur), d))) % q)
    e = [1]
    for k in range(1, d + 1):
        val = sum(((1 if (i - 1) % 2 == 0 else -1) * e[k - i] * psums[i]) for i in range(1, k + 1)) % q
        e.append((val * pow(k, q - 2, q)) % q)
    return [((1 if (d - k) % 2 == 0 else -1) * e[d - k]) % q for k in range(d + 1)]


def _rur_to_eliminant(P_co, D_co, N_co, c_mul):
    """Exact primitive integer eliminant of x = -N(t) / (c_mul * D(t)) mod P(t) via modular traces and CRT."""
    import flint
    d = len(P_co) - 1
    p_cand = (1 << 61) - 1
    M = 1
    acc = [0] * (d + 1)
    prev_probe = None
    primes_used = 0
    while True:
        while not flint.fmpz(p_cand).is_prime():
            p_cand -= 2
        q = p_cand
        p_cand -= 2
        if int(P_co[-1]) % q == 0 or c_mul % q == 0:
            continue
        try:
            cp = _charpoly_mod_p(P_co, D_co, N_co, c_mul, q)
        except ZeroDivisionError:
            continue
        if M == 1:
            acc, M = cp, q
        else:
            inv_M = pow(M % q, q - 2, q)
            acc = [acc[i] + M * (((cp[i] - acc[i]) * inv_M) % q) for i in range(d + 1)]
            M *= q
        primes_used += 1
        if primes_used >= 20 and primes_used % 10 == 0:
            probe = (_rat_recon(acc[0], M), _rat_recon(acc[d - 1], M), _rat_recon(acc[d // 2], M))
            if probe == prev_probe:
                rats = [_rat_recon(v, M) for v in acc]
                den = 1
                for r in rats:
                    den = lcm(den, r.denominator)
                co = [int(r * den) for r in rats]
                g = reduce(gcd, co)
                co = [c // g for c in co]
                if co[-1] < 0:
                    co = [-c for c in co]
                # Verify against a fresh prime to guard against premature probe stabilization.
                while True:
                    while not flint.fmpz(p_cand).is_prime():
                        p_cand -= 2
                    qf = p_cand
                    p_cand -= 2
                    if int(P_co[-1]) % qf == 0 or c_mul % qf == 0 or co[-1] % qf == 0:
                        continue
                    try:
                        cp_chk = _charpoly_mod_p(P_co, D_co, N_co, c_mul, qf)
                        lc_inv = pow(co[-1] % qf, qf - 2, qf)
                        if [(c * lc_inv) % qf for c in co] == cp_chk:
                            return co
                        break
                    except ZeroDivisionError:
                        continue
            prev_probe = probe


def msolve_input(n, seed=None):
    """Polynomial system for the 3x3x3 scaling equations p_ijk = t_ijk a_ij b_ik c_jk with the two-way margins of n.
    Gauge a_1j = a_i1 = b_1k = 1, so c_11 = p_111 / t_111 is the last variable (eliminant() rescales by t_111). Returns msolve input text."""
    cells = list(itertools.product(range(3), repeat=3)); N = dict(zip(cells, n))
    T = {c: 1 for c in cells} if seed is None else dict(zip(cells, seed))
    A = {(i, j): f"a{i+1}{j+1}" for i in range(1, 3) for j in range(1, 3)}
    B = {(i, k): f"b{i+1}{k+1}" for i in range(1, 3) for k in range(3)}
    C = {(j, k): f"c{j+1}{k+1}" for j in range(3) for k in range(3)}

    def P(c):
        i, j, k = c; f = [str(T[c])] if T[c] != 1 else []
        if (i, j) in A: f.append(A[(i, j)])
        if (i, k) in B: f.append(B[(i, k)])
        f.append(C[(j, k)]); return "*".join(f)
    rhs = lambda v: f"-{v}" if v >= 0 else f"+{-v}"
    eqs = []
    for (i, j) in A:
        eqs.append("+".join(P(c) for c in cells if c[:2] == (i, j)) + rhs(sum(N[c] for c in cells if c[:2] == (i, j))))
    for (i, k) in B:
        sel = [c for c in cells if (c[0], c[2]) == (i, k)]; eqs.append("+".join(map(P, sel)) + rhs(sum(N[c] for c in sel)))
    for (j, k) in C:
        sel = [c for c in cells if c[1:] == (j, k)]; eqs.append("+".join(map(P, sel)) + rhs(sum(N[c] for c in sel)))
    vs = list(A.values()) + list(B.values()) + [v for key, v in C.items() if key != (0, 0)] + [C[(0, 0)]]
    return ",".join(vs) + "\n0\n" + ",\n".join(eqs) + "\n"


def move_cell(v, cell):
    """Reorder a 27-vector so that the 0-based cell (i, j, k) moves to (0, 0, 0) (swap level 0 with level i, etc.).
    The swap is an involution, and the fitted value of `cell` for v is the fitted value of 111 for the result."""
    sw = lambda x, t: t if x == 0 else (0 if x == t else x)
    cells = list(itertools.product(range(3), repeat=3))
    return [v[9 * sw(i, cell[0]) + 3 * sw(j, cell[1]) + sw(k, cell[2])] for (i, j, k) in cells]


def eliminant(n, seed=None, msolve=None, threads=8, cell=(0, 0, 0)):
    """Integer coefficients [c_0, ..., c_d] of the eliminant of the fitted value of `cell` (0-based; default 111)
    over all complex solutions, primitive with a positive leading coefficient. seed = start table (positive integers)."""
    if cell != (0, 0, 0):
        n = move_cell(n, cell); seed = None if seed is None else move_cell(seed, cell)
    exe = msolve or shutil.which("msolve")
    if not exe:
        raise FileNotFoundError("msolve not found")
    with tempfile.TemporaryDirectory() as d:
        fi, fo = Path(d) / "in.ms", Path(d) / "out.ms"
        fi.write_text(msolve_input(n, seed))
        subprocess.run([exe, "-f", str(fi), "-P", "1", "-o", str(fo), "-t", str(threads)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        txt = fo.read_text()
    m_hdr = re.search(r"(\[\s*'[^']*'(?:\s*,\s*'[^']*')*\s*\])\s*,\s*\[([-0-9,\s]*)\]\s*,\s*\[1,", txt)
    vars_out = ast.literal_eval(m_hdr.group(1))
    lin = [int(v) for v in m_hdr.group(2).split(",")]
    if vars_out[-1] == "c11" and not any(lin[:-1]) and lin[-1] == 1:
        m = re.compile(r"\[\[(\d+),\s*\[([^\]]*)\]\]").search(txt, txt.index("[1,"))
        co = [int(x) for x in m.group(2).replace("\n", "").split(",")]
    else:
        idx_rur = txt.index("[1,", m_hdr.start(2))
        m_sols = re.search(r"\]\]\]\]\s*,\s*\[1,", txt[idx_rur:])
        rur_end = idx_rur + m_sols.start() + 3 if m_sols else txt.rindex("]]]") + 3
        rur_data = ast.literal_eval(txt[idx_rur:rur_end])
        P_co, D_co, params = rur_data[1][0][1], rur_data[1][1][1], rur_data[1][2]
        if vars_out[-1] != "c11":
            idx = vars_out.index("c11")
            co = _rur_to_eliminant(P_co, D_co, params[idx][0][1], params[idx][1])
        else:
            if lin[-1] == 0:
                raise RuntimeError("c_11 cannot be recovered from separating linear form")
            # t = sum(lin[i] * v_i) + lin[-1] * c_11, with v_i = -N_i / (c_i * D).
            c_lcm = 1
            for i, coef in enumerate(lin[:-1]):
                if coef:
                    c_lcm = lcm(c_lcm, params[i][1])
            deg_n = max(len(D_co) + 1, *(len(params[i][0][1]) for i, coef in enumerate(lin[:-1]) if coef))
            N_co = [0] * deg_n
            # c_11 = -(c_lcm * t * D + sum lin[i]*(c_lcm/c_i)*N_i) / (-lin[-1]*c_lcm * D)
            for k, dc in enumerate(D_co):
                N_co[k + 1] -= c_lcm * dc
            for i, coef in enumerate(lin[:-1]):
                if coef:
                    scale = coef * (c_lcm // params[i][1])
                    for k, nc in enumerate(params[i][0][1]):
                        N_co[k] -= scale * nc
            while len(N_co) > 1 and N_co[-1] == 0:
                N_co.pop()
            co = _rur_to_eliminant(P_co, D_co, N_co, lin[-1] * c_lcm)
    t = 1 if seed is None else seed[0]
    if t != 1:   # roots in c_11 = p_111 / t_111; substitute c_11 = x / t and clear denominators
        d = len(co) - 1; co = [c * t ** (d - k) for k, c in enumerate(co)]
    g = reduce(gcd, co); co = [c // g for c in co]
    if co[-1] < 0:
        co = [-c for c in co]
    return co
