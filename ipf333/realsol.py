"""Real solutions of the scaling equations (msolve's certified isolating boxes) and their sign patterns."""
import re, shutil, subprocess, tempfile
from fractions import Fraction
from pathlib import Path
from .eliminant import msolve_input


def real_solutions(n, msolve=None, threads=8, seed=None):
    """List of real solutions, each a list of (lo, hi) Fraction intervals in the variable order of msolve_input."""
    exe = msolve or shutil.which("msolve")
    with tempfile.TemporaryDirectory() as d:
        fi, fo = Path(d) / "in.ms", Path(d) / "out.ms"
        fi.write_text(msolve_input(n, seed))
        subprocess.run([exe, "-f", str(fi), "-o", str(fo), "-t", str(threads)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        txt = fo.read_text()
    num = r"(-?\d+)(?:\s*/\s*2\^(\d+))?"
    val = lambda m: Fraction(int(m[0]), 2 ** int(m[1]) if m[1] else 1)
    body = txt[txt.index("[[[") + 1:]
    sols = []
    for sol in re.findall(r"\[(\[[^\[\]]*\](?:,\s*\[[^\[\]]*\])*)\]", body):
        ivs = re.findall(r"\[\s*" + num + r"\s*,\s*" + num + r"\s*\]", sol)
        sols.append([(val(iv[:2]), val(iv[2:])) for iv in ivs])
    return sols


def sign_summary(sols):
    """Number of real solutions, and how many have every coordinate positive (that is, a positive table)."""
    pos = [s for s in sols if all(lo > 0 for lo, hi in s)]
    return len(sols), len(pos)
