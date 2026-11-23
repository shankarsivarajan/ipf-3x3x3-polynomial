"""The 103 linear forms in the end coefficients of the 3x3x3 eliminant (Table 1 of the paper) and their exponents."""
import json
from math import prod
from pathlib import Path

FORMS = json.load(open(Path(__file__).with_name("forms.json")))["forms"]


def form_values(n):
    """[(class, exponent, value)] for the 103 forms on the table n."""
    return [(f["class"], f["exponent"], sum(l * x for l, x in zip(f["lambda"], n))) for f in FORMS]


def end_factorization(n, co):
    """Exact quotients c_0 / prod_{e>0} l^e and c_d / prod_{e<0} l^|e| (Conjecture 3.1: both are constants
    independent of n up to sign; for the uniform seed both equal +-1 times a fixed integer content)."""
    vals = form_values(n)
    pc = prod(v ** e for _, e, v in vals if e > 0); pl = prod(v ** -e for _, e, v in vals if e < 0)
    from fractions import Fraction
    return Fraction(co[0], pc), Fraction(co[-1], pl), sum(e for _, e, _ in vals if e > 0), sum(-e for _, e, _ in vals if e < 0)
