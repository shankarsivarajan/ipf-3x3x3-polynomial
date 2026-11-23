"""Forced factor F_k(n) of the coefficient c_k (Proposition 6.1 of the paper)."""
from math import prod
from .forms import form_values


def order(cls, e, k):
    if cls == "N":
        m = 511 - k
        return 17 if m == 0 else max(0, 16 - m)
    return max(0, e - k) if e > 0 else max(0, -e - (511 - k))


def forced(n, k):
    return prod(v ** order(c, e, k) for c, e, v in form_values(n))
