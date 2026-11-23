# The polynomial behind the IPF limit of a 3x3x3 table

Companion code for

> Shankar Sivarajan, *The boundary structure of the maximum likelihood eliminant of 3x3x3 tables* (2026).

## Overview

Iterative proportional fitting (IPF, also called raking) repeatedly rescales the uniform $3\times 3\times 3$ table until its three two-way margins match those of an observed table $n$. In each sweep, the algorithm matches $n_{ij+}$, then $n_{i+k}$, and then $n_{+jk}$, and returns the cumulative adjustment factors. The limit $\hat p$ is the unique positive table of the form $\hat p_{ijk} = a_{ij} b_{ik} c_{jk}$ with the two-way margins of $n$, and is the maximum likelihood estimate under the no-three-way interaction model. Although the limit $\hat p$ does not depend on sweep order, the adjustment factors are determined only up to the gauge transformation $a_{ij} \to a_{ij} u_i v_j$, $b_{ik} \to b_{ik} w_k / u_i$, $c_{jk} \to c_{jk} / (v_j w_k)$, so another implementation may report different, equally valid factors.

In general $\hat p$ has no closed form, but each fitted cell is a root of a univariate integer polynomial $E(x; n)$ of degree $511$ whose roots give that cell's value across all complex solutions of the scaling equations. This package provides:

* **IPF solver (`ipf`):** Computes the IPF limit to 55 digits together with the cumulative adjustment factors $(a, b, c)$.
* **Exact eliminants (`eliminant`):** Builds the scaling equations for [msolve](https://msolve.lip6.fr) and computes the primitive integer polynomial $E(x; n)$ of any chosen cell. It handles variable reordering in `msolve -P 1` output and reconstructs the target cell's univariate polynomial from `msolve`'s rational univariate representation (RUR) when `msolve` selects a non-coordinate separating linear form. Both `ipf` and `eliminant` also accept an arbitrary positive start table (`start=` and `seed=`), replacing the model with $t_{ijk} a_{ij} b_{ik} c_{jk}$, for which the generic eliminant has degree 2,457.
* **Real solutions (`real_solutions`, `sign_summary`):** Isolates all real solutions of the scaling equations via `msolve` and classifies their signs, verifying that exactly one is a positive table (the IPF limit).

The eliminant gives an exact algebraic description of the IPF limit and all complex solutions, not a faster numerical replacement for IPF. Starting from the uniform table, computing a degree-511 eliminant typically takes about a minute on 8 to 16 threads and uses under 1 GB of memory; generic start tables are substantially heavier.

## Quickstart

From Python, `ipf` computes the fitted table $\hat p$ and `eliminant` computes the exact integer polynomial $E(x; n) \in \mathbb{Z}[x]$ for any chosen cell:

```python
from ipf333 import ipf, eliminant

n = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3, 8]
p_hat, (a, b, c), sweeps = ipf(n)        # dicts keyed by 0-based indices
E = eliminant(n, cell=(0, 0, 0))         # [c_0, ..., c_511], constant term first
```

From the command line, `ipf_polynomial.py` accepts 27 nonnegative integers separated by commas or whitespace (either inline or in a file), ordered lexicographically as $111, 112, 113, 121, \ldots, 333$ with the last index varying fastest:

```bash
python ipf_polynomial.py 3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8 --out E.txt
python ipf_polynomial.py table.txt --cell 213 --real
```

Output coefficients are written one per line, constant term first, normalized to be coprime with a positive leading coefficient. Passing `--real` also isolates all real solutions of the scaling equations and checks how many are positive tables.

## Paper-to-code guide

| Paper | Entry point / file | What it computes |
|---|---|---|
| Section 1 | `ipf`, `factorization_residual` (`ipf333/ipf.py`) | Arbitrary-precision IPF limit and adjustment factors $(a, b, c)$ |
| Section 1, Proposition 7.3 | `eliminant`, `msolve_input` (`ipf333/eliminant.py`), `ipf_polynomial.py` | Exact primitive eliminant $E(x; n) \in \mathbb{Z}[x]$ for uniform or arbitrary positive start tables (`start=`, `seed=`) |
| Section 2, Table 1 | `FORMS`, `form_values` (`ipf333/forms.py`, `ipf333/forms.json`) | The 103 active boundary forms and their exponents |
| Conjecture 3.1 | `end_factorization` (`ipf333/forms.py`) | Checks factorization of $c_0$ and $c_{511}$ into the 103 boundary forms |
| Proposition 6.1 | `forced`, `order` (`ipf333/forced.py`) | Forced factor $F_k(n)$ dividing each coefficient $c_k(n)$ |
| Remark 6.16 | `data/line/` | Line of tables, mod-prime interpolants of all 512 coefficients, and 325 hyperplane crossings |
| Appendix A | `real_solutions`, `sign_summary` (`ipf333/realsol.py`), `data/base_table_E.txt`, `data/base_table_example.json` | Certified real-solution isolation, positivity counts, and Appendix A reference values |
| Sections 3, 6, Appendix A | `paper_checks.py` | Verification suite for the uniform-start statements |

## Reproducing the paper

Requirements: Python 3.10+, `mpmath`, `python-flint`, `numpy`, and the [`msolve`](https://msolve.lip6.fr) binary (tested with version 0.11.0, either on `PATH` or passed with `--msolve`).

```bash
pip install -e .                                        # or: pip install mpmath python-flint numpy
python paper_checks.py appendix --real --json out.json
python paper_checks.py line --t 1000 -45
python paper_checks.py crossings
```

The three subcommands check the paper's uniform-start results:

* **`appendix`:** Reproduces the worked example of Appendix A, verifies that the end coefficients factor into the 103 boundary forms of Table 1 (Conjecture 3.1) and that every coefficient is divisible by its forced factor (Proposition 6.1), and counts the real roots and real solutions.
* **`line`:** Checks the full eliminant along the line of tables in Section 3. The `data/line` directory stores, modulo two primes, the polynomial interpolants in $t$ for all 512 coefficients (after dividing out their forced factors), reconstructed from 759 exact eliminants with $0 \le t \le 760$. For each requested value of $t$, the script computes the eliminant with `msolve` and compares all 512 coefficients against the interpolants.
* **`crossings`:** Checks factorization at the isolated crossings of the Section 3 line with the 325 boundary hyperplanes (`data/line/crossings.json`). At each crossing, the interpolated polynomial is divided by the `msolve` eliminant evaluated at that point; the quotient must be $c x^e$, a constant, or a polynomial of degree $|e|$ for constant-side, lead-side, and neutral forms, respectively. Most crossings finish in under a minute, though margin crossings take longer.

## Related

The case where one variable is binary ($2\times J\times K$ tables) is treated in the companion repository [`ipf-binary-eliminant`](https://github.com/shankarsivarajan/ipf-binary-eliminant) and the paper

> Shankar Sivarajan, *An explicit determinantal formula for limits of iterative proportional fitting with a binary variable*, [hexagon:2611.00705](https://hexagonmath.org/2611.00705) (2026).

The case of one-dimensional margins (matrix and tensor scaling) is treated in [`sinkhorn-eliminant`](https://github.com/shankarsivarajan/sinkhorn-eliminant) and the paper

> Shankar Sivarajan, *An explicit determinantal formula for Sinkhorn limits of matrices and tensors*, [hexagon:2611.00433](https://hexagonmath.org/2611.00433) (2026).
