"""
Floating-point parameter scan for
    N = n1 + n2,  n1, n2 sums of two squares,  Omega(n1 n2) <= k,     N = 2 mod 4.
Heuristic main-term model (NOT a proof, NOT certified).

Setting: n <= N, n = 1 mod 4 (then N - n = 1 mod 4).  X = N/4.  For a prime q not dividing N,
q | n and q | N-n each have density 1/q, and never both.  "Inert" = q = 3 mod 4, "split" = q = 1 mod 4.
The level of distribution is 1 (trivial remainder), shared by all sieve components: sum of theta_i <= 1.

Each side X in {A: n, B: N-n} must have no inert factor (=> sum of two squares, given X = 1 mod 4) and no
prime factor < z = N^alpha.  Inert primes are sifted on side X up to y_X = N^(1/2 - delta_X).  If
delta_X = 0 the side is clean; otherwise survivors on side X may have two inert factors >= y_X
("bad"), which are removed by a cheap switching bound (Section 'bad' below).

Vector sieve (Brudern-Fouvry, k components):  prod a_i >= sum_i a_i^- prod_{j!=i} a_j^+ - (k-1) prod a_j^+.

Schemes
  S1  2 linear components: inert on n(N-n) up to y (both sides, delta), split on n(N-n) up to z.
  S2  3 components: A-inert semi-linear up to sqrt N (A clean), B-inert semi-linear up to N^(1/2-delta),
      split on n(N-n) linear up to z.
  S3  4 semi-linear components: A-inert, B-inert (both delta), A-split, B-split.

Normalisation: everything is divided by X V_A V_B (times the common singular-series correction).
  main     = L - lam R
  R        = int_alpha^beta1 (1 - u/beta1) U(1-u) du/u      (2 sides x density 1/2)
  bad_X    = 2 c(delta_X, alpha) e^gamma sqrt(alpha (1/2 - delta_X)) U_other
             (upper sieve on {N - t}, t = q1 q2 m, for the other side's conditions at level 1/2)
Richert weight on both sides jointly:  w = 1 - lam sum_{q | n1 n2, q split, z <= q < N^beta1} (1 - log q/log y1)
  => Omega(n1 n2) < 1/lam + 2/beta1,  and always Omega(n1 n2) < 2/alpha.
"""
import itertools
import numpy as np
from sieve_core import F_half, f_half, F_lin, f_lin, c_switch, EG

_simplex_cache = {}


def simplex(k, step):
    key = (k, step)
    if key not in _simplex_cache:
        m = int(round(1 / step))
        pts = [c for c in itertools.product(range(1, m), repeat=k - 1) if sum(c) < m]
        P = np.array(pts, dtype=float) * step
        _simplex_cache[key] = np.hstack([P, 1 - P.sum(axis=1, keepdims=True)])
    return _simplex_cache[key]


def funcs(kind):
    return (F_half, f_half) if kind == "half" else (F_lin, f_lin)


def vec_lower(comps, total, step):
    W = simplex(len(comps), step) * total
    Fs, fs = [], []
    for (kind, e), th in zip(comps, W.T):
        F, f = funcs(kind)
        Fs.append(F(th / e)); fs.append(f(th / e))
    Fs = np.array(Fs); fs = np.array(fs)
    prodF = np.prod(Fs, axis=0)
    L = sum(fs[i] * prodF / Fs[i] for i in range(len(comps))) - (len(comps) - 1) * prodF
    return float(L.max())


def vec_upper(comps, total, step):
    W = simplex(len(comps), step) * total
    prod = np.ones(W.shape[0])
    for (kind, e), th in zip(comps, W.T):
        prod *= funcs(kind)[0](th / e)
    return float(prod.min())


def side_upper(delta, alpha):
    # one side's conditions (inert up to N^(1/2-delta), all primes up to z) at level 1/2
    return vec_upper([("half", 0.5 - delta), ("half", alpha)], 0.5, 0.005)


def scheme(name, delta, alpha):
    if name == "S1":
        comps = [("lin", 0.5 - delta), ("lin", alpha)]
        bad = 2 * (2 * c_switch(delta, alpha) * EG * np.sqrt(alpha * (0.5 - delta)) * side_upper(delta, alpha))
    elif name == "S2":
        comps = [("half", 0.5), ("half", 0.5 - delta), ("lin", alpha)]
        bad = 2 * c_switch(delta, alpha) * EG * np.sqrt(alpha * (0.5 - delta)) * side_upper(0.0, alpha)
    else:
        comps = [("half", 0.5 - delta), ("half", 0.5 - delta), ("half", alpha), ("half", alpha)]
        bad = 2 * (2 * c_switch(delta, alpha) * EG * np.sqrt(alpha * (0.5 - delta)) * side_upper(delta, alpha))
    return comps, bad


def evaluate(name, delta, alpha, betas=np.arange(0.10, 0.51, 0.02)):
    comps, bad = scheme(name, delta, alpha)
    stepL = 0.005 if len(comps) <= 3 else 0.01
    L = vec_lower(comps, 1.0, stepL)
    slack = L - bad
    if slack <= 0:
        return []
    out = [(int(np.ceil(2 / alpha - 1e-9)) - 1, 0.0, None, slack)]          # no weights
    u = np.linspace(alpha, 0.5, 40)
    stepU = 0.02 if len(comps) <= 3 else 0.04
    U = np.array([vec_upper(comps, 1 - ui, stepU) for ui in u])
    for b1 in betas:
        if b1 <= alpha:
            continue
        m = u <= b1
        if m.sum() < 2:
            continue
        R = np.trapz((1 - u[m] / b1) * U[m] / u[m], u[m])
        lam = 0.95 * slack / R
        bound = min(2 / alpha, 1 / lam + 2 / b1)
        out.append((int(np.ceil(bound - 1e-9)) - 1, lam, b1, 0.05 * slack))
    return out


if __name__ == "__main__":
    import sys
    deltas = np.arange(0.02, 0.245, 0.02)
    alphas = np.arange(0.004, 0.06, 0.004)
    for name in ("S1", "S2", "S3"):
        best = {}
        for d in deltas:
            for a in alphas:
                for k, lam, b1, net in evaluate(name, d, a):
                    key = (k, lam > 0)
                    if key not in best or net > best[key][-1]:
                        best[key] = (d, a, lam, b1, net)
        for weighted in (False, True):
            ks = sorted(k for (k, w) in best if w == weighted)
            if not ks:
                print(f"{name} {'Richert' if weighted else 'no weights':10s}: no positive configuration")
                continue
            k = ks[0]; d, a, lam, b1, net = best[(k, weighted)]
            print(f"{name} {'Richert' if weighted else 'no weights':10s}: smallest Omega(n1 n2) bound k = {k:3d}"
                  f"  (delta={d:.2f}, alpha={a:.3f}, lam={lam:.3f}, beta1={b1}, kept={net:+.4f})")
        sys.stdout.flush()
