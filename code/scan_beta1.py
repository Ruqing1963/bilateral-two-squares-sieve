"""
Scheme S1 with the Richert range extended beyond beta1 = 1/2 (allowed here because the level of
distribution is 1, so A_q with q = N^u still has level 1 - u > 0 for u < 1).

For w_n = 1 - lam (T_A + T_B), T_X = sum_{q | n_X, q split, z <= q < N^beta1} (1 - log q/log y1) >= 0:
   w_n > 0  =>  Omega(n1 n2) < 1/lam + 2/beta1        (joint bound, as in the paper)
            =>  Omega(n_X)   < 1/lam + 1/beta1        (each side, since T_other >= 0)
Floating point, heuristic main-term model as in scan.py.
"""
import numpy as np
from scan import scheme, vec_lower, vec_upper


def scan(deltas, alphas, betas):
    best_joint, best_side = {}, {}
    for d in deltas:
        for a in alphas:
            comps, bad = scheme("S1", d, a)
            L = vec_lower(comps, 1.0, 0.005)
            slack = L - bad
            if slack <= 0:
                continue
            u = np.linspace(a, 0.97, 120)
            U = np.array([vec_upper(comps, 1 - x, 0.02) for x in u])
            for b1 in betas:
                m = u <= b1
                R = np.trapz((1 - u[m] / b1) * U[m] / u[m], u[m])
                lam = 0.95 * slack / R
                kj = int(np.ceil(1 / lam + 2 / b1 - 1e-9)) - 1
                ks = min(int(np.ceil(1 / lam + 1 / b1 - 1e-9)) - 1, kj - 1)
                rec = (d, a, round(lam, 4), round(b1, 2))
                best_joint.setdefault(kj, rec); best_side.setdefault(ks, rec)
    return best_joint, best_side


if __name__ == "__main__":
    deltas = [0.18, 0.20, 0.22, 0.24, 0.245]
    alphas = np.arange(0.006, 0.04, 0.002)
    for label, betas in (("beta1 <= 0.5", np.arange(0.30, 0.501, 0.02)),
                         ("beta1 <= 0.95", np.arange(0.30, 0.951, 0.02))):
        j, s = scan(deltas, alphas, betas)
        kj, ks = min(j), min(s)
        print(f"{label:14s}: joint Omega(n1 n2) <= {kj} (delta, alpha, lam, beta1 = {j[kj]});"
              f"  max side Omega(n_i) <= {ks} ({s[ks]})")
