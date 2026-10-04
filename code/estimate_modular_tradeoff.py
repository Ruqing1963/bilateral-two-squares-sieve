"""
Trade-off model: "modular route" (weights r(n) r(N-n), level theta) versus the paper's route (sift inert
primes on n(N-n), level 1, switching).  Floating point, heuristic main-term model, NOT certified.

Modular route.  With the weight r(n) r(N-n), positivity already forces both sides to be sums of two squares,
so no inert sieve and no switching are needed.  Under this weight a split prime q divides n with density
~ 2/q and an inert prime with density ~ 1/q^2.  Each side is therefore a linear sieve (kappa = 1) in the split
primes, and we sift all primes < z = N^alpha on both sides:
    vector sieve of two linear components, levels theta_A + theta_B = theta,  s_X = theta_X / alpha,
    L = f1(s_A) F1(s_B) + F1(s_A) f1(s_B) - F1(s_A) F1(s_B)
    R = 2 int_alpha^beta1 (1 - u/beta1) U(theta - u) du/u      (split primes, density 2/q on each side)
    positivity: L - lam R > 0;  Omega(n1 n2) < 1/lam + 2/beta1, max Omega(n_i) < 1/lam + 1/beta1,
    and Omega(n1 n2) < 2/alpha.
The level theta of r(n) r(N-n) in arithmetic progressions is the unknown; we tabulate theta.
A genuine two-dimensional (Diamond-Halberstam-Richert) sieve would do somewhat better than the vector of two
linear sieves; the same crude vector sieve is used for the paper's route, so the comparison is like for like.
"""
import numpy as np
from scan import vec_lower, vec_upper


def modular_k(theta, alphas, betas):
    best_joint, best_side, best_now = None, None, None
    for a in alphas:
        comps = [("lin", a), ("lin", a)]
        L = vec_lower(comps, theta, 0.005)
        if L <= 0:
            continue
        k0 = int(np.ceil(2 / a - 1e-9)) - 1
        best_now = k0 if best_now is None else min(best_now, k0)
        u = np.linspace(a, min(0.95, theta - 0.01), 80)
        U = np.array([vec_upper(comps, theta - x, 0.01) for x in u])
        for b1 in betas:
            if b1 <= a or b1 > u[-1]:
                continue
            m = u <= b1
            if m.sum() < 4:
                continue
            R = 2 * np.trapz((1 - u[m] / b1) * U[m] / u[m], u[m])
            lam = 0.95 * L / R
            kj = int(np.ceil(1 / lam + 2 / b1 - 1e-9)) - 1
            ks = min(int(np.ceil(1 / lam + 1 / b1 - 1e-9)) - 1, kj - 1)
            best_joint = kj if best_joint is None else min(best_joint, kj)
            best_side = ks if best_side is None else min(best_side, ks)
    return best_now, best_joint, best_side


if __name__ == "__main__":
    alphas = np.arange(0.005, 0.30, 0.005)
    betas = np.arange(0.10, 0.951, 0.02)
    print("paper's route (level 1, inert sieve + switching):  no weights 27,  Richert: joint 11, max side 8")
    print("modular route, level theta of r(n)r(N-n):")
    print(" theta | no weights | Richert joint | Richert max side")
    for theta in (1/3, 0.40, 0.45, 0.50, 0.55, 0.60, 0.75, 1.00):
        k0, kj, ks = modular_k(theta, alphas, betas)
        print(f" {theta:5.3f} | {k0!s:>10} | {kj!s:>13} | {ks!s:>16}")
