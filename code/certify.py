"""
Rigorous (outward-rounded) check of the MAIN-TERM inequality of scheme S1 for
    N = n1 + n2,  n1, n2 sums of two squares,  N = 2 mod 4,
with the fixed parameters
    delta = 0.22, alpha = 0.010, lam = 1/8, beta1 = 2/5   (1/lam + 2/beta1 = 13  =>  Omega(n1 n2) <= 12).

    net = L - lam*R - bad,  in units of X V_I(y) V_S(z)  (see README / paper, Sections 2-4):
      L   = f1(s_I) F1(s_S) + F1(s_I) f1(s_S) - F1(s_I) F1(s_S),  s_I = theta_I/(1/2-delta), s_S = (1-theta_I)/alpha
      R   = int_alpha^beta1 (1 - u/beta1) F1(t/(1/2-delta)) F1((1-u-t)/alpha) du/u     (any split t = t(u))
      bad = 2 * 2 c(delta, alpha) e^gamma sqrt(alpha (1/2-delta)) * Fh(t'/(1/2-delta)) Fh((1/2-t')/alpha)
    F1, f1: linear sieve (kappa = 1, beta = 2);  Fh: semi-linear upper function (kappa = 1/2, beta = 1).
This certifies the numerical inequality only, not the sieve lemmas or the o(1) terms.

Assumptions (standard): F decreasing, f increasing; closed forms
    kappa=1:   F = 2e^g/s on (0,3],  f = 2e^g log(s-1)/s on [2,4];   (sF)' = f(s-1),  (sf)' = F(s-1)
    kappa=1/2: F = A/sqrt s on (0,2], f = A arccosh(sqrt s)/sqrt s on [1,3], A = 2 sqrt(e^g/pi);
               (sqrt s F)' = f(s-1)/(2 sqrt s), (sqrt s f)' = F(s-1)/(2 sqrt s).
Rounding: closed forms, logs and constants in mpmath interval arithmetic (80 bits); recurrences in IEEE
doubles with one outward nextafter step after each operation.
"""
import math
import numpy as np
from mpmath import iv
iv.prec = 80

dn = lambda x: math.nextafter(x, -math.inf)
up = lambda x: math.nextafter(x, math.inf)
EG_iv = iv.exp(iv.euler)

H_DEN = 1000                       # grid s = i/1000


def table(kappa, s_top):
    """Rigorous lower/upper enclosures of F and f on the grid i/H_DEN, 0 <= i <= s_top*H_DEN."""
    n = s_top * H_DEN
    beta = 1 if kappa == 0.5 else 2
    if kappa == 0.5:
        A = 2 * iv.sqrt(EG_iv / iv.pi)
        Fc = lambda s: A / iv.sqrt(s)
        fc = lambda s: A * iv.log(iv.sqrt(s) + iv.sqrt(s - 1)) / iv.sqrt(s)
        pw = iv.sqrt
    else:
        A = 2 * EG_iv
        Fc = lambda s: A / s
        fc = lambda s: A * iv.log(s - 1) / s
        pw = lambda s: s
    Flo = np.zeros(n + 1); Fhi = np.full(n + 1, math.inf); flo = np.zeros(n + 1); fhi = np.zeros(n + 1)
    pl = [0.0] * (n + 1); ph = [0.0] * (n + 1)          # enclosures of s^kappa on the grid
    for i in range(n + 1):
        v = pw(iv.mpf([i, i]) / H_DEN); pl[i], ph[i] = dn(float(v.a)), up(float(v.b))
    iF, if0, if1 = (beta + 1) * H_DEN, beta * H_DEN, (beta + 2) * H_DEN
    for i in range(1, if1 + 1):
        s = iv.mpf([i, i]) / H_DEN
        if i <= iF:
            v = Fc(s); Flo[i], Fhi[i] = dn(float(v.a)), up(float(v.b))
        if if0 < i <= if1:
            v = fc(s); flo[i], fhi[i] = max(0.0, dn(float(v.a))), up(float(v.b))
    gF_lo, gF_hi = dn(pl[iF] * Flo[iF]), up(ph[iF] * Fhi[iF])
    gf_lo, gf_hi = dn(pl[if1] * flo[if1]), up(ph[if1] * fhi[if1])
    k = float(kappa)
    for i in range(iF, n):
        j = i - H_DEN
        w_lo, w_hi = dn(pl[i + 1] - ph[i]), up(ph[i + 1] - pl[i])      # int of kappa t^(kappa-1) dt
        gF_lo = dn(gF_lo + dn(flo[j] * w_lo)); gF_hi = up(gF_hi + up(fhi[j + 1] * w_hi))
        Flo[i + 1] = dn(gF_lo / ph[i + 1]); Fhi[i + 1] = up(gF_hi / pl[i + 1])
        if i >= if1:
            gf_lo = dn(gf_lo + dn(Flo[j + 1] * w_lo)); gf_hi = up(gf_hi + up(Fhi[j] * w_hi))
            flo[i + 1] = dn(gf_lo / ph[i + 1]); fhi[i + 1] = up(gf_hi / pl[i + 1])
    return dict(Flo=Flo, Fhi=Fhi, flo=flo, fhi=fhi, n=n, A_hi=up(float(A.b)), kappa=kappa, beta=beta)


def F_up(T, s):
    """Upper bound for F(s), s > 0 (F decreasing: use the grid point at or below s)."""
    if s <= T["beta"] + 1 - 1e-9:                     # closed form A s^-kappa, rounded up
        return up(T["A_hi"] / dn(s ** T["kappa"]) * (1 + 1e-15))
    i = int(math.floor(s * H_DEN)) - 1
    return T["Fhi"][min(i, T["n"])]


def F_dn(T, s):
    i = int(math.ceil(s * H_DEN)) + 1
    return T["Flo"][i] if i <= T["n"] else 0.0


def f_dn(T, s):
    i = int(math.floor(s * H_DEN)) - 1
    return T["flo"][min(i, T["n"])] if i >= 0 else 0.0


lin = table(1.0, 110)
half = table(0.5, 60)
print(f"linear:      F(3) in [{lin['Flo'][3000]:.6f},{lin['Fhi'][3000]:.6f}], f(4) in [{lin['flo'][4000]:.6f},"
      f"{lin['fhi'][4000]:.6f}], F(110) <= {lin['Fhi'][-1]:.6f}, f(110) >= {lin['flo'][-1]:.6f}")
print(f"semi-linear: F(2.5) in [{half['Flo'][2500]:.6f},{half['Fhi'][2500]:.6f}], f(60) >= {half['flo'][-1]:.6f}")



def certify(delta_s, alpha_s, lam, beta1, label):
    """delta_s, alpha_s: exact decimal strings; lam, beta1: exactly representable binary fractions or decimals
    (only used through 1/lam + 2/beta1, checked in exact rational arithmetic, and as float factors)."""
    from fractions import Fraction
    delta, alpha = float(delta_s), float(alpha_s)
    ey = 0.5 - delta                                       # log y / log N
    ey_s = str(Fraction(1, 2) - Fraction(delta_s))
    ey_iv = iv.mpf(Fraction(ey_s).numerator) / Fraction(ey_s).denominator
    al_iv = iv.mpf(Fraction(alpha_s).numerator) / Fraction(alpha_s).denominator
    bound = 1 / Fraction(str(lam)) + 2 / Fraction(str(beta1))
    k = int(-(-bound // 1)) - 1                            # Omega < bound  =>  Omega <= ceil(bound) - 1
    print(f"\n=== {label}: delta={delta_s}, alpha={alpha_s}, lam={lam}, beta1={beta1};"
          f" 1/lam + 2/beta1 = {bound} => Omega(n1 n2) <= {k} ===")

    # ---- L (lower) ----
    best = None
    for tI in np.arange(0.60, 0.995, 0.0005):
        sI, sS = tI / ey, (1 - tI) / alpha
        v = f_dn(lin, sI) * F_dn(lin, sS) + F_up(lin, sI) * (f_dn(lin, sS) - F_up(lin, sS))
        if best is None or v > best[0]:
            best = (v, tI)
    tI = best[1]; sI, sS = tI / ey, (1 - tI) / alpha
    L_lo = dn(dn(f_dn(lin, sI) * F_dn(lin, sS)) + dn(F_up(lin, sI) * dn(f_dn(lin, sS) - F_up(lin, sS))))

    # ---- R (upper) ----
    NU = 4000
    R_hi = 0.0
    for kk in range(NU):
        ua, ub = dn(alpha + (beta1 - alpha) * kk / NU), up(alpha + (beta1 - alpha) * (kk + 1) / NU)
        lev = 1 - ub
        bestU = min(F_up(lin, t / ey) * F_up(lin, dn((1 - ub - t) / alpha))
                    for t in np.linspace(lev * 0.5, lev * 0.995, 80))
        weight = up(up(1 - ua / beta1) / ua)
        R_hi = up(R_hi + up(up(weight * bestU) * up(ub - ua)))

    # ---- c (upper) ----
    # I(mu) = log((a-mu)/b)/(1-mu), a = 1/2+delta, b = 1/2-delta. Its derivative has the sign of
    # log((a-mu)/b) - (1-mu)/(a-mu); the first term decreases and the second increases in mu (a < 1),
    # so I is decreasing on [0, 2 delta] as soon as log(a/b) < 1/a.
    a_iv = 1 - ey_iv
    assert float(iv.log(a_iv / ey_iv).b) < float((1 / a_iv).a), "monotonicity of I(mu) not established"
    i_al = int(round(alpha * 2000)); i_top = int(round(2 * delta * 2000))
    assert abs(i_al - alpha * 2000) < 1e-9 and abs(i_top - 2 * delta * 2000) < 1e-9
    mass = np.zeros(i_top + 1)
    for kq in range(i_al, i_top):
        mass[kq] = up(float((iv.mpf('0.5') * iv.log(iv.mpf([kq + 1, kq + 1]) / kq)).b))
    total = np.zeros(i_top + 1); total[0] = 1.0; term = total.copy()
    for j in range(1, i_top // i_al + 1):
        term = np.convolve(term, mass)[: i_top + 1] / j
        total += term
    total *= (1 + 1e-9)
    c_hi = 0.0
    for kq in range(i_top + 1):
        if total[kq] == 0:
            continue
        mu = iv.mpf([kq, kq]) / 2000
        hi_ = (1 - mu) / 2
        if hi_.a <= ey_iv.b:
            continue
        aa = 1 - mu
        I = (iv.log(hi_ / ey_iv) - iv.log((aa - hi_) / (aa - ey_iv))) / aa
        c_hi = up(c_hi + up(total[kq] * up(0.25 * float(I.b))))

    # ---- bad (upper) ----
    Us = min(F_up(half, t / ey) * F_up(half, dn((0.5 - t) / alpha)) for t in np.linspace(0.15, 0.495, 700))
    pref = 4 * EG_iv * iv.sqrt(al_iv * ey_iv)
    bad_hi = up(up(float(pref.b) * c_hi) * Us)

    net_lo = dn(dn(L_lo - up(lam * R_hi)) - bad_hi)
    scale = iv.exp(-2 * iv.euler) / (2 * al_iv * ey_iv)
    print(f"theta_I = {tI:.4f} (s_I = {sI:.4f}, s_S = {sS:.2f})")
    print(f"L   >= {L_lo:.6f}")
    print(f"R   <= {R_hi:.6f}")
    print(f"c   <= {c_hi:.6f}")
    print(f"U_side <= {Us:.6f}")
    print(f"bad <= {bad_hi:.6f}")
    print(f"certified:  net >= {net_lo:+.6f}   ->  {'POSITIVE' if net_lo > 0 else 'NOT certified'}")
    if net_lo > 0:
        print(f"lower bound constant: #reps >= {dn(net_lo * float(scale.a)):.5f} * S(N) N/(log N)^2 * (1+o(1))")
    return net_lo


if __name__ == "__main__":
    certify("0.22", "0.010", 0.125, 0.40, "k = 12")
    certify("0.245", "0.016", 0.125, 0.50, "k = 11")
