"""
Local conditions for N = Q(x,y) + Q(z,w), Q the principal form of a fundamental discriminant D < 0 with
h(D) = 1. An integer n is represented by Q iff every inert prime has even exponent.

We need residue classes of n such that n = e*u, N - n = e'*u', where e, e' are powers of the ramified
prime r of D (r = 2 for D = -4, -8, r = -D otherwise), (u u', 2D) = 1, and chi_D(u) = chi_D(u') = 1
(=> an even number of inert prime factors on each side). For odd D, n and N - n must also be odd.
cost = Omega(e) + Omega(e') is the number of forced extra prime factors, so k_D(N) = 11 + cost.

chi_D(u) depends only on u modulo |D| (u odd), so the condition is local at the single prime r (plus
"n odd" at 2 for odd D, always satisfiable since N is even). We enumerate n modulo r^E for each class of N
modulo r^(E-2) and record the minimal cost. For large r the answer is predicted by the count
#{t mod r : chi(t) = chi(N' - t) = 1} = (r - 1 - 2 chi(N'))/4  (N' = N or N/r),
which we also check directly modulo r.
"""
DISCS = [-3, -4, -7, -8, -11, -19, -43, -67, -163]


def chi(D, u):
    if D == -4:
        return 1 if u % 4 == 1 else -1
    if D == -8:
        return 1 if u % 8 in (1, 3) else -1
    q = -D
    return 1 if pow(u % q, (q - 1) // 2, q) == 1 else -1


def side(D, r, x, mod):
    """(cost, u) for x modulo `mod`, or None if the valuation / chi(u) is not determined."""
    c = 0
    while x % r == 0:
        x //= r; mod //= r; c += 1
        if x == 0 or mod == 1:
            return None
    need = 8 if D in (-4, -8) else -D
    if mod % need:
        return None
    return c, x


def cost_table(D, E):
    r = 2 if D in (-4, -8) else -D
    mod = r ** E
    period = r ** (E - 2)
    res = {}
    for N in range(period):
        if D in (-4, -8) and N % 2:
            continue                                   # N even
        best = None
        for n in range(mod):
            a = side(D, r, n, mod); b = side(D, r, (N - n) % mod, mod)
            if a is None or b is None:
                continue
            if chi(D, a[1]) == 1 and chi(D, b[1]) == 1:
                best = a[0] + b[0] if best is None else min(best, a[0] + b[0])
        res[N] = best
    return r, period, res


def check_large(q):
    """Verify the count formula modulo q for q nmid N' (both chi = +1)."""
    ok = True
    for Np in range(1, q):
        cnt = sum(1 for t in range(1, q) if (Np - t) % q and chi(-q, t) == 1 and chi(-q, Np - t) == 1)
        ok &= cnt == (q - 1 - 2 * chi(-q, Np)) // 4 and cnt >= 1
    return ok


if __name__ == "__main__":
    for D in DISCS:
        if D in (-4, -8):
            r, period, res = cost_table(D, 8)
            even = {N: c for N, c in res.items()}
            groups = {}
            for N, c in even.items():
                groups.setdefault(c, []).append(N)
            desc = "; ".join(f"{'none' if c is None else 'k=' + str(11 + c)}: N = {sorted(v)} mod {period}"
                             for c, v in sorted(groups.items(), key=lambda kv: (kv[0] is None, kv[0] or 0)))
            print(f"D={D:5d} (r=2): {desc}")
        elif -D <= 11:
            r, period, res = cost_table(D, 4)
            groups = {}
            for N, c in res.items():
                groups.setdefault(c, []).append(N)
            desc = "; ".join(f"{'none' if c is None else 'k=' + str(11 + c)}: {len(v)} classes of N mod {period}"
                             + (f" {sorted(v)[:8]}" if c and len(v) < 12 else "")
                             for c, v in sorted(groups.items(), key=lambda kv: (kv[0] is None, kv[0] or 0)))
            print(f"D={D:5d} (r={r}): {desc}")
        else:
            q = -D
            print(f"D={D:5d} (r={q}): count formula verified mod {q}: {check_large(q)}  "
                  f"=> k=11 if {q} does not divide N, k=13 if {q} || N (both sides divisible by {q})")
