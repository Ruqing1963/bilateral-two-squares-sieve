# Bilateral sieve for sums of two squares with few prime factors (work in progress)

**Author:** Ruqing Chen, GUT Geoservice Inc., Montreal, Canada (ruqing@hotmail.com)

Companion project to
[almost-prime-sum-of-two-squares](https://github.com/Ruqing1963/almost-prime-sum-of-two-squares)
(doi:10.5281/zenodo.23136112).

## Problem

For a large integer N ≡ 2 (mod 4), count representations

N = n₁ + n₂, n₁ = x² + y², n₂ = z² + w², Ω(n₁ n₂) ≤ k.

For N ≢ 2 (mod 4) the prime 2 forces extra factors of 2 (for N = 2^j the only representation is
n₁ = n₂ = 2^{j−1}), so the general statement should count the odd parts only.

## Status: draft with a certified main-term inequality

With delta = 0.22, alpha = 0.01, lambda = 1/8, beta1 = 2/5, `code/certify.py` proves (interval arithmetic)
net >= +0.039808, giving Omega(n1 n2) <= 12 and at least 2.24 S(N) N/(log N)^2 representations
(main-term level; the sieve lemmas are in the draft but not refereed).

## Exploration

`code/scan.py` evaluates a **heuristic main-term model** in floating point.
The level of distribution of {n(N − n)} is 1, but it is shared by all sieve components. Sifting the
inert primes (≡ 3 mod 4) up to √N with the linear sieve on n(N − n) lands exactly on s = 2, where
f₁(2) = 0. So the inert primes are sifted only up to N^{1/2−δ}, and survivors with two large inert factors
are removed by a switching upper bound at the Bombieri–Vinogradov level.

| scheme | description | smallest k for Ω(n₁n₂) |
|---|---|---|
| S1 | 2 linear sieves on n(N−n): inert up to N^{1/2−δ}, split up to N^α | 27 (no weights), **11** (Richert, δ = 0.245), 12 (δ = 0.22) |
| S2 | asymmetric: n clean up to √N, N−n truncated (3-vector) | no positive configuration |
| S3 | 4 semi-linear sieves (each side, inert and split) | 71 (no weights), 23 (Richert) |

Caveats of the model: the singular-series corrections are assumed to cancel between the main term and
the switching term; the vector sieve is the crude Brüdern–Fouvry inequality; nothing is interval-certified;
the literature check (Hooley, Indlekofer, Brüdern–Fouvry, …) is still to be done.

## Layout

```
code/sieve_core.py   beta-sieve functions (kappa = 1/2 and 1) and the switching density c(delta, alpha)
code/scan.py         schemes S1-S3 and parameter scan
results/             scan outputs
paper/               draft skeleton
```

Reproduce: `cd code && python scan.py` (about 40 s).
