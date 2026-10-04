# Even integers as sums of two sums of two squares with few prime factors

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23137626.svg)](https://doi.org/10.5281/zenodo.23137626)

**Author:** Ruqing Chen, GUT Geoservice Inc., Montreal, Canada (ruqing@hotmail.com)

Companion project to
[almost-prime-sum-of-two-squares](https://github.com/Ruqing1963/almost-prime-sum-of-two-squares)
(doi:10.5281/zenodo.23136112).

> **Theorem.** For every sufficiently large N ≡ 2 (mod 4), the number of n ≤ N, n ≡ 1 (mod 4), such that
> n and N − n are both sums of two squares and Ω(n(N − n)) ≤ 11 is at least (3/2)·𝔖(N)·N/(log N)²,
> where 𝔖(N) is the binary Goldbach singular series.

For N ≢ 2 (mod 4) the prime 2 forces extra factors of 2 (for N = 2^j the only representation is
n₁ = n₂ = 2^{j−1}), so a general statement should count the prime factors of the odd parts only.

## Status

Draft preprint (paper/, 5 pages). The final numerical inequality is certified with outward-rounded
interval arithmetic. The sieve lemmas are written out in the draft but have **not yet been independently
refereed**. A web literature search (Hooley, Indlekofer, Blomer, Brüdern–Fouvry, Blomer–Grimmelt–Li–Rydin
Myerson) found no overlapping result; a database (MathSciNet/zbMATH) search has not been done.

## Certified bounds (`code/certify.py`)

| δ | α | λ | β₁ | k = Ω(n₁n₂) bound | L ≥ | R ≤ | c ≤ | U ≤ | **net ≥** | constant |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.245 | 0.016 | 1/8 | 1/2 | **11** | 0.944016 | 2.830154 | 1.067736 | 1.134620 | **+0.038949** | 1.504 |
| 0.22 | 0.010 | 1/8 | 2/5 | 12 | 0.913373 | 3.144708 | 1.092388 | 1.166735 | **+0.039811** | 2.241 |

Here k = ⌈1/λ + 2/β₁⌉ − 1, computed in exact rational arithmetic; "constant" is the lower bound for
(number of n)/(𝔖(N) N/(log N)²).

## Method in brief

The level of distribution of {n(N − n)} is 1 (trivial remainder |r_d| ≤ ρ(d)), shared by all sieve
components. The inert primes (≡ 3 mod 4) form a sieve of dimension 1 on n(N − n); sifting them up to √N at
level N lands on s = 2, where the linear lower bound function vanishes. They are therefore sifted only up to
N^{1/2−δ}. Since n ≡ N − n ≡ 1 (mod 4), a survivor that is not a pair of sums of two squares has two large
inert factors on the **same** side, and these are removed by a switching upper bound at the
Bombieri–Vinogradov level. The local factors of the main term and the switching term coincide, so the
singular series cancels exactly. Richert's weights act on n(N − n) jointly.

## Exploration (floating point, `code/scan.py`)

| scheme | description | smallest k for Ω(n₁n₂) |
|---|---|---|
| S1 | 2 linear sieves on n(N−n): inert up to N^{1/2−δ}, split up to N^α | 27 (no weights), **11** (Richert) |
| S2 | asymmetric: n clean up to √N, N−n truncated (3-vector) | no positive configuration |
| S3 | 4 semi-linear sieves (each side, inert and split) | 71 (no weights), 23 (Richert) |

## Layout

```
paper/Chen2026b_BilateralTwoSquares.tex / .pdf   draft paper
code/certify.py        interval-arithmetic certification of both parameter sets (~15 s)
code/sieve_core.py     floating-point sieve functions and switching density
code/scan.py           schemes S1-S3 and parameter scan (~40 s)
results/               outputs of certify.py and scan.py
```

Requirements: Python ≥ 3.9, numpy, mpmath. Reproduce: `cd code && python certify.py`.

## Citation

```bibtex
@misc{Chen2026BilateralTwoSquares,
  author = {Chen, Ruqing},
  title  = {Even integers as sums of two sums of two squares with few prime factors},
  year   = {2026},
  doi    = {10.5281/zenodo.23137626},
  url    = {https://doi.org/10.5281/zenodo.23137626}
}
```

## AI assistance

Parts of the computations and of the draft were produced with the help of the Claude language model.
The author is responsible for the content.

## License

Code: MIT (see `LICENSE`). Paper: CC BY 4.0.
