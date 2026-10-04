# Even integers as sums of two sums of two squares with few prime factors

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23137626.svg)](https://doi.org/10.5281/zenodo.23137626)

**Author:** Ruqing Chen, GUT Geoservice Inc., Montreal, Canada (ruqing@hotmail.com)

Companion project to
[almost-prime-sum-of-two-squares](https://github.com/Ruqing1963/almost-prime-sum-of-two-squares)
(doi:10.5281/zenodo.23136112).

> **Theorem.** For every sufficiently large N ≡ 2 (mod 4), the number of n ≤ N, n ≡ 1 (mod 4), such that
> n and N − n are both sums of two squares and Ω(n(N − n)) ≤ 11 is at least (3/2)·𝔖(N)·N/(log N)²,
> where 𝔖(N) is the binary Goldbach singular series.

> **Corollary.** For the same N, at least (3/5)·𝔖(N)·N/(log N)² such n satisfy
> max(Ω(n), Ω(N − n)) ≤ 8 (and Ω(n(N − n)) ≤ 12).

> **Theorem 1.2.** The same holds for the principal form Q_Δ of each of the nine imaginary quadratic fields
> of class number one, with Ω(n(N − n)) ≤ k_Δ(N) and count ≫_Δ 𝔖(N)N/(log N)², N even:

| Δ | Q_Δ | k_Δ(N) = 11 | k_Δ(N) = 12, 13 |
|---|---|---|---|
| −3 | x² + xy + y² | N ≡ 2 (mod 3) | 12: N ≡ 1 (mod 3); 13: N ≡ 6 (mod 9) |
| −4 | x² + y² | N ≡ 2 (mod 4) | 13: N ≡ 4 (mod 8) |
| −7 | x² + xy + 2y² | 7 ∤ N | 13: 7 ∥ N |
| −8 | x² + 2y² | v₂(N) ≤ 2 | 13: v₂(N) = 3 |
| −11, −19, −43, −67, −163 | x² + xy + ((1−Δ)/4)y² | r ∤ N | 13: r ∥ N |

The extra prime factors are forced powers of the ramified prime r (code/local_conditions.py).
The main term and the switching term carry the same factor for every modulus M, so the certified
inequality applies to all nine Δ (Section 5 of the paper).

For N ≢ 2 (mod 4) the prime 2 forces extra factors of 2 (for N = 2^j the only representation is
n₁ = n₂ = 2^{j−1}), so a general statement should count the prime factors of the odd parts only.

## Status

Draft preprint (paper/, 8 pages). The final numerical inequality is certified with outward-rounded
interval arithmetic. The sieve lemmas are written out in the draft but have **not yet been independently
refereed**. A web literature search (Hooley, Indlekofer, Blomer, Brüdern–Fouvry, Blomer–Grimmelt–Li–Rydin
Myerson) found no overlapping result; a database (MathSciNet/zbMATH) search has not been done.

## Certified bounds (`code/certify.py`)

| δ | α | λ | β₁ | k = Ω(n₁n₂) bound | L ≥ | R ≤ | c ≤ | U ≤ | **net ≥** | constant |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.245 | 0.016 | 1/8 | 1/2 | **11** | 0.944016 | 2.830154 | 1.067736 | 1.134620 | **+0.038949** | 1.504 |
| 0.22 | 0.010 | 1/8 | 2/5 | 12 | 0.913373 | 3.144708 | 1.092388 | 1.166735 | **+0.039811** | 2.241 |

| 0.24 | 0.022 | 89/500 | 3/10 | 12 (max side **8**) | 0.919835 | 1.939896 | 0.878693 | 1.166966 | **+0.022031** | 0.607 |

Here k = ⌈1/λ + 2/β₁⌉ − 1 and the side bound is ⌈1/λ + 1/β₁⌉ − 1, computed in exact rational arithmetic; "constant" is the lower bound for
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
paper/Chen2026b_BilateralTwoSquares.tex / .pdf   draft paper (8 pages)
paper/figures/figure1_linear_sieve.pdf, figure2_switching_domain.pdf   Figures 1-2 (vector PDF)
code/plot_figures.py   regenerates Figures 1-2 (illustrative, floating point)
code/certify.py        interval-arithmetic certification of the three parameter sets (~20 s)
code/local_conditions.py  local conditions for the nine discriminants (Table 3 of the paper)
code/scan_beta1.py     joint and per-side Omega bounds, Richert range beta1 up to 0.95 (no gain)
code/estimate_modular_tradeoff.py  exploratory: weights r(n)r(N-n) at an assumed level theta (not part of the paper)
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
