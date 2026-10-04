"""
Figures for the paper (illustrative, floating point; not part of the certified computation).
  paper/figures/figure1_linear_sieve.pdf      F_1, f_1 with the threshold s = 2 and the working points s_I, s_S
  paper/figures/figure2_switching_domain.pdf  the switching region in the (u1, u2) plane, q_i = N^{u_i}
Run from the repository root:  python code/plot_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from sieve_core import F_lin, f_lin

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams["pdf.fonttype"] = 42

# certified parameter set for k = 11 (see code/certify.py)
DELTA, ALPHA, THETA_I = 0.245, 0.016, 0.9220
S_I = THETA_I / (0.5 - DELTA)

# ---- Figure 1 ---------------------------------------------------------------------------------------------
s = np.linspace(0.5, 4.5, 800)
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.plot(s, F_lin(s), label=r"$F_1(s)$")
ax.plot(s, f_lin(s), label=r"$f_1(s)$")
ax.axhline(1, color="0.6", lw=0.8)
ax.axvline(2, color="0.4", lw=0.8, ls=":")
ax.annotate(r"$f_1(2)=0$", xy=(2, 0), xytext=(2.25, 0.35), arrowprops=dict(arrowstyle="->", lw=0.7))
ax.axvline(S_I, color="C2", lw=1.0, ls="--")
ax.annotate(rf"$s_I\approx{S_I:.2f}$", xy=(S_I, float(f_lin(S_I))), xytext=(S_I - 1.05, 1.9),
            arrowprops=dict(arrowstyle="->", lw=0.7, color="C2"), color="C2")
ax.set_xlabel("s"); ax.set_ylim(0, 3.2); ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "figure1_linear_sieve.pdf"), format="pdf", bbox_inches="tight")

# ---- Figure 2 ---------------------------------------------------------------------------------------------
lo = 0.5 - DELTA
fig, ax = plt.subplots(figsize=(5.2, 4.2))
# m = 1: lo <= u1 <= u2, u1 + u2 <= 1
tri = Polygon([(lo, lo), (0.5, 0.5), (lo, 1 - lo)], closed=True, fc="C0", alpha=0.25, ec="C0",
              label=r"$m=1$:  $u_1+u_2\leq 1$")
ax.add_patch(tri)
# m > 1 (m >= N^alpha): u1 + u2 <= 1 - alpha
tri2 = Polygon([(lo, lo), (0.5 - ALPHA / 2, 0.5 - ALPHA / 2), (lo, 1 - ALPHA - lo)], closed=True,
               fc="C1", alpha=0.35, ec="C1", label=r"$m>1$:  $u_1+u_2\leq 1-\mu$, $\mu\geq\alpha$")
ax.add_patch(tri2)
x = np.linspace(0.2, 0.8, 10)
ax.plot(x, x, color="0.5", lw=0.7, ls=":")
ax.plot(x, 1 - x, color="0.5", lw=0.7, ls=":")
ax.axvline(lo, color="0.4", lw=0.8, ls="--")
ax.text(lo + 0.004, 0.27, r"$u_1=\frac{1}{2}-\delta$", fontsize=9)
ax.text(0.39, 0.70, r"$q_1q_2m$ on one side:" "\n" r"$\mathcal{E}_1$ (in $n$) or $\mathcal{E}_2$ (in $N-n$)",
        fontsize=8.5)
ax.set_xlim(0.22, 0.55); ax.set_ylim(0.22, 0.8)
ax.set_xlabel(r"$u_1=\log q_1/\log N$"); ax.set_ylabel(r"$u_2=\log q_2/\log N$")
ax.legend(frameon=False, loc="lower right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "figure2_switching_domain.pdf"), format="pdf", bbox_inches="tight")
print("written:", os.path.join(OUT, "figure1_linear_sieve.pdf"), os.path.join(OUT, "figure2_switching_domain.pdf"))
