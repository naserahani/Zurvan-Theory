#!/usr/bin/env python3
"""
Generates Table 1 and Figure 1 of the paper: the predicted w(z) history of the
c_H = 1 model.

    dOmega_L/dln a = Omega_L (1 - Omega_L)(1 + 2 sqrt(Omega_L))
    w(z)           = -1/3 - (2/3) sqrt(Omega_L(z))

Starting value Omega_L0 = 0.7282 is the combined BAO + Pantheon+ best fit
obtained in p6_2_pantheon.py, not an external input.

Outputs, written next to this file:
    wz_table.tex   - the numeric rows of Table 1
    wz_figure.pdf  - Figure 1
"""

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OL0 = 0.7282
HERE = Path(__file__).parent

rhs = lambda t, y: [y[0] * (1 - y[0]) * (1 + 2 * np.sqrt(np.clip(y[0], 1e-30, 1)))]
w = lambda OL: -1.0 / 3.0 - (2.0 / 3.0) * np.sqrt(OL)


def omega_L(z):
    if np.isinf(z):
        return 0.0
    if z == 0.0:
        return OL0
    s = solve_ivp(rhs, (0.0, np.log(1.0 / (1.0 + z))), [OL0],
                  rtol=1e-12, atol=1e-15)
    return float(s.y[0, -1])


# ---- Table 1 -------------------------------------------------------------
rows = []
print(f"{'z':>8}{'Omega_L':>11}{'w':>10}")
for z in [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 10.0, np.inf]:
    OL = omega_L(z)
    lab = r"$\infty$" if np.isinf(z) else f"${z:.2f}$"
    val = r"$-1/3$" if np.isinf(z) else f"${w(OL):.3f}$"
    rows.append(f"{lab} & ${OL:.4f}$ & {val} \\\\")
    print(f"{'inf' if np.isinf(z) else format(z,'.2f'):>8}"
          f"{OL:>11.4f}{w(OL):>+10.4f}")
(HERE / "wz_table.tex").write_text("\n".join(rows) + "\n")

# ---- Figure 1 ------------------------------------------------------------
zg = np.linspace(0.0, 3.0, 400)
s = solve_ivp(rhs, (0.0, np.log(1.0 / 4.0)), [OL0], rtol=1e-12, atol=1e-15,
              dense_output=True)
wg = np.array([w(float(s.sol(np.log(1.0 / (1.0 + z)))[0])) for z in zg])

fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(zg, wg, lw=2.2, color="#b5533c",
        label=r"$c_H=1$:  $w=-\frac{1}{3}-\frac{2}{3}\sqrt{\Omega_\Lambda(z)}$")
ax.axhline(-1, color="0.35", lw=1.4, ls="--", label=r"$\Lambda$CDM")
ax.plot(zg, -0.42 - 1.75 * zg / (1 + zg), lw=1.5, color="#4a7c59", ls=":",
        label=r"DESI CPL fit ($w_0=-0.42,\ w_a=-1.75$)")
ax.scatter([0], [w(OL0)], s=34, color="#b5533c", zorder=5)
ax.annotate(rf"$w_0={w(OL0):.3f}$", (0.06, w(OL0)), fontsize=9, va="center")
ax.set_xlabel(r"$z$"); ax.set_ylabel(r"$w(z)$")
ax.set_xlim(0, 3); ax.set_ylim(-1.85, -0.35)
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
ax.grid(alpha=0.15, lw=0.6)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig(HERE / "wz_figure.pdf")
print("\nwrote wz_table.tex and wz_figure.pdf")
