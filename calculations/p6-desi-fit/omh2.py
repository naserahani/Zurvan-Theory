import numpy as np
from scipy.optimize import brentq
# DESI fitting formula for the drag-epoch sound horizon
def rd(om_h2, ob_h2=0.02236, Neff=3.04):
    return 147.05*(ob_h2/0.02236)**-0.13*(om_h2/0.1432)**-0.23*(Neff/3.04)**-0.1

print("Self-consistent physical matter density implied by BAO + each model")
print("  (BAO fixes X = h*r_d; the model fixes Omega_m0; r_d depends on omega_m)\n")
for lab, Om, X in (("flat LCDM", 0.2973, 101.54), ("HDE c_H=1", 0.2717, 99.75)):
    # solve  om_h2 = Om * (X/rd(om_h2))^2
    f = lambda w: Om*(X/rd(w))**2 - w
    w = brentq(f, 0.05, 0.30, xtol=1e-8)
    h = X/rd(w)
    print(f"  {lab:<11} omega_m = {w:.4f}   h = {h:.4f}  ->  H0 = {100*h:.2f}"
          f"   r_d = {rd(w):.2f} Mpc")
print(f"\n  Planck 2018  omega_m = 0.1430 +/- 0.0011   (early-universe, "
      f"largely independent of late-time DE)")
for lab, Om, X in (("flat LCDM", 0.2973, 101.54), ("HDE c_H=1", 0.2717, 99.75)):
    f = lambda w: Om*(X/rd(w))**2 - w
    w = brentq(f, 0.05, 0.30, xtol=1e-8)
    print(f"  {lab:<11} deviation = {(w-0.1430)/0.0011:+.1f} sigma  (naive)")
