from pathlib import Path
import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
from p6_3_joint import chi2_all
res={}
for model,lab,x0 in (("lcdm","flat LCDM",[0.300,0.685,0.02253]),
                     ("hde","HDE c_H=1",[0.320,0.642,0.02326])):
    r=minimize(lambda p: chi2_all(*p,model,("cmb","bao","sn")), x0,
               method="Nelder-Mead", options=dict(xatol=1e-7,fatol=1e-8,maxiter=6000))
    Om,h,ob=r.x
    parts={k:chi2_all(Om,h,ob,model,(k,)) for k in ("cmb","bao","sn")}
    res[lab]=(Om,h,ob,r.fun,parts)
    print(f"\n{lab}")
    print(f"   Omega_m0={Om:.4f}  H0={100*h:.2f}  omega_b={ob:.5f}  omega_m={Om*h**2:.5f}")
    print(f"   chi2 = {r.fun:.2f}   CMB {parts['cmb']:.2f}/3 | BAO {parts['bao']:.2f}/13 | SN {parts['sn']:.2f}/1580")
a=res["HDE c_H=1"]; b=res["flat LCDM"]
print(f"\nDelta chi2 total (HDE-LCDM) = {a[3]-b[3]:+.2f}  [3 params each]")
for k in ("cmb","bao","sn"):
    print(f"   from {k.upper():4s}: {a[4][k]-b[4][k]:+9.2f}")
print(f"\nlate-time-only HDE fit gave Omega_m0 = 0.2718; joint pushes it to {a[0]:.4f}")
np.save(Path(__file__).with_name("p63.npy"), res, allow_pickle=True)
