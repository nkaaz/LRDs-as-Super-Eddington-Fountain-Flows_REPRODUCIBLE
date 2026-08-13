"""1-D marginally-unbound ("dilute_redd") wind models -> Sirocco imports.

The mu_* family: spherical wind coasting at exactly the local escape velocity,
    v(r) = v_esc(r) = sqrt(2 G M_BH / r),   rho(r) = Mdot / (4 pi r^2 v(r)),
from a central mass M_BH = 1e6 Msun radiating L = L_Edd (kappa_es = 0.34).

Conventions (these reproduce the original production imports; running this
script for Mdot=10 regenerates runs/mu_m10/outflow.import.txt to <1e-4):
  * r_in  = 1.6e15 cm * (Mdot / 5 Msun/yr)      (the Eddington radius scaling)
  * r_out = 4.5771e19 cm * (Mdot / 10)**2        (~10x the tau_es=1 photosphere;
                                                  tau_es ~ Mdot/sqrt(r) => r_ph ~ Mdot^2)
  * tau_es(r_in) = 37.82 * sqrt(Mdot / 10)       (family-exact scaling)
  * source: diluted blackbody with W = 1/tau_es(r_in) and color temperature
    T_color = T_eff * tau_es(r_in)**0.25, where sigma T_eff^4 = L_Edd/(4 pi r_in^2)
    (grey-interior scaling T_rad ~ T_eff tau^{1/4}; L is conserved via W).
  * grid: 100 log cells (101 edges) from r_in to r_out + 1 inner ghost edge at
    r_in/1.1; densities are cell-centered (geometric mean of edges), with the
    ghost row repeating the first center; t_e = T_color at the ghost, 1e4 K in
    the wind. Import columns: i, r [cm], v_r [cm/s], rho [g/cm^3], t_e [K].

Usage:  python generate_model_1d_mu.py [Mdot1 Mdot2 ...]   (default: full grid)
Writes sirocco_imports_f/outflow_mu_Mdot{X}.import.txt next to this script.
"""
import sys
import numpy as np
from pathlib import Path

G, Msun, yr, sigma_SB = 6.674e-8, 1.989e33, 3.156e7, 5.670e-5
M_BH  = 1e6 * Msun
GM    = G * M_BH
L_Edd = 1.439e44
kappa = 0.34
TAU_REF, MDOT_REF = 37.82, 10.0      # tau_es(r_in) of the Mdot=10 production run

HERE = Path(__file__).parent
OUT  = HERE / 'sirocco_imports_f'
OUT.mkdir(exist_ok=True)

MDOT_GRID = [2.5, 5.0, 10.0, 15.0, 20.0]


def build(Mdot_msun):
    r_in  = 1.6e15 * (Mdot_msun / 5.0)
    r_out = 4.5771e19 * (Mdot_msun / 10.0) ** 2
    tau   = TAU_REF * np.sqrt(Mdot_msun / MDOT_REF)
    T_eff = (L_Edd / (4 * np.pi * r_in**2 * sigma_SB)) ** 0.25
    T_col = T_eff * tau ** 0.25
    W     = 1.0 / tau
    Mdot  = Mdot_msun * Msun / yr
    v_esc = lambda r: np.sqrt(2 * GM / r)
    rho   = lambda r: Mdot / (4 * np.pi * r**2 * v_esc(r))

    edges = np.logspace(np.log10(r_in), np.log10(r_out), 101)
    r_all = np.concatenate([[r_in / 1.1], edges])
    cen   = np.sqrt(r_all[:-1] * r_all[1:])

    outfile = OUT / f'outflow_mu_Mdot{Mdot_msun:g}.import.txt'
    with outfile.open('w') as fh:
        fh.write("# 1D spherical import for Sirocco — pure escape velocity\n")
        fh.write(f"# M_bh=1e6 Msun, L=L_Edd={L_Edd:.3e}, Mdot={Mdot_msun:g} Msun/yr\n")
        fh.write(f"# r_in={r_in:.4e} (= 1.6e15 * Mdot/5), tau_es(r_in)={tau:.2f}\n")
        fh.write(f"# T_eff={T_eff:.1f}, T_color={T_col:.1f}, W={W:.6f}\n")
        fh.write(f"# r_out={r_out:.4e}, n_cells=100\n")
        fh.write("# i   r(cm)   v_r(cm/s)   mass_rho(g/cm^3)   t_e(K)\n")
        for i in range(len(r_all)):
            rc = cen[0] if i == 0 else cen[i - 1]
            te = T_col if i == 0 else 10000.0
            fh.write(f"{i}  {r_all[i]:.8e}  {v_esc(r_all[i]):.8e}  "
                     f"{rho(rc):.8e}  {te:.1f}\n")
    print(f"Mdot={Mdot_msun:5g}: r_in={r_in:.3e}  tau={tau:6.2f}  "
          f"T_eff={T_eff:6.0f}  T_color={T_col:6.0f}  W={W:.6f}  -> {outfile.name}")
    print(f"           pf: Central_object.radius={r_in:.4e}  temp={T_col:.0f}  "
          f"dilution_factor={W:.6f}  Wind.radmax={r_out:.4e}")


if __name__ == '__main__':
    grid = [float(a) for a in sys.argv[1:]] or MDOT_GRID
    for md in grid:
        build(md)
