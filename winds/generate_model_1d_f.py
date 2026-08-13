"""1-D physical models from the f-parameterized dimensionless runs (runs_f/).

Uses the new tiring parameterization  m_sph = 0.2/Γ₀  (Ṁ_out = Ṁ, h = 0.3),
which fixes the heating radius at  R = Ṁ G M_BH / (0.2 Γ₀ L_o)
= Ṁ G M_BH / (0.2 L_Edd), and a set of tiring fractions f = m/m_max.
f ∈ {0.60, 0.90} are the rounded stand-ins for the old solver's Γ₀ grid
(Γ₀ = 1.60, 1.40 had f = 0.592, 0.889); f = 0.94 stands in for Γ₀ = 1.375
(f = 0.948, past the solver's reach — see solve_f_grid.py).

For each (f, Ṁ) write:
  * runs_f/models/f{F}_Mdot{X}.npz — physical profile arrays
    (r, v, rho, tau) plus scalars, for plotting.
  * sirocco_imports_f/outflow_f{F}_Mdot{X}.import.txt with
        columns: i, r [cm], v_r [cm/s], ρ [g/cm³], T_e [K]
        rows:    1 inner-ghost cell + (n_cells + 1) shell edges

Radial domain of the Sirocco import:
    r_in  = radius where τ_es = 100 (full ionization), or R if τ_base < 100
    r_ph  = radius where τ_es = 1
    r_out = 10 × r_ph
Inner-ghost: r_in / 1.1; T_e = 10000 K in the wind, T_core at the ghost.
"""
import numpy as np
from pathlib import Path
from scipy.integrate import cumulative_trapezoid

# --- Constants (CGS) and fixed parameters ---
G_cgs    = 6.674e-8
Msun     = 1.989e33
yr       = 3.156e7
c_cgs    = 2.998e10
sigma_SB = 5.670e-5

M_BH   = 1e6 * Msun
kappa  = 0.34
L_Edd  = 4 * np.pi * G_cgs * M_BH * c_cgs / kappa

HERE        = Path(__file__).parent
RUNS        = HERE / 'runs_f'
MODELS_OUT  = RUNS / 'models'
SIROCCO_OUT = HERE / 'sirocco_imports_f'
MODELS_OUT.mkdir(exist_ok=True)
SIROCCO_OUT.mkdir(exist_ok=True)

F_SET     = [0.60, 0.90, 0.94]
N_CELLS   = 100
T_E_WIND  = 10000.0
Mdot_grid = [2.5, 5.0, 10.0, 15.0, 20.0]   # M☉/yr


def load_dimensionless(npz_path):
    """Augmented dimensionless grid (x, w) and optical-depth integral.

    The saved x-grid is uniform, so the outer region (1−x ≪ 1) is coarse;
    augment it with log-spaced points in 1−x (w is asymptotically constant
    there, so log-w vs log(1−x) interpolation is accurate). Returns the
    scalars dict, arrays (x, w), and τ(x)/τ★ (decreasing in x)."""
    d = np.load(npz_path)
    x_s, w_s = d['x'], np.maximum(d['w'], 1e-30)

    tail = 1 - np.geomspace(1e-9, 1 - x_s[-2], 400)[::-1]
    x = np.unique(np.concatenate([x_s, tail]))
    w = np.exp(np.interp(np.log(1 - x), np.log(1 - x_s)[::-1],
                         np.log(w_s)[::-1]))

    cum = cumulative_trapezoid(1 / np.sqrt(w), x, initial=0)
    tau_dimless = cum[-1] - cum
    scalars = {k: float(d[k]) for k in
               ('f', 'Gamma_0', 'm', 'w1', 'q0', 'w_sg')}
    scalars['m_max'] = float(d['m_max']) if 'm_max' in d.files else \
        1 - 1/scalars['Gamma_0']
    return scalars, x, w, tau_dimless


def x_at_tau(x, tau_dimless, t):
    return float(np.interp(-t, -tau_dimless, x))


def process(sc, x, w, tau_dimless, Mdot_Msun):
    f       = sc['f']
    Gamma_0 = sc['Gamma_0']
    Mdot    = Mdot_Msun * Msun / yr
    R       = Mdot * G_cgs * M_BH / (0.2 * L_Edd)
    v_esc   = np.sqrt(2 * G_cgs * M_BH / R)
    L_o     = Gamma_0 * L_Edd
    tau_star = kappa * Mdot / (4 * np.pi * R * v_esc)
    tau_base = tau_star * tau_dimless[0]

    # --- physical profile npz -------------------------------------------
    r   = R / (1 - x)
    v   = v_esc * np.sqrt(w)
    rho = Mdot / (4 * np.pi * r**2 * np.maximum(v, 1e-300))
    tau = tau_star * tau_dimless

    model_file = MODELS_OUT / f"f{f:.2f}_Mdot{Mdot_Msun:.1f}.npz"
    np.savez_compressed(
        model_file,
        M_BH=M_BH, kappa=kappa, Mdot=Mdot, Mdot_Msun_per_yr=Mdot_Msun,
        f=f, Gamma_0=Gamma_0, m=sc['m'], m_max=sc['m_max'],
        w1=sc['w1'], q0=sc['q0'], w_sg=sc['w_sg'],
        R=R, v_esc=v_esc, L_Edd=L_Edd, L_o=L_o,
        tau_star=tau_star, tau_base=tau_base,
        x=x, w=w, r=r, v=v, rho=rho, tau=tau,
    )

    # --- Sirocco import --------------------------------------------------
    log1mx = np.log(1 - x)
    logw   = np.log(np.maximum(w, 1e-30))

    def vr(rq):
        wq = np.exp(np.interp(np.log(R / rq), log1mx[::-1], logw[::-1]))
        return v_esc * np.sqrt(wq)

    def rho_of(rq):
        return Mdot / (4 * np.pi * rq**2 * vr(rq))

    if tau_base > 100.0:
        r_in = R / (1 - x_at_tau(x, tau_dimless, 100.0 / tau_star))
        inner_label = 'τ_es=100'
    else:
        r_in = R
        inner_label = f'R (τ_es={tau_base:.0f}<100 throughout)'
    r_ph  = R / (1 - x_at_tau(x, tau_dimless, 1.0 / tau_star))
    r_out = 10.0 * r_ph
    T_core = (L_o / (4 * np.pi * sigma_SB * r_in**2))**0.25

    r_edges    = np.logspace(np.log10(r_in), np.log10(r_out), N_CELLS + 1)
    r_ghost_in = r_in / 1.1
    r_all      = np.concatenate([[r_ghost_in], r_edges])
    r_centers  = np.sqrt(r_all[:-1] * r_all[1:])

    outfile = SIROCCO_OUT / f'outflow_f{f:.2f}_Mdot{Mdot_Msun:.1f}.import.txt'
    with outfile.open('w') as fh:
        fh.write(f"# 1D spherical import for Sirocco — Owocki+2017 super-Eddington wind\n")
        fh.write(f"# f=m/m_max={f}, Γ₀={Gamma_0:.6f}, m={sc['m']:.6f}, "
                 f"w_sg=0 (tiring_solver, runs_f/)\n")
        fh.write(f"# m_sph = 0.2/Γ₀ (Ṁ_out=Ṁ, h=0.3) → R = ṀGM/(0.2 L_Edd)\n")
        fh.write(f"# M_bh={M_BH:.3e} g, Ṁ={Mdot_Msun} M☉/yr\n")
        fh.write(f"# L_o={L_o:.3e} erg/s,  κ_es={kappa} cm²/g\n")
        fh.write(f"# Heating radius R = {R:.4e} cm; τ_es(R) ≈ {tau_base:.0f}\n")
        fh.write(f"# r_in  ({inner_label})  = {r_in:.4e} cm\n")
        fh.write(f"# r_ph  (τ_es=1)    = {r_ph:.4e} cm\n")
        fh.write(f"# r_out (=10×r_ph)  = {r_out:.4e} cm\n")
        fh.write(f"# n_cells={N_CELLS}, T_core={T_core:.0f} K\n")
        fh.write(f"# i  r(cm)  v_r(cm/s)  mass_rho(g/cm^3)  t_e(K)\n")
        for i in range(len(r_all)):
            r_c = r_centers[0] if i == 0 else r_centers[i - 1]
            t_e = T_core if i == 0 else T_E_WIND
            fh.write(f"{i}  {r_all[i]:.8e}  {vr(r_all[i]):.8e}  "
                     f"{rho_of(r_c):.8e}  {t_e:.1f}\n")

    return dict(model=model_file.name, sirocco=outfile.name, R=R,
                tau_base=tau_base, r_in=r_in, r_ph=r_ph, r_out=r_out,
                T_core=T_core, v_in=vr(r_in), v_ph=vr(r_ph))


def main():
    for F in F_SET:
        run = RUNS / f"f{F:.2f}.npz"
        sc, x, w, tau_dimless = load_dimensionless(run)
        print(f"{run.name}: f={sc['f']}, Γ₀={sc['Gamma_0']:.6f}, "
              f"m={sc['m']:.6f}, ∫dx/√w={tau_dimless[0]:.2f}")
        for Mdot in Mdot_grid:
            r = process(sc, x, w, tau_dimless, Mdot)
            print(f"  Ṁ={Mdot:5.1f}: τ_base={r['tau_base']:7.0f}  R={r['R']:.3e}  "
                  f"r_in={r['r_in']:.3e} ({r['r_in']/r['R']:.2f}R)  "
                  f"r_ph={r['r_ph']:.3e}  T_core={r['T_core']:.0f} K")


if __name__ == '__main__':
    main()
