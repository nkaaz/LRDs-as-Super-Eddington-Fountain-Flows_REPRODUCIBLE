"""Solve the wind at fixed tiring fractions f = m/m_max using tiring_solver.

New-model parameterization (Ṁ_out/Ṁ = 1, h = 0.3):
    m = 0.2/Γ₀  and  m_max = 1 − 1/Γ₀
so a target f = m/m_max fixes both:
    Γ₀ = 1 + 0.2/f,   m = 0.2/Γ₀ = f·m_max.

Outputs: runs_f/f{F}.npz with the dimensionless profile on a fine mesh.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / 'tiring_solver'))
from wind import solve, m_max, ode_residual

OUT = HERE / 'runs_f'
OUT.mkdir(exist_ok=True)

# 0.95 is out of reach at these low Gamma_0 (continuation stalls at f~0.947,
# and a 250k-node mesh cap does not help); 0.94 is the closest point the
# solver can validate. See tiring_solver README.
# 0.60 and 0.90 are the rounded stand-ins for the old solver's Γ₀ grid
# (Γ₀ = 1.60, 1.40, 1.375 had f = 0.592, 0.889, 0.948 -> 0.60, 0.90, 0.94).
F_GRID = [0.60, 0.75, 0.85, 0.90, 0.94, 0.95]


def main():
    # Optional CLI: solve only the listed f values, e.g. `solve_f_grid.py 0.6 0.9`
    targets = [float(a) for a in sys.argv[1:]] or F_GRID
    for f in targets:
        Gamma_0 = 1 + 0.2 / f
        mm = m_max(Gamma_0)
        m = 0.2 / Gamma_0
        print(f"f={f}: Gamma_0={Gamma_0:.6f}, m={m:.6f}, m_max={mm:.6f} "
              f"(check f={m/mm:.6f})", flush=True)
        try:
            s = solve(m, Gamma_0, verbose=True)
        except (RuntimeError, ValueError) as e:
            print(f"  FAILED: {e}", flush=True)
            continue

        x = np.linspace(1e-9, 1 - 1e-9, 4000)
        w, q = s.sol(x)
        p = q * (1 - x)**2
        res = ode_residual(s)
        outfile = OUT / f"f{f:.2f}.npz"
        np.savez_compressed(
            outfile,
            f=f, Gamma_0=Gamma_0, m=s.m, m_max=mm,
            w1=s.w1, q0=s.q0, w_sg=s.w_sg, ode_residual=res,
            x=x, w=w, q=q, p=p,
        )
        print(f"  -> {outfile.name}: w1={s.w1:.4e}, q0={s.q0:.4e}, "
              f"v_inf/v_esc={np.sqrt(s.w1):.4f}, residual={res:.2e}", flush=True)


if __name__ == '__main__':
    main()
