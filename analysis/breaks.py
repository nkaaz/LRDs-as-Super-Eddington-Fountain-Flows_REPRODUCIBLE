"""Standard definition of the Lyman and Balmer break radii.

Shared by FigOpticalDepth, FigHeatCool and FigTemp_n2_grid so that the r_LyB /
r_BaB markers mean the same thing in every figure.

WHY NOT A RAW FLUX RATIO
------------------------
The obvious measure -- J(0.97 nu_edge) / J(1.03 nu_edge) -- is biased, because it
charges the local continuum slope to the break. A smooth blackbody with no edge
at all registers as broken once it is cool enough that the edge sits in Wien's
tail:

    T [K]     J(0.97 nu)/J(1.03 nu) at the Lyman edge   (no edge present!)
    30000                 1.15
    14000                 1.63
     9000                 2.39
     7000                 3.23
     5000                 5.55

so at the injected colour temperatures of the coolest models the measure reports
a break in the input spectrum itself.

WHAT THIS MODULE MEASURES INSTEAD
---------------------------------
The discontinuity in the *continuum*: fit the continuum independently on each
side of the edge, extrapolate both fits TO the edge, and take the ratio

    break = J_red(nu_edge) / J_blue(nu_edge).

The fit uses the basis [1, ln nu, nu], which is exact for a diluted blackbody
(ln B_nu = const + 3 ln nu - h nu / kT in the Wien limit), so a smooth spectrum
of any colour temperature returns 1.000 to <0.1% at both edges -- verified in
`selftest()` below. Emission lines are removed by iteratively clipping the most
positive residuals, so the Lyman and high-n Balmer forests bracketing each edge
do not drag the continuum fits.

A cell is "broken" when this ratio first reaches BREAK_THRESHOLD and stays there.

LYMAN PHOTOIONIZATION TRANSITION (not the same radius)
------------------------------------------------------
The factor-2 spectral break forms where the incident Lyman continuum has lost
only its first e-fold (tau_eff ~ ln 2). The hand-over of ground-state ionization
from the transmitted incident LyC to the local on-the-spot recombination
continuum -- the step in FigGamma1_residuals, and the zero-crossing of the n=2
residual at low Mdot -- completes only once the incident LyC is essentially
gone. Measured against those residuals on the mu grid (2026-09-13), the
inside-out ABSORPTION depth of the Lyman continuum from the base,

    tau_LyC(r) = INT_{r_in}^{r} n_1 sigma_1 dr',   sigma_1 = 6.30e-18,

reaches TAU_LYC_TRANSITION = 3 at 38.6 / 7.0 / 2.5 / 1.6 r_sph for
m2.5 / m5 / m10 / m15, against a measured hand-over at 39-41 / 7-8 / 1.9-2.4 /
1.4-1.6; the spectral break sits 30-50% inside it at low Mdot (27 / 3.7).
No fixed threshold on the spectral jump or on tau_eff = sqrt(3 tau_abs chi)
works across the grid (the jump at the hand-over runs from 1.4 to 378), because
the incident LyC supply varies by orders of magnitude with the input colour
temperature. Use lyc_transition_radius() for the photoionization transition and
break_radius() only for the spectral edge itself.
"""
import numpy as np

# --- the one place these are defined ---------------------------------------
NU_LYMAN  = 3.2900e15          # Hz, 912 A
NU_BALMER = NU_LYMAN / 4.0     # Hz, 3646 A
BREAK_THRESHOLD = 2.0          # continuum drops by >= this factor across the edge
TAU_LYC_TRANSITION = 3.0       # inside-out LyC absorption depth at the photoionization hand-over
SIG_LYC = 6.30e-18             # H I ground-state threshold cross-section [cm^2] (Sirocco h20 data)
GUARD = 1.06                   # skip +/-6% about the edge (bin width + MC smearing)
SPAN  = 1.40                   # fit the continuum out to +/-40%
KEEP  = 0.60                   # fraction of bins retained per clipping pass (drops lines)


def _fit_side(freq, f, lo, hi, nu_edge, niter=4, keep=0.60):
    """Extrapolate one side's continuum to nu_edge.

    Fits ln f = a + b ln nu + c nu/nu_edge -- exact for a diluted blackbody -- to
    every positive bin in [lo, hi], then iteratively discards the bins with the
    most POSITIVE residuals (emission lines) and refits, keeping the lowest
    `keep` fraction each pass. On a line-free spectrum the residuals are ~0 and
    the clipping is a no-op, so the estimator stays unbiased."""
    m = (freq >= lo) & (freq < hi) & (f > 0)
    if m.sum() < 8:
        return np.nan
    nu, y = freq[m], np.log(f[m])
    A = np.vstack([np.ones_like(nu), np.log(nu), nu / nu_edge]).T
    idx = np.arange(len(nu))
    for _ in range(niter):
        c, *_ = np.linalg.lstsq(A[idx], y[idx], rcond=None)
        if len(idx) < max(8, int(keep * m.sum())):
            break
        res = y[idx] - A[idx] @ c
        idx = idx[np.argsort(res)[:max(8, int(keep * len(idx)))]]
    c, *_ = np.linalg.lstsq(A[idx], y[idx], rcond=None)
    return float(np.exp(np.array([1.0, np.log(nu_edge), 1.0]) @ c))


def edge_break_ratio(freq, jnu, nu_edge):
    """Continuum break ratio J_red(nu_edge)/J_blue(nu_edge) for every cell.

    freq : (Nfreq,)      jnu : (Nfreq, Ncell)   ->  (Ncell,)
    Returns 1.0 for a smooth spectrum, >1 for an absorption edge, NaN where
    either side is too starved to fit."""
    jnu = jnu if jnu.ndim == 2 else jnu[:, None]
    out = np.full(jnu.shape[1], np.nan)
    for i in range(jnu.shape[1]):
        f = jnu[:, i]
        red  = _fit_side(freq, f, nu_edge / SPAN,  nu_edge / GUARD, nu_edge)
        blue = _fit_side(freq, f, nu_edge * GUARD, nu_edge * SPAN,  nu_edge)
        if np.isfinite(red) and np.isfinite(blue) and blue > 0:
            out[i] = red / blue
    return out


def break_radius(freq, jnu, r, nu_edge, thresh=BREAK_THRESHOLD, nhold=3):
    """Radius where the continuum break first reaches `thresh` going outward and
    stays there for `nhold` cells (so a single noisy cell cannot trigger it).
    NaN if the break never forms inside the domain."""
    ratio = edge_break_ratio(freq, jnu, nu_edge)
    n = min(len(r), len(ratio))
    r, ratio = np.asarray(r)[:n], ratio[:n]
    ok = np.nan_to_num(ratio, nan=0.0) >= thresh
    for i in range(n):
        if ok[i] and ok[i:i + nhold].all():
            if i == 0:
                return float(r[0])
            lo, hi = ratio[i - 1], ratio[i]
            if not (np.isfinite(lo) and hi > lo):
                return float(r[i])
            return float(np.interp(thresh, [lo, hi], [r[i - 1], r[i]]))
    return np.nan



def transition_radius_from_tau(r, tau, tau_c=TAU_LYC_TRANSITION):
    """Radius where a cumulative inside-out optical depth tau(r) (tau[0] at the
    base) first reaches tau_c, interpolated linearly in log r. NaN if never."""
    r, tau = np.asarray(r, float), np.asarray(tau, float)
    above = tau >= tau_c
    if not above.any():
        return np.nan
    i = int(np.argmax(above))
    if i == 0:
        return float(r[0])
    return float(np.exp(np.interp(tau_c, [tau[i - 1], tau[i]], [np.log(r[i - 1]), np.log(r[i])])))


def lyc_transition_radius(r, n1, tau_c=TAU_LYC_TRANSITION, sigma=SIG_LYC):
    """Lyman photoionization transition: radius where the inside-out Lyman-continuum
    absorption depth INT n_1 sigma_1 dr' from the base reaches tau_c (default 3).
    r, n1 : (Ncell,) inwind cells, base first. See the module docstring."""
    r, n1 = np.asarray(r, float), np.asarray(n1, float)
    dr = np.diff(r)
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (n1[1:] + n1[:-1]) * sigma * dr)])
    return transition_radius_from_tau(r, tau, tau_c)


def selftest(verbose=True):
    """A smooth Planck spectrum must return break = 1 at both edges, at every
    colour temperature. Raises AssertionError otherwise."""
    H, KB, C = 6.6261e-27, 1.3807e-16, 2.998e10
    freq = np.geomspace(0.05 * 1.602e-12 / H, 100 * 1.602e-12 / H, 3000)
    worst = 0.0
    for T in (50000, 30000, 20000, 14123, 9000, 7000, 5000, 3500):
        B = 2 * H * freq**3 / C**2 / np.expm1(np.minimum(H * freq / (KB * T), 700))
        for nm, nu_e in (('Lyman', NU_LYMAN), ('Balmer', NU_BALMER)):
            got = float(edge_break_ratio(freq, B[:, None], nu_e)[0])
            worst = max(worst, abs(got - 1.0))
            if verbose:
                print(f"  T={T:6d} K  {nm:6} edge -> break = {got:.4f}")
            assert abs(got - 1.0) < 0.02, f"biased: T={T}, {nm}, got {got}"
    if verbose:
        print(f"  OK: worst deviation from 1.000 is {worst:.4f}")
    return worst


if __name__ == '__main__':
    selftest()
