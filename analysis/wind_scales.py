"""Spherization radius, in one place, so every figure normalises radii the same way.

r_sph is where the vertical radiation force on electrons balances gravity in a
thin disc: kappa_es D(r)/c = Omega^2 H with D = (3/8pi) Mdot Omega^2, giving

    r_sph = 3 kappa_es Mdot / (8 pi c h)  =  Mdot G M_BH / ((2h/3) L_Edd)

which for h = 0.3 is the  Mdot G M_BH / (0.2 L_Edd)  used to set the heating
radius R of the photon-tired models (winds/generate_model_1d_f.py), so for those
runs r/r_sph = 1 really is the base of the wind.

NOTE this is Mdot-DEPENDENT: r_sph = 2.844e14 cm * (Mdot / [Msun/yr]), i.e.
1.42e15 cm at Mdot = 5. FigHeatCool_TeXi and FigHaHbNetEmission previously
divided every run by the single constant 1.6e15 cm, which is r_sph only at
Mdot = 5 and distorts the other panels by up to 3x.

CONVENTION: the manuscript's Eq. 2 currently quotes 1.6e15 cm (Mdot/5), which
corresponds to kappa_es = 0.383 at h = 0.3, or to h = 0.267 at kappa_es = 0.34.
If that prefactor is kept rather than corrected, set KAPPA_ES = 0.3826 below and
every figure follows; nothing else needs to change.
"""
import numpy as np

G_CGS   = 6.674e-8
MSUN    = 1.989e33
YR      = 3.156e7
C_CGS   = 2.998e10
M_BH    = 1e6 * MSUN

KAPPA_ES = 0.34      # cm^2/g, electron scattering
H_OVER_R = 0.30      # disc aspect ratio at the spherization radius

L_EDD = 4 * np.pi * G_CGS * M_BH * C_CGS / KAPPA_ES


def r_sph(Mdot_Msun_per_yr):
    """Spherization radius [cm] for a mass-loss rate in Msun/yr."""
    Mdot = np.asarray(Mdot_Msun_per_yr, dtype=float) * MSUN / YR
    return Mdot * G_CGS * M_BH / ((2.0 * H_OVER_R / 3.0) * L_EDD)


def mdot_from_key(key):
    """'m2.5' / 'MU_m10' / 'pt_m15' -> 2.5 / 10.0 / 15.0 (Msun/yr)."""
    return float(key.rsplit('m', 1)[-1])


if __name__ == '__main__':
    print(f"kappa_es = {KAPPA_ES},  h = {H_OVER_R},  L_Edd = {L_EDD:.4e} erg/s")
    for md in (2.5, 5.0, 10.0, 15.0):
        print(f"  Mdot = {md:4} Msun/yr  ->  r_sph = {float(r_sph(md)):.4e} cm"
              f"   (old flat normalisation 1.6e15 was off by {1.6e15/float(r_sph(md)):.2f}x)")
