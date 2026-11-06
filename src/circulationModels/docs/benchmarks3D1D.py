#!/usr/bin/env python3
"""
Reference solutions for verifying 3D-1D coupling (see BENCHMARKS3D1D.md).

    benchmarks3D1D.py rigid    closed form of the rigid segment benchmark,
                               checked against a direct ODE integration
    benchmarks3D1D.py mms      symbolic check of the Dobroserdova, Olshanskii
                               & Simakov (2016) manufactured solution, the
                               corrected 3D pressure and the source terms

Needs numpy and scipy (rigid) and sympy (mms).
"""
import sys

import numpy as np


# --------------------------------------------------------------------------
# Rigid segment: a rigid straight pipe of length l between two semi-infinite
# 1D lines with wave speed c. Linear, inviscid; the pipe is an inertance
# rho l / A in series, so with tau = l/(2c)
#     Gamma = i w tau/(1 + i w tau),    T = 1/(1 + i w tau),
#     tau dp_t/dt + p_t = p_inc,        p_r = p_inc - p_t
# --------------------------------------------------------------------------

def reflection(omega, tau):
    """Pressure reflection and transmission coefficients"""
    iwt = 1j*omega*tau
    return iwt/(1 + iwt), 1/(1 + iwt)


def gaussian(T, P, sigma):
    """Incident pulse at the interface, peak at T = 0"""
    return P*np.exp(-T**2/(2*sigma**2))


def transmitted(T, P, sigma, tau):
    """Transmitted pressure for the Gaussian incident pulse:
    p_t = P sigma/tau sqrt(pi/2) exp(-T^2/(2 sigma^2)) erfcx(a),
    a = (sigma^2/tau - T)/(sqrt(2) sigma); erfcx keeps it finite as
    tau -> 0, and the erfc form is used where a < 0"""
    from scipy.special import erfc, erfcx

    T = np.asarray(T, dtype=float)
    a = (sigma**2/tau - T)/(np.sqrt(2)*sigma)
    k = P*sigma/tau*np.sqrt(np.pi/2)
    pos = a >= 0
    out = np.empty_like(T)
    out[pos] = k*np.exp(-T[pos]**2/(2*sigma**2))*erfcx(a[pos])
    # Only where a < 0, so that the exponential is not evaluated (and does
    # not overflow) where it is not used
    neg = ~pos
    out[neg] = k*np.exp(sigma**2/(2*tau**2) - T[neg]/tau)*erfc(a[neg])
    return out


def checkRigid():
    P, sigma = 1.0, 5e-3
    print("rigid segment: Gaussian pulse, sigma = %g s" % sigma)
    for tau in (5e-4, 5e-3, 5e-2):
        # Reference: RK4 on tau dp/dt = p_inc - p from rest
        dt = min(tau, sigma)/2000
        T = np.arange(-10*sigma, 10*sigma + 20*tau, dt)

        def rhs(Ti, p):
            return (gaussian(Ti, P, sigma) - p)/tau

        p = np.zeros_like(T)
        for n in range(len(T) - 1):
            h = dt
            k1 = rhs(T[n], p[n])
            k2 = rhs(T[n] + h/2, p[n] + h/2*k1)
            k3 = rhs(T[n] + h/2, p[n] + h/2*k2)
            k4 = rhs(T[n] + h, p[n] + h*k3)
            p[n + 1] = p[n] + h/6*(k1 + 2*k2 + 2*k3 + k4)

        exact = transmitted(T, P, sigma, tau)
        pr = gaussian(T, P, sigma) - exact
        # Energy: incident = reflected + transmitted (same impedance)
        Ei = np.trapezoid(gaussian(T, P, sigma)**2, T)
        Ert = np.trapezoid(pr**2 + exact**2, T)
        print("  tau = %-7g max|p_t - ODE| = %.2e   peak p_t = %.4f"
              "   peak |p_r| = %.4f   energy (r+t)/i = %.6f"
              % (tau, np.max(np.abs(exact - p)), exact.max(),
                 np.abs(pr).max(), Ert/Ei))

    print("\n  |Gamma| for a 5 cm rigid segment, c = 5 m/s (tau = 5 ms):")
    for f in (1, 2, 5, 10, 20):
        G, _ = reflection(2*np.pi*f, 5e-3)
        print("    %2d Hz  %.3f" % (f, abs(G)))


# --------------------------------------------------------------------------
# Dobroserdova, Olshanskii & Simakov (2016), Section 4.1
# --------------------------------------------------------------------------

def dosMMS():
    """The manufactured solution as printed, with the corrected 3D
    pressure. Returns a dict of sympy expressions."""
    import sympy as sp

    t, x, y, z = sp.symbols("t x y z", real=True)
    al, be, c0, Sh, rho, nu = sp.symbols(
        "alpha beta c0 Shat rho nu", positive=True)

    c = sp.cos(2*be*sp.pi*t)
    s = sp.sin(2*be*sp.pi*t)
    E = sp.exp(1 - c)

    m = dict(t=t, x=x, y=y, z=z, alpha=al, beta=be, c0=c0, Shat=Sh,
             rho=rho, nu=nu)

    # 1D on x in [0, 5] (uniform in x), coupled at x = b = 5
    m["S"] = Sh*sp.exp(al*(c - 1))
    m["u1D"] = al*sp.exp(-al*(c - 1))*(1 - sp.exp(c - 1))
    m["pbar"] = c0**2*al*(c - 1)

    # 0D: I = R0 = 0, V0 = 0, C = 1/c0^2, time-dependent resistances
    m["C"] = 1/c0**2
    m["V"] = al*(c - 1) + al/c0**2*(E - 1)
    m["Vdot"] = 2*be*sp.pi*s*al*(E/c0**2 - 1)
    m["R1D0D"] = -E/Sh
    m["R0D3D"] = E/Sh

    # 3D: x in (-1, 1), y^2 + z^2 < 1, Gamma at x = -1
    X = Sh*al*(1 - sp.exp(c - 1)) - m["Vdot"]
    m["u3D"] = 2/sp.pi*X*(1 - y**2 - z**2)
    m["p3D_printed"] = (c0**2*al*(c - 1)
                        - 2*be*sp.pi/Sh*al*s*al*E*(E/c0**2 - 1)
                        + 10*(-1 - x))
    m["p3D"] = (c0**2*al*(c - 1)
                + 2*be*sp.pi/Sh*al*s*E*(E/c0**2 - 1)
                + 10*(-1 - x))
    return m


def checkMMS():
    import sympy as sp

    m = dosMMS()
    t, x, y, z = m["t"], m["x"], m["y"], m["z"]
    al, c0, Sh, rho, nu = m["alpha"], m["c0"], m["Shat"], m["rho"], m["nu"]
    S, u1, pbar = m["S"], m["u1D"], m["pbar"]
    p0D = m["V"]/m["C"]

    r, th = sp.symbols("r theta", nonnegative=True)
    Q1D = S*u1
    Q3D = sp.integrate(sp.integrate(
        m["u3D"].subs({y: r*sp.cos(th), z: r*sp.sin(th)})*r,
        (r, 0, 1)), (th, 0, 2*sp.pi))

    def check(name, e):
        print("  %-48s %s" % (name, "ok" if sp.simplify(e) == 0 else "FAILS"))

    print("Dobroserdova et al. (2016) MMS, Shat = pi (3D radius 1), rho = 1:")
    check("(2)  pbar = rho c0^2 ln(S/Shat)", pbar - c0**2*sp.log(S/Sh))
    check("     Vdot as printed = dV/dt", m["Vdot"] - sp.diff(m["V"], t))
    check("(7)  p0D = (V - V0)/C, I = R0 = V0 = 0", p0D - c0**2*m["V"])
    check("(8)  dV/dt = Q1D - Q3D", sp.diff(m["V"], t) - (Q1D - Q3D))
    check("(9)  pbar - p0D = R1D0D Q1D", pbar - p0D - m["R1D0D"]*Q1D)
    check("(10) p0D - p = R0D3D Q3D, p as printed",
          p0D - m["p3D_printed"].subs(x, -1) - m["R0D3D"]*Q3D)
    check("(10) p0D - p = R0D3D Q3D, p corrected",
          p0D - m["p3D"].subs(x, -1) - m["R0D3D"]*Q3D)
    check("     div u = 0", sp.diff(m["u3D"], x))
    check("     d u/dn = 0 on Gamma (viscous traction)",
          sp.diff(m["u3D"], x).subs(x, -1))

    print("\nSources the solution needs (rho = 1 in the paper):")
    phi = sp.diff(S, t) + sp.diff(S*u1, x)
    psi = sp.diff(u1, t) + sp.diff(u1**2/2 + pbar/rho, x)
    fx = (rho*sp.diff(m["u3D"], t)
          - nu*(sp.diff(m["u3D"], y, 2) + sp.diff(m["u3D"], z, 2))
          + sp.diff(m["p3D"], x))
    print("  phi  =", sp.simplify(phi))
    print("  psi  =", sp.simplify(psi))
    print("  f_x  =", fx)
    print("  outlet traction (x = 1): -p n, p =",
          sp.simplify(m["p3D"].subs(x, 1)))


if __name__ == "__main__":
    what = sys.argv[1:] or ["rigid", "mms"]
    for w in what:
        {"rigid": checkRigid, "mms": checkMMS}[w]()
        print()
