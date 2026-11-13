#!/usr/bin/env python3
"""
Reference solutions for verifying 3D-1D coupling (see BENCHMARKS3D1D.md).

    benchmarks3D1D.py rigid    closed form of the rigid segment benchmark,
                               checked against a direct ODE integration
    benchmarks3D1D.py mms      symbolic check of the Dobroserdova, Olshanskii
                               & Simakov (2016) manufactured solution, the
                               corrected 3D pressure and the source terms
    benchmarks3D1D.py csv [DIR]
                               the reference time series and profiles of the
                               rigid-segment benchmark as specified in
                               BENCHMARKS3D1D.md (section 1, "Specification"),
                               as CSV files in DIR (default: reference/)

Needs numpy and scipy (rigid, csv) and sympy (mms).
"""
import math
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


# --------------------------------------------------------------------------
# The rigid-segment benchmark as specified (BENCHMARKS3D1D.md, section 1):
# the reference a code is compared with, written as CSV
# --------------------------------------------------------------------------

class RigidSegment:
    """The benchmark's data, SI. The 1D vessels: the single-pulse vessel of
    Boileau et al. (2015), Table II, with the sqrtElastic tube law
    linearised at the diastolic area (c0 = 6.1721 m/s). The pipe: length l,
    cross-section A3D. The inflow: Boileau Eq. 14, Q = Qpeak
    exp(-1e4 (t - 0.05)^2), i.e. a Gaussian of standard deviation
    sigma = 1/sqrt(2e4) s peaking at t = 0.05 s, times scale."""

    rho = 1050.0
    R = 0.01                    # vessel radius at the diastolic pressure
    h = 1.5e-3                  # wall thickness
    E = 400e3                   # Young's modulus
    L1 = 1.0                    # proximal vessel (flow inlet to pipe)
    L2 = 1.0                    # distal vessel (pipe to absorbing end)
    l = 0.08                    # pipe length
    Qpeak = 1e-6
    t_peak = 0.05
    sigma = 1/math.sqrt(2e4)

    def __init__(self, A3D=None, scale=1.0):
        # The pipe's cross-section: the discrete patch area of the 3D mesh
        # (the tutorial's 32-sided polygon by default), or pi R^2 for a 1D
        # reduction of the pipe
        self.A3D = A3D if A3D is not None else 0.5*32*self.R**2*math.sin(2*math.pi/32)
        self.scale = scale
        self.A0 = math.pi*self.R**2
        beta = 4/3*math.sqrt(math.pi)*self.E*self.h
        self.c0 = math.sqrt(beta/(2*self.rho*self.A0))*self.A0**0.25
        self.Zc = self.rho*self.c0/self.A0
        self.L = self.rho*self.l/self.A3D            # inertance
        self.tau = self.L/(2*self.Zc)
        self.P = self.Zc*self.Qpeak*scale           # incident peak
        self.t_arrive = self.t_peak + self.L1/self.c0

    def pinc(self, T):
        return gaussian(T, self.P, self.sigma)

    def pt(self, T):
        return transmitted(T, self.P, self.sigma, self.tau)

    def pipe(self, t):
        """Q through the pipe, P at its inlet and outlet"""
        T = np.asarray(t, dtype=float) - self.t_arrive
        p = self.pt(T)
        return p/self.Zc, 2*self.pinc(T) - p, p

    def profile(self, vessel, x, t):
        """P, Q and A along a vessel ("proximal" or "distal", x from its
        upstream end) at time t: the incident, reflected and, at the flow
        inlet (Gamma = +1 in pressure), re-reflected waves"""
        x = np.asarray(x, dtype=float)
        c0, Zc = self.c0, self.Zc
        if vessel == "proximal":
            d = self.L1 - x
            pr = lambda T: self.pinc(T) - self.pt(T)
            fwd = self.pinc(t - x/c0 - self.t_peak) + pr(t - self.t_arrive - (self.L1 + x)/c0)
            bwd = pr(t - self.t_arrive - d/c0)
        else:
            fwd = self.pt(t - self.t_arrive - x/c0)
            bwd = np.zeros_like(x)
        P = fwd + bwd
        Q = (fwd - bwd)/Zc
        A = self.A0*(1 + P/(self.rho*c0**2))
        return P, Q, A

    def profile_avg(self, vessel, x, dx, t, n=5):
        """Cell averages of profile() over [x - dx/2, x + dx/2] by n-point
        Gauss-Legendre quadrature: what a finite-volume code's cell values
        should be compared with"""
        g, w = np.polynomial.legendre.leggauss(n)
        acc = None
        for gi, wi in zip(g, w):
            vals = self.profile(vessel, x + 0.5*dx*gi, t)
            acc = [wi*v for v in vals] if acc is None else [a + wi*v for a, v in zip(acc, vals)]
        return [a/2 for a in acc]


def writeCSV(outdir="reference"):
    import os
    os.makedirs(outdir, exist_ok=True)
    b = RigidSegment()
    hdr = (f"# rigid-segment benchmark (BENCHMARKS3D1D.md): rho {b.rho}, R {b.R}, "
           f"h {b.h}, E {b.E}, c0 {b.c0:.10g} m/s, Zc {b.Zc:.10g} Pa s/m3, "
           f"pipe l {b.l} m, A3D {b.A3D:.10g} m2 (32-gon), L {b.L:.10g} kg/m4, "
           f"tau {b.tau:.10g} s, inflow Qpeak {b.Qpeak} m3/s at t = {b.t_peak} s, "
           f"sigma {b.sigma:.10g} s, incident peak {b.P:.10g} Pa\n")
    t = np.round(np.arange(0, 0.35 + 1e-9, 1e-4), 6)
    Q, P1, P2 = b.pipe(t)
    with open(os.path.join(outdir, "rigidSegment_pipe.csv"), "w") as f:
        f.write(hdr + "t [s],Q [m3/s],P1 pipe inlet [Pa],P2 pipe outlet [Pa]\n")
        for row in zip(t, Q, P1, P2):
            f.write(",".join(f"{v:.12g}" for v in row) + "\n")
    for vessel, L in (("proximal", b.L1), ("distal", b.L2)):
        x = np.linspace(0, L, 2001)
        P, Qx, A = b.profile(vessel, x, 0.35)
        with open(os.path.join(outdir, f"rigidSegment_{vessel}_t0.35.csv"), "w") as f:
            f.write(hdr + "x [m],P [Pa],Q [m3/s],A [m2]   (point values at t = 0.35 s;"
                    " compare cell averages with profile_avg())\n")
            for row in zip(x, P, Qx, A):
                f.write(",".join(f"{v:.12g}" for v in row) + "\n")
    print(f"written to {outdir}/: rigidSegment_pipe.csv, "
          "rigidSegment_proximal_t0.35.csv, rigidSegment_distal_t0.35.csv")
    print(hdr.strip("# \n"))


if __name__ == "__main__":
    what = sys.argv[1:] or ["rigid", "mms"]
    if what[0] == "csv":
        writeCSV(*what[1:2])
        sys.exit()
    for w in what:
        {"rigid": checkRigid, "mms": checkMMS}[w]()
        print()
