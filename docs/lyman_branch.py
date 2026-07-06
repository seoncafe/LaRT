#!/usr/bin/env python3
"""Hydrogen radiative cascade after Lyman-series absorption.

Computes Einstein A coefficients A(nl -> n'l') from exact analytic hydrogenic
radial dipole integrals, then:
  - Lyman-series absorption oscillator strengths f(1s -> np)
  - direct branching ratios of each np upper level
  - full radiative cascade after 1s->np absorption:
      probability the de-excitation terminates via the 2s two-photon continuum
      vs. a direct np'->1s Lyman photon, and expected photons emitted per line.

Atomic units (a0 = 1, energy in Hartree).  The A-coefficient prefactor is
calibrated on the exactly-known A(2p->1s) and cross-checked on several others.
"""
import numpy as np
from scipy.special import genlaguerre, factorial
from scipy.integrate import quad

# ---- hydrogen radial wavefunction (atomic units, a0 = 1) ------------------
def R_nl(n, l):
    norm = np.sqrt((2.0/n)**3 * factorial(n-l-1) / (2.0*n*factorial(n+l)))
    lag = genlaguerre(n-l-1, 2*l+1)
    def f(r):
        rho = 2.0*r/n
        return norm * np.exp(-r/n) * rho**l * lag(rho)
    return f

def radial_integral(n, l, npr, lpr):
    """I = \int_0^inf R_nl R_n'l' r^3 dr  (in bohr)."""
    Ra, Rb = R_nl(n, l), R_nl(npr, lpr)
    integrand = lambda r: Ra(r)*Rb(r)*r**3
    rmax = 8.0*max(n, npr)**2
    val, _ = quad(integrand, 0.0, rmax, limit=400)
    return val

def dE_hartree(n, npr):
    """E_n - E_n'  (n upper, n' lower); positive for emission."""
    return 0.5*(1.0/npr**2 - 1.0/n**2)

# ---- oscillator strength (absorption, lower nl -> upper n'l') -------------
def f_absorption(nl, nu):
    (nlow, llow), (nup, lup) = nl, nu
    lmax = max(llow, lup)
    I = radial_integral(nup, lup, nlow, llow)
    dE = dE_hartree(nup, nlow)          # Hartree
    return (2.0/3.0)*dE*(lmax/(2.0*llow+1.0))*I**2

# ---- Einstein A (emission, upper nl -> lower n'l') ------------------------
# A = C * omega^3 * (lmax/(2l_up+1)) * I^2 ; omega in Hartree/hbar (a.u.).
# Calibrate C on A(2p->1s) = 6.2649e8 s^-1.
def A_kernel(n, l, npr, lpr):
    lmax = max(l, lpr)
    I = radial_integral(n, l, npr, lpr)
    w = dE_hartree(n, npr)              # a.u. angular frequency (E in Hartree, hbar=1)
    return w**3 * (lmax/(2.0*l+1.0)) * I**2

A_2p_1s_ref = 6.2649e8
C = A_2p_1s_ref / A_kernel(2, 1, 1, 0)

def A(n, l, npr, lpr):
    """Einstein A(nl -> n'l') in s^-1.  Requires n>n' and |l-l'|=1."""
    if npr >= n or abs(l-lpr) != 1 or lpr < 0:
        return 0.0
    return C * A_kernel(n, l, npr, lpr)

# =========================================================================
# 1. validation
# =========================================================================
print("=== validation (computed vs. literature, s^-1) ===")
checks = [
    ("A(2p->1s) Ly-a", A(2,1,1,0), 6.2649e8),
    ("A(3p->1s) Ly-b", A(3,1,1,0), 1.6725e8),
    ("A(3s->2p) Ha",   A(3,0,2,1), 6.3132e6),
    ("A(3p->2s) Ha",   A(3,1,2,0), 2.2449e7),
    ("A(3d->2p) Ha",   A(3,2,2,1), 6.4651e7),
    ("A(4p->1s) Ly-g", A(4,1,1,0), 6.8186e7),
    ("A(4d->2p) Hb",   A(4,2,2,1), 2.0625e7),
    ("A(4s->3p) Pa",   A(4,0,3,1), 1.8290e6),
    ("A(4f->3d) Pa",   A(4,3,3,2), 1.3785e7),
]
for name, comp, lit in checks:
    print(f"  {name:16s} {comp:.4e}  lit {lit:.4e}  ratio {comp/lit:.4f}")

print("\n=== Lyman-series absorption oscillator strengths f(1s->np) ===")
flit = {2:0.4162, 3:0.07910, 4:0.02899, 5:0.01394, 6:0.007799,
        7:0.004814, 8:0.003183, 9:0.002216, 10:0.001605}
for n in range(2, 11):
    fc = f_absorption((1,0),(n,1))
    print(f"  1s->{n}p  f = {fc:.5f}   lit {flit[n]:.5f}   ratio {fc/flit[n]:.4f}"
          f"   (Ly-a ratio {fc/f_absorption((1,0),(2,1)):.5f})")

# =========================================================================
# 2. direct branching ratios of each np upper level
# =========================================================================
def lower_channels(n, l):
    """allowed dipole lower states (n',l') with n'<n, l'=l+-1."""
    out = []
    for npr in range(1, n):
        for lpr in (l-1, l+1):
            if 0 <= lpr <= npr-1:
                out.append((npr, lpr))
    return out

def term(n, l):
    s = {0:'s',1:'p',2:'d',3:'f',4:'g',5:'h',6:'i',7:'k',8:'l',9:'m'}[l]
    return f"{n}{s}"

print("\n=== direct branching ratios of np levels ===")
for n in range(2, 11):
    ch = lower_channels(n, 1)
    Atot = sum(A(n,1,npr,lpr) for npr,lpr in ch)
    print(f"\n {term(n,1)}  (A_tot = {Atot:.4e} s^-1):")
    for npr,lpr in sorted(ch):
        a = A(n,1,npr,lpr)
        print(f"    -> {term(npr,lpr):3s}  A={a:.4e}  b={a/Atot:.5f}")

# =========================================================================
# 3. full radiative cascade after 1s -> np absorption
#    Terminal states: 1s (ground) and 2s (-> two-photon continuum).
#    Track: probability of each terminal + expected photons per line.
# =========================================================================
NMAX = 10
states = [(n,l) for n in range(1,NMAX+1) for l in range(0,n)]
TERMINAL = {(1,0), (2,0)}   # 1s ground; 2s -> two-photon

def branchings(n,l):
    ch = lower_channels(n,l)
    Atot = sum(A(n,l,npr,lpr) for npr,lpr in ch)
    return [((npr,lpr), A(n,l,npr,lpr)/Atot) for npr,lpr in ch]

from functools import lru_cache

# probability that cascade starting at (n,l) terminates in 2s (two-photon)
@lru_cache(maxsize=None)
def p_twophoton(n,l):
    if (n,l) == (2,0):   # 2s
        return 1.0
    if (n,l) == (1,0):   # 1s ground
        return 0.0
    p = 0.0
    for (npr,lpr), b in branchings(n,l):
        p += b * p_twophoton(npr,lpr)
    return p

# expected number of times a given line (nu->nl transition) is emitted,
# starting from (n0,l0); accumulate visitation probabilities.
@lru_cache(maxsize=None)
def visit_prob(n0,l0,n,l):
    """probability the cascade from (n0,l0) passes through state (n,l)."""
    if (n0,l0) == (n,l):
        return 1.0
    if (n0,l0) in TERMINAL:
        return 0.0
    p = 0.0
    for (npr,lpr), b in branchings(n0,l0):
        p += b * visit_prob(npr,lpr,n,l)
    return p

print("\n=== cascade after 1s->np absorption ===")
print(f"{'level':6s} {'P(2gamma via 2s)':>16s} {'P(direct Lyman)':>16s}")
for n in range(2,11):
    p2 = p_twophoton(n,1)
    print(f" {term(n,1):5s} {p2:16.5f} {1-p2:16.5f}")

# For Ly-b (3p) specifically, expected photons per line
print("\n=== Ly-beta (3p) cascade: expected photons per line ===")
n0,l0 = 3,1
lines = {}
for (n,l) in states:
    pv = visit_prob(n0,l0,n,l)
    if pv <= 0 or (n,l) in TERMINAL: 
        # still may emit from terminal? 2s emits two-photon (not a line); 1s none
        continue
    for (npr,lpr), b in branchings(n,l):
        lines[(n,l,npr,lpr)] = lines.get((n,l,npr,lpr),0) + pv*b
for k in sorted(lines, key=lambda t:(-lines[t])):
    n,l,npr,lpr = k
    print(f"  {term(n,l)}->{term(npr,lpr):3s}  <N> = {lines[k]:.5f}")
# two-photon terminus
print(f"  2s->1s two-photon  <N> = {visit_prob(n0,l0,2,0):.5f}")

# =========================================================================
# 4. full cascade inventory for every Lyman upper level  (data dump)
# =========================================================================
def cascade_lines(n0,l0):
    lines = {}
    for (n,l) in states:
        pv = visit_prob(n0,l0,n,l)
        if pv <= 0 or (n,l) in TERMINAL:
            continue
        for (npr,lpr), b in branchings(n,l):
            lines[(n,l,npr,lpr)] = lines.get((n,l,npr,lpr),0)+pv*b
    return lines

print("\n\n############ CASCADE INVENTORY (expected photons per line) ############")
for n0 in range(2,11):
    print(f"\n#### start = {term(n0,1)}  (Ly-{ {2:'alpha',3:'beta',4:'gamma',5:'delta',6:'epsilon',7:'6',8:'7',9:'8',10:'9'}[n0] })")
    lines = cascade_lines(n0,1)
    tot = 0.0
    for k in sorted(lines, key=lambda t:(-lines[t])):
        n,l,npr,lpr = k
        print(f"    {term(n,l)}->{term(npr,lpr):3s}  <N>={lines[k]:.6f}")
        tot += lines[k]
    p2 = visit_prob(n0,1,2,0)
    print(f"    2s->1s two-photon  <N>={p2:.6f}")
    print(f"    [check] total line photons+2gamma end = {tot:.4f};  P(2gamma)={p2:.5f}")
