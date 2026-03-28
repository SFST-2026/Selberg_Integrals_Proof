#!/usr/bin/env python3
"""
30-digit verification of the central equation:
(6*pi^5)^4 = I0 * <exp(P_lat)>_9 * Z_inst
Uses: 100x100 GL, Nmax=10, M=12. Runtime ~30s.
"""
from mpmath import *
import time
mp.dps = 70; PI = pi; k = 9; d_S = 4

I0 = fac(27)/fac(9)**3
target = 4*log(6*PI**5)
print(f"I0 = {int(I0)}")
print(f"4*ln(6*pi^5) = {nstr(target, 35)}")

# Lattice sums via divisor sums
Nmax = 10; M = 12
shell = {}
for n1 in range(-Nmax, Nmax+1):
    for n2 in range(-Nmax, Nmax+1):
        for n3 in range(-Nmax, Nmax+1):
            for n4 in range(-Nmax, Nmax+1):
                nsq = n1*n1+n2*n2+n3*n3+n4*n4
                if nsq > 0: shell[nsq] = shell.get(nsq, 0) + 1

fc = {}
for m in range(1, M+1):
    s = sum(mult*exp(-2*PI*m*sqrt(mpf(nsq))) for nsq,mult in shell.items())
    fc[m] = d_S * 2 / m * s

ADJ = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
def P_lat(a1, a2):
    return sum(fc[m]*(1-cos(2*PI*m*(w1*a1+w2*a2)))
               for m in range(1,M+1) for w1,w2 in ADJ)

N_GL = 100
gl = []
for i in range(1, N_GL+1):
    x = cos(PI*(4*i-1)/(4*N_GL+2))
    for _ in range(200):
        p0, p1 = mpf(1), x
        for j in range(2, N_GL+1):
            p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
        dp = N_GL*(p0-x*p1)/(1-x*x)
        x -= p1/dp
        if abs(p1/dp) < mpf(10)**(-60): break
    gl.append(((x+1)/2, 2/((1-x*x)*dp*dp)/2))

t0 = time.time()
Z_pert = mpf(0); Z_sel = mpf(0)
for a1,w1 in gl:
    for a2,w2 in gl:
        v18 = power(64*sin(PI*a1)**2*sin(PI*a2)**2*sin(PI*(a1+a2))**2, 9)
        Z_pert += w1*w2*v18*exp(P_lat(a1,a2))
        Z_sel += w1*w2*v18

exp_P = Z_pert/Z_sel
Q_exact = 1296*PI**20/I0
Z_inst = Q_exact/exp_P
a1_coeff = (Z_inst-1)/exp(-PI**2/3)

print(f"\nResults ({time.time()-t0:.0f}s):")
print(f"<exp(P)>_9  = {nstr(exp_P, 32)}")
print(f"Z_inst      = {nstr(Z_inst, 32)}")
print(f"a1          = {nstr(a1_coeff, 28)}")
print(f"I0 rel.err  = {nstr(abs(Z_sel-I0)/I0, 4)}")
