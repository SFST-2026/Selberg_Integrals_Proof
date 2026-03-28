#!/usr/bin/env python3
"""
Verify Dixon identity I0 = C(3k;k,k,k) = (3k)!/(k!)^3 for k=1..15.
Also verify Selberg => Dixon connection.
"""
from mpmath import *
mp.dps = 40

for k in range(1, 16):
    # Dixon sum
    dixon = sum((-1)**j * binomial(2*k, k+j)**3 for j in range(-k, k+1))
    trinomial = fac(3*k) / fac(k)**3
    # Numerical Selberg integral
    N = 40
    gl = []
    for i in range(1, N+1):
        x = cos(pi*(4*i-1)/(4*N+2))
        for _ in range(50):
            p0, p1 = mpf(1), x
            for j in range(2, N+1):
                p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
            dp = N*(p0-x*p1)/(1-x*x)
            x -= p1/dp
        gl.append(((x+1)/2, 2/((1-x*x)*dp*dp)/2))
    I0 = sum(w1*w2*power(64*sin(pi*a1)**2*sin(pi*a2)**2*sin(pi*(a1+a2))**2, k)
             for a1,w1 in gl for a2,w2 in gl)
    err = abs(I0 - trinomial)/trinomial
    print(f"k={k:>2d}: Dixon={int(dixon):>15d}, C(3k;k,k,k)={int(trinomial):>15d}, "
          f"GL rel.err={nstr(err,4)}")
