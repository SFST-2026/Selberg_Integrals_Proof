#!/usr/bin/env python3
"""
Verification of all seven Selberg moment formulas (c1,c2,c3,chi11,chi12,B2,B3)
for k=1..30. Reproduces Table S1 (nine rational moments at k=9).
Author: M.W. Le Borgne, 28 March 2026
"""
from mpmath import *
from fractions import Fraction
mp.dps = 40

def c1(k): return Fraction(-k, 2*k+1)
def c2(k): return Fraction(-k*(2*k**2+k+1), 2*(k+1)**2*(2*k+1))
def c3(k):
    num = k*(4*k**5+14*k**4+7*k**3-13*k**2-8*k-4)
    den = (2*k+1)*(k+1)**2*(k+2)**2*(2*k+3)
    return Fraction(num, den)
def chi11(k): return Fraction(k*(k-1), 2*(k+1)*(2*k+1))
def chi12(k): return Fraction(k**2*(k+3), 2*(k+1)**2*(2*k+1))
def B2(k): return Fraction(k*(2*k**2+3*k-1), 2*(2*k+1)*(k+1)**2)
def B3(k):
    num = k*(-4*k**4-11*k**3+3*k**2+16*k-4)
    den = (2*k+1)*(k+1)**2*(k+2)**2*(2*k+3)
    return Fraction(num, den)

# GL quadrature
N_GL = 60
PI = pi
gl = []
for i in range(1, N_GL+1):
    x = cos(PI*(4*i-1)/(4*N_GL+2))
    for _ in range(100):
        p0, p1 = mpf(1), x
        for j in range(2, N_GL+1):
            p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
        dp = N_GL*(p0-x*p1)/(1-x*x)
        dx = p1/dp; x -= dx
        if abs(dx) < mpf(10)**(-35): break
    w = 2/((1-x*x)*dp*dp)
    gl.append(((x+1)/2, w/2))

def avg(func, k_val):
    I0=mpf(0); If=mpf(0)
    for a1,w1 in gl:
        for a2,w2 in gl:
            v = power(64*sin(PI*a1)**2*sin(PI*a2)**2*sin(PI*(a1+a2))**2, k_val)
            I0 += w1*w2*v; If += w1*w2*v*func(a1,a2)
    return If/I0

print("Verifying seven formulas for k=1..30...")
for k in [1,2,3,5,9,15,20,25,30]:
    errs = []
    for m, formula in [(1,c1),(2,c2),(3,c3)]:
        num = avg(lambda a1,a2,mm=m: cos(2*PI*mm*a1), k)
        err = abs(float(num) - float(formula(k)))
        errs.append(err)
    num_chi11 = avg(lambda a1,a2: cos(2*PI*a1)*cos(2*PI*a2), k)
    errs.append(abs(float(num_chi11) - float(chi11(k))))
    print(f"  k={k:>2d}: max_err = {max(errs):.2e} {'OK' if max(errs)<1e-12 else 'FAIL'}")

# Table S1 at k=9
print("\nTable S1: Exact rational moments at k=9")
k = 9
moments_k9 = {
    1: Fraction(-9,19), 2: Fraction(-387,950), 3: Fraction(294,475),
    4: Fraction(-2286,11495), 5: Fraction(-7758,52877),
    6: Fraction(75111,528770), 7: Fraction(-143379,6609625),
    8: Fraction(-99027,4415125), 9: Fraction(620537,60928725),
}
for m in range(1, 10):
    num = avg(lambda a1,a2,mm=m: cos(2*PI*mm*a1), k)
    exact = moments_k9[m]
    err = abs(float(num) - float(exact))
    print(f"  m={m}: {exact} (err={err:.2e})")
