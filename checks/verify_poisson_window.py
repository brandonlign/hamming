#!/usr/bin/env python3
"""Recompute the vertex and edge Poisson-window constants, including binary orbits."""
from math import ceil, exp, floor, log, pi, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln, i0e


def kernel(h):
    return -log(i0e(exp(h) / (2*pi))) - max(h, 0)/2


def parameters(q, j0, j1):
    p = 1 - 1/q
    v = p*(1-p)
    a = (3*log(q)/sqrt(2*v))**(2/3)
    r = sqrt(2*v*a)
    g = 2 if q == 2 else 1
    t = (g*a/2 + r**5*(1+2*v)/(360*v**3) - 2*v*j0/r)/r
    u = g/(4*r)
    d = r**3/(6*v) + 2*v*v*j1/r**3
    return p, v, a, r, t, u, d


def centered_shell_cost(n, p, v, a, r, t, u, edge=False):
    s = n**(1/3)
    radii = np.arange(max(0, floor(n*p-2*r*s*s)),
                      min(n, ceil(n*p+2*r*s*s))+1)
    log_pmf = (gammaln(n+1) - gammaln(radii+1) - gammaln(n-radii+1)
               + radii*log(p) + (n-radii)*log(1-p))
    h = .5*log(2*pi*v*n) + a*s + t/s + u*log(n)/s**2 + log_pmf
    if edge:
        h += np.log(2-radii/(n*p))
    total = -np.log(i0e(np.exp(h)/(2*pi))).sum()
    g = 2 if p == .5 else 1
    return total + 2*n*log(1-p) - g*a*s/2 - (g/4+2*edge)*log(n)


def main():
    j0 = quad(kernel, -50, 50, epsabs=1e-12, points=[0])[0]
    j1 = quad(lambda h: h*kernel(h), -50, 50,
              epsabs=1e-12, points=[0])[0]
    assert abs(j0 + .1276510534987) < 1e-10
    assert abs(j1 + .8031368846847) < 1e-10
    print(f"J0={j0:.14f}, J1={j1:.14f}")
    for q in (2, 3, 5):
        p, v, a, r, t, u, d = parameters(q, j0, j1)
        mean = ((2*pi*v)**.5/16 if q == 2 else (2*pi*v)**.25/2)*exp(-d)
        if q == 2:
            assert abs(mean - .0431364146940171) < 1e-12
        print(f"q={q}: t*={t:.12f}, u*={u:.12f}, D={d:.12f}, "
              f"Lambda(0)={mean:.12f}")
        values = [centered_shell_cost(n, p, v, a, r, t, u)
                  for n in (2000, 16000, 128000, 1024000)]
        assert abs(values[-1]-d) < abs(values[0]-d)
        print(" centered shell costs:", ", ".join(f"{z:.8f}" for z in values))
        if q >= 3:
            delta = (q-3)*log(q)/(q-1)
            assert abs(delta-r**3*(3*p-2)/(6*v*p)) < 1e-13
            edge_mean = (q-1)**2/8 * (2*pi*v)**.25 * exp(-d-delta)
            edge_values = [centered_shell_cost(n, p, v, a, r, t, 9/(4*r), edge=True)
                           for n in (16000, 128000, 1024000)]
            assert abs(edge_values[-1]-d-delta) < abs(edge_values[0]-d-delta)
            print(f" edge Delta={delta:.12f}, Lambda_e(0)={edge_mean:.12f}")
            print(" centered edge costs:", ", ".join(f"{z:.8f}" for z in edge_values))



if __name__ == "__main__":
    main()
