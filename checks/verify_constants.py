#!/usr/bin/env python3
"""Recompute the binary-face constants and window multipliers to 60 digits.

Python 3.11 or newer, no dependencies.
"""

from decimal import Decimal as D, getcontext

getcontext().prec = 75
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286")
ONE, TWO = D(1), D(2)


def entropy(rho):
    return -rho * rho.ln() - (1 - rho) * (1 - rho).ln() + rho * TWO.ln()


def stationarity(rho):
    dh = (TWO * (1 - rho) / rho).ln()
    return 2 * rho * dh - entropy(rho)


def main():
    lo, hi = D("0.01"), D("0.99")
    assert stationarity(lo) > 0 > stationarity(hi)
    for _ in range(210):
        mid = (lo + hi) / 2
        if stationarity(mid) > 0:
            lo = mid
        else:
            hi = mid
    rho = (lo + hi) / 2
    v = rho / 4
    c = ((3 * entropy(rho) / (rho / 2).sqrt()).ln() * TWO / 3).exp()
    r = (2 * v * c).sqrt()
    noiseless = v.sqrt() / (2 * (2 * PI).sqrt())
    noisy = noiseless * (4 * ONE.exp() / r) ** 2
    lower_gain_noiseless = c / 4 - D("1.23") * v / r
    lower_gain_noisy = c / 4 - D("3.68") * v / r
    assert hi - lo < D("1e-62")
    assert abs(stationarity(rho)) < D("1e-62")
    assert abs(c - D("3.4417465206384")) < D("1e-12")
    assert abs(noiseless - D("0.0602978769")) < D("1e-10")
    assert abs(noisy - D("11.333347057")) < D("1e-9")
    assert min(lower_gain_noiseless, lower_gain_noisy) > 0

    b = PI / 4
    gaussian_boundary = PI / 12 + (-b).exp() / (b * (1 - (-b).exp()))
    assert gaussian_boundary < D("1.34") < 3 * TWO.ln()
    for q in range(2, 16):
        mu = 1 - ONE / q
        w = mu * (1 - mu)
        c3 = (2 * mu - 1) / (6 * w**2)
        c4 = (1 - 3 * w) / (12 * w**3)
        assert abs(c3**2 * w / 2 - c4 / 5 + (1 + 2*w)/(360*w**3)) < D("1e-65")

    for label, value in (
        ("classical endpoint boundary-cost upper bound", gaussian_boundary),
        ("rho_star", rho), ("C", c), ("R_star", r),
        ("kappa_noiseless", noiseless), ("kappa_noisy", noisy),
        ("noiseless endpoint log-mean gain bound", lower_gain_noiseless),
        ("noisy endpoint log-mean gain bound", lower_gain_noisy),
    ):
        print(f"{label}: {value:.18f}")


if __name__ == "__main__":
    main()
