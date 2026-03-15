from __future__ import annotations

import secrets
from .poly import Poly
from .utils import prime_factors

def random_monic_poly(p: int, deg: int) -> Poly:
    if deg < 1:
        raise ValueError("deg must be >= 1")
    coeffs = [secrets.randbelow(p) for _ in range(deg)]
    coeffs.append(1)  # monic
    return Poly(tuple(coeffs), p)

def rabin_is_irreducible(g: Poly) -> bool:
    """
    Rabin test (as in handout):
      1) g | x^{p^k} - x   <=>  x^{p^k} mod g == x
      2) for each prime divisor d of k: gcd(g, x^{p^{k/d}} - x) == 1
    """
    p = g.p
    k = g.degree()
    x = Poly.x(p)

    if Poly.pow_mod(x, p**k, g) != x.mod(g):
        return False

    for d in prime_factors(k):
        exp = p ** (k // d)
        h = Poly.pow_mod(x, exp, g) - x
        if Poly.gcd(g, h).degree() >= 1:
            return False
    return True

def find_irreducible_monic(p: int, deg: int, max_tries: int = 50000) -> Poly:
    for _ in range(max_tries):
        g = random_monic_poly(p, deg)
        if rabin_is_irreducible(g):
            return g
    raise RuntimeError("failed to find irreducible polynomial")
