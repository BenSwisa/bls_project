from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple
from .utils import inv_mod

def _trim(a: List[int]) -> List[int]:
    while len(a) > 0 and a[-1] == 0:
        a.pop()
    return a

@dataclass(frozen=True)
class Poly:
    """Polynomial over F_p represented by coeff list [c0,c1,...] (low to high)."""
    coeffs: Tuple[int, ...]
    p: int

    @staticmethod
    def zero(p: int) -> "Poly":
        return Poly((), p)

    @staticmethod
    def one(p: int) -> "Poly":
        return Poly((1 % p,), p)

    @staticmethod
    def x(p: int) -> "Poly":
        return Poly((0, 1), p)

    def degree(self) -> int:
        return len(self.coeffs) - 1

    def lc(self) -> int:
        return self.coeffs[-1] if self.coeffs else 0

    def __repr__(self) -> str:
        if not self.coeffs:
            return "0"
        terms = []
        for i,c in enumerate(self.coeffs):
            c %= self.p
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                terms.append(f"{c}*x")
            else:
                terms.append(f"{c}*x^{i}")
        return " + ".join(terms) if terms else "0"

    def __add__(self, other: "Poly") -> "Poly":
        assert self.p == other.p
        n = max(len(self.coeffs), len(other.coeffs))
        out = [0] * n
        for i in range(n):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = other.coeffs[i] if i < len(other.coeffs) else 0
            out[i] = (a + b) % self.p
        return Poly(tuple(_trim(out)), self.p)

    def __sub__(self, other: "Poly") -> "Poly":
        assert self.p == other.p
        n = max(len(self.coeffs), len(other.coeffs))
        out = [0] * n
        for i in range(n):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = other.coeffs[i] if i < len(other.coeffs) else 0
            out[i] = (a - b) % self.p
        return Poly(tuple(_trim(out)), self.p)

    def __neg__(self) -> "Poly":
        out = [(-c) % self.p for c in self.coeffs]
        return Poly(tuple(_trim(out)), self.p)

    def __mul__(self, other: "Poly") -> "Poly":
        assert self.p == other.p
        if not self.coeffs or not other.coeffs:
            return Poly.zero(self.p)
        out = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i,a in enumerate(self.coeffs):
            a %= self.p
            if a == 0:
                continue
            for j,b in enumerate(other.coeffs):
                b %= self.p
                if b == 0:
                    continue
                out[i+j] = (out[i+j] + a*b) % self.p
        return Poly(tuple(_trim(out)), self.p)

    def monic(self) -> "Poly":
        if not self.coeffs:
            return self
        inv = inv_mod(self.lc(), self.p)
        return self.scale(inv)

    def scale(self, c: int) -> "Poly":
        c %= self.p
        out = [(c * a) % self.p for a in self.coeffs]
        return Poly(tuple(_trim(out)), self.p)

    def divmod(self, other: "Poly") -> Tuple["Poly", "Poly"]:
        """Polynomial long division over F_p. Returns (q,r)."""
        assert self.p == other.p
        p = self.p
        if not other.coeffs:
            raise ZeroDivisionError("poly div by 0")
        a = list(self.coeffs)
        b = list(other.coeffs)
        _trim(a); _trim(b)
        if not a:
            return Poly.zero(p), Poly.zero(p)
        deg_a = len(a) - 1
        deg_b = len(b) - 1
        if deg_a < deg_b:
            return Poly.zero(p), Poly(tuple(a), p)
        inv_lc_b = inv_mod(b[-1], p)
        q = [0] * (deg_a - deg_b + 1)
        r = a[:]
        for k in range(deg_a - deg_b, -1, -1):
            if len(r) - 1 < deg_b + k:
                continue
            coeff = (r[deg_b + k] * inv_lc_b) % p
            q[k] = coeff
            if coeff != 0:
                for j in range(deg_b + 1):
                    r[j + k] = (r[j + k] - coeff * b[j]) % p
        return Poly(tuple(_trim(q)), p), Poly(tuple(_trim(r)), p)

    def mod(self, m: "Poly") -> "Poly":
        return self.divmod(m)[1]

    @staticmethod
    def gcd(a: "Poly", b: "Poly") -> "Poly":
        assert a.p == b.p
        p = a.p
        x, y = a, b
        while y.coeffs:
            x, y = y, x.mod(y)
        return x.monic() if x.coeffs else Poly.zero(p)

    @staticmethod
    def pow_mod(base: "Poly", exp: int, mod_poly: "Poly") -> "Poly":
        """Fast exponentiation of polys modulo mod_poly."""
        assert base.p == mod_poly.p
        p = base.p
        result = Poly.one(p)
        b = base.mod(mod_poly)
        e = exp
        while e > 0:
            if e & 1:
                result = (result * b).mod(mod_poly)
            e >>= 1
            if e:
                b = (b * b).mod(mod_poly)
        return result

def poly_ext_gcd(a: Poly, b: Poly) -> Tuple[Poly, Poly, Poly]:
    """Extended gcd for polynomials: returns (g, s, t) with s*a + t*b = g."""
    assert a.p == b.p
    p = a.p
    r0, r1 = a, b
    s0, s1 = Poly.one(p), Poly.zero(p)
    t0, t1 = Poly.zero(p), Poly.one(p)
    while r1.coeffs:
        q, r2 = r0.divmod(r1)
        r0, r1 = r1, r2
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if r0.coeffs:
        inv = inv_mod(r0.lc(), p)
        r0 = r0.scale(inv)
        s0 = s0.scale(inv)
        t0 = t0.scale(inv)
    return r0, s0, t0
