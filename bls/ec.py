from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import secrets

from .fields import Fp, Fpk, sqrt_ts, field_zero_like

FieldElt = object  # Fp or Fpk

@dataclass(frozen=True)
class ECPoint:
    curve: "EllipticCurve"
    x: Optional[FieldElt]
    y: Optional[FieldElt]

    def is_infinity(self) -> bool:
        return self.x is None and self.y is None

    def __neg__(self) -> "ECPoint":
        if self.is_infinity():
            return self
        return ECPoint(self.curve, self.x, -self.y)

    def __add__(self, other: "ECPoint") -> "ECPoint":
        return self.curve.add(self, other)

    def __sub__(self, other: "ECPoint") -> "ECPoint":
        return self + (-other)

    def __rmul__(self, k: int) -> "ECPoint":
        return self.curve.mul(self, k)

    def __repr__(self) -> str:
        if self.is_infinity():
            return "O"
        return f"({self.x}, {self.y})"

@dataclass(frozen=True)
class EllipticCurve:
    field_kind: str  # "Fp" or "Fpk"
    p: int
    A: FieldElt
    B: FieldElt
    modulus: Optional[object] = None  # Poly for Fpk

    def infinity(self) -> ECPoint:
        return ECPoint(self, None, None)

    def point(self, x: FieldElt, y: FieldElt) -> ECPoint:
        P = ECPoint(self, x, y)
        if not self.is_on_curve(P):
            raise ValueError("point is not on curve")
        return P

    def is_on_curve(self, P: ECPoint) -> bool:
        if P.is_infinity():
            return True
        x, y = P.x, P.y
        return (y * y) == (x * x * x + self.A * x + self.B)

    def add(self, P: ECPoint, Q: ECPoint) -> ECPoint:
        if P.is_infinity():
            return Q
        if Q.is_infinity():
            return P

        x1, y1 = P.x, P.y
        x2, y2 = Q.x, Q.y

        if x1 == x2 and y1 == -y2:
            return self.infinity()

        if P != Q:
            if x1 == x2:
                return self.infinity()
            lam = (y2 - y1) / (x2 - x1)
        else:
            if y1 == field_zero_like(y1):
                return self.infinity()
            lam = (x1 * x1 * 3 + self.A) / (y1 * 2)

        x3 = lam * lam - x1 - x2
        y3 = lam * (x1 - x3) - y1
        return ECPoint(self, x3, y3)

    def mul(self, P: ECPoint, k: int) -> ECPoint:
        if k < 0:
            return self.mul(-P, -k)
        R = self.infinity()
        addend = P
        n = k
        while n > 0:
            if n & 1:
                R = self.add(R, addend)
            n >>= 1
            if n:
                addend = self.add(addend, addend)
        return R

    def random_field_element(self):
        if self.field_kind == "Fp":
            return Fp(self.p, secrets.randbelow(self.p))
        else:
            assert self.modulus is not None
            coeffs = tuple(secrets.randbelow(self.p) for _ in range(self.modulus.degree()))
            return Fpk(self.p, self.modulus, coeffs)

    def sqrt(self, a):
        return sqrt_ts(a)

    def random_point(self, max_tries: int = 5000) -> ECPoint:
        for _ in range(max_tries):
            x = self.random_field_element()
            rhs = x * x * x + self.A * x + self.B
            y = self.sqrt(rhs)
            if y is None:
                continue
            return ECPoint(self, x, y)
        raise RuntimeError("failed to sample random point (too many tries)")

def embed_fp_to_fpk(x: Fp, modulus) -> Fpk:
    return Fpk.from_int(x.p, modulus, int(x))

def lift_point_to_extension(P: ECPoint, curve_ext: EllipticCurve) -> ECPoint:
    if P.is_infinity():
        return curve_ext.infinity()
    x = embed_fp_to_fpk(P.x, curve_ext.modulus)
    y = embed_fp_to_fpk(P.y, curve_ext.modulus)
    return ECPoint(curve_ext, x, y)
