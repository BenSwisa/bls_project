from __future__ import annotations

from .fields import Fp, sqrt_ts
from .ec import EllipticCurve, ECPoint

def message_to_x(m: str, p: int, encoding: str = "cp1255") -> int:
    b = m.encode(encoding, errors="strict")
    x = 0
    for byte in b:
        x = (x * 256 + byte) % p
    return x

def increment_and_try(curve: EllipticCurve, x0: int) -> ECPoint:
    assert curve.field_kind == "Fp"
    p = curve.p
    x = x0 % p
    while True:
        X = Fp(p, x)
        rhs = X*X*X + curve.A*X + curve.B
        y = sqrt_ts(rhs)
        if y is not None:
            return curve.point(X, y)
        x = (x + 1) % p

def hash_to_point(curve: EllipticCurve, m: str, cofactor: int, encoding: str = "cp1255") -> ECPoint:
    p = curve.p
    x0 = message_to_x(m, p, encoding=encoding)
    x = x0
    while True:
        Ptemp = increment_and_try(curve, x)
        H = cofactor * Ptemp
        if not H.is_infinity():
            return H
        x = (x + 1) % p
