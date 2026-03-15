from __future__ import annotations

from .ec import ECPoint, EllipticCurve

def _line_eval(curve: EllipticCurve, P: ECPoint, Q: ECPoint, R: ECPoint):
    """
    Evaluate l_{P,Q}(R). For vertical line, returns R.x - P.x.
    For regular line: (R.y - yP) - lambda*(R.x - xP)
    """
    if P.is_infinity() or Q.is_infinity() or R.is_infinity():
        raise ValueError("line evaluation with infinity not supported")

    x1, y1 = P.x, P.y
    x2, y2 = Q.x, Q.y
    xR, yR = R.x, R.y

    if P != Q:
        if x1 == x2:
            return xR - x1
        lam = (y2 - y1) / (x2 - x1)
    else:
        if y1 == 0 * y1:
            return xR - x1
        lam = (x1 * x1 * 3 + curve.A) / (y1 * 2)

    return (yR - y1) - lam * (xR - x1)

def _vertical_eval(R: ECPoint, Q: ECPoint):
    if R.is_infinity() or Q.is_infinity():
        raise ValueError("vertical evaluation with infinity not supported")
    return Q.x - R.x


def miller_function(curve: EllipticCurve, P: ECPoint, Q: ECPoint, r: int):
    if r <= 0:
        raise ValueError("r must be positive")
    if P.is_infinity() or Q.is_infinity():
        raise ValueError("P and Q must be finite")

    bits = bin(r)[2:]
    R = P

    one = Q.x * 0 + 1  # 1 in the underlying field
    f = one

    for bit in bits[1:]:
        # Doubling step: f <- f^2 * l_{R,R}(Q) / v_{2R}(Q)
        if R.is_infinity():
            # Should not happen for r odd except possibly at the very end, but keep safe:
            f = f * f
            R2 = R
            v = one
            l = one
        else:
            l = _line_eval(curve, R, R, Q)
            R2 = R + R
            v = one if R2.is_infinity() else _vertical_eval(R2, Q)
        if v == 0 * v:
            raise ZeroDivisionError("bad Q: Miller denominator is 0")
        f = (f * f) * (l / v)
        R = R2

        if bit == "1":
            # Addition step: f <- f * l_{R,P}(Q) / v_{R+P}(Q)
            l = _line_eval(curve, R, P, Q) if not R.is_infinity() else one
            R3 = R + P
            v = one if R3.is_infinity() else _vertical_eval(R3, Q)
            if v == 0 * v:
                raise ZeroDivisionError("bad Q: Miller denominator is 0")
            f = f * (l / v)
            R = R3
    return f


def reduced_tate_pairing(curve: EllipticCurve, P: ECPoint, Q: ECPoint, r: int, p: int, k: int):
    f = miller_function(curve, P, Q, r)
    exp = (p**k - 1) // r
    return f ** exp
