from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
from .fields import Fp, Fpk
from .ec import EllipticCurve, ECPoint, lift_point_to_extension
from .utils import largest_prime_factor, multiplicative_order
from .irreducible import find_irreducible_monic
from .hash_to_point import hash_to_point
from .pairing import reduced_tate_pairing

def count_points_fp(curve: EllipticCurve) -> int:
    """Count #E(F_p) by enumerating x and checking quadratic residues (O(p))."""
    assert curve.field_kind == "Fp"
    p = curve.p
    n = 1  # point at infinity
    for x in range(p):
        X = Fp(p, x)
        rhs = X*X*X + curve.A*X + curve.B
        ls = rhs.legendre()
        if ls == 0:
            n += 1
        elif ls == 1:
            n += 2
    return n

def find_point_of_order_r(curve: EllipticCurve, n: int, r: int) -> ECPoint:
    """Pick random R and take P = (n/r)R until P has order r."""
    assert n % r == 0
    h = n // r
    while True:
        R = curve.random_point()
        P = h * R
        if P.is_infinity():
            continue
        if (r * P).is_infinity():
            return P

def frobenius(point: ECPoint, p: int) -> ECPoint:
    """Frobenius: (x,y)->(x^p,y^p) on E(F_{p^k})."""
    if point.is_infinity():
        return point
    return ECPoint(point.curve, point.x ** p, point.y ** p)

def compute_E_Fpk_order(p: int, n_fp: int, k: int) -> int:
    """Compute #E(F_{p^k}) using trace recurrence from the lecture/handout."""
    t = p + 1 - n_fp
    a0, a1 = 2, t
    if k == 0:
        ak = a0
    elif k == 1:
        ak = a1
    else:
        a_prev2, a_prev1 = a0, a1
        for _ in range(2, k+1):
            a = t * a_prev1 - p * a_prev2
            a_prev2, a_prev1 = a_prev1, a
        ak = a_prev1
    return p**k + 1 - ak


def find_Q_of_order_r(curve_ext: EllipticCurve, p: int, n_fp: int, r: int, k: int, max_outer: int = 5000) -> ECPoint:
    """
    Generalized method (works also for the sample supersingular case):

      Nk = #E(F_{p^k})
      write Nk = r^s * c  with gcd(c,r)=1
      pick random T in E(F_{p^k})
      S = c*T
      Q = Frobenius(S) - S

    For ordinary curves usually s=1, so c = Nk/r and this matches the handout.
    """
    Nk = compute_E_Fpk_order(p, n_fp, k)

    # r-adic valuation s and cofactor c
    s = 0
    tmp = Nk
    while tmp % r == 0:
        tmp //= r
        s += 1
    c = Nk // (r ** s)

    for _ in range(max_outer):
        T = curve_ext.random_point()
        S = c * T
        if S.is_infinity():
            continue
        Q = frobenius(S, p) - S
        if Q.is_infinity():
            continue
        if not (r * Q).is_infinity():
            continue
        return Q

    raise RuntimeError("failed to find Q of order r")

@dataclass
class SystemParams:
    p: int
    A: int
    B: int
    n: int
    r: int
    k: int
    f_modulus_coeffs: Tuple[int, ...]
    P: ECPoint
    Q: ECPoint
    curve_fp: EllipticCurve
    curve_ext: EllipticCurve

def setup_system(p: int, A: int, B: int) -> SystemParams:
    if p <= 3 or p % 4 != 3:
        raise ValueError("Project requires prime p>3 with p ≡ 3 (mod 4).")
    A_fp = Fp(p, A)
    B_fp = Fp(p, B)
    curve_fp = EllipticCurve("Fp", p, A_fp, B_fp)

    disc = (A_fp*A_fp*A_fp) * Fp(p, 4) + (B_fp*B_fp) * Fp(p, 27)
    if disc.is_zero():
        raise ValueError("Curve is singular: discriminant is 0 mod p.")

    n = count_points_fp(curve_fp)
    r = largest_prime_factor(n)
    k = multiplicative_order(p, r)

    f = find_irreducible_monic(p, k)
    A_ext = Fpk.from_int(p, f, A % p)
    B_ext = Fpk.from_int(p, f, B % p)
    curve_ext = EllipticCurve("Fpk", p, A_ext, B_ext, modulus=f)

    P = find_point_of_order_r(curve_fp, n, r)
    Q = find_Q_of_order_r(curve_ext, p, n, r, k)

    return SystemParams(
        p=p, A=A%p, B=B%p, n=n, r=r, k=k,
        f_modulus_coeffs=f.coeffs,
        P=P, Q=Q,
        curve_fp=curve_fp, curve_ext=curve_ext
    )

def sign(params: SystemParams, priv_a: int, message: str, encoding: str = "cp1255"):
    a = priv_a % params.r
    if a == 0:
        raise ValueError("private key a must be non-zero mod r")
    cofactor = params.n // params.r
    Hm = hash_to_point(params.curve_fp, message, cofactor, encoding=encoding)
    sigma = a * Hm
    pk = a * params.Q
    return Hm, sigma, pk

def verify(params: SystemParams, Hm: ECPoint, sigma: ECPoint, pk: ECPoint):
    Hm_ext = lift_point_to_extension(Hm, params.curve_ext)
    sigma_ext = lift_point_to_extension(sigma, params.curve_ext)

    left = reduced_tate_pairing(params.curve_ext, sigma_ext, params.Q, params.r, params.p, params.k)
    right = reduced_tate_pairing(params.curve_ext, Hm_ext, pk, params.r, params.p, params.k)
    return left == right, left, right
