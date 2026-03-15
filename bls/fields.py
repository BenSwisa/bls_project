from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
from .utils import inv_mod
from .poly import Poly, poly_ext_gcd

# ---------- Prime field F_p ----------

@dataclass(frozen=True)
class Fp:
    p: int
    n: int  # 0..p-1

    def __post_init__(self):
        object.__setattr__(self, "n", self.n % self.p)

    @staticmethod
    def zero(p: int) -> "Fp":
        return Fp(p, 0)

    @staticmethod
    def one(p: int) -> "Fp":
        return Fp(p, 1)

    def __int__(self) -> int:
        return self.n

    def __repr__(self) -> str:
        return str(self.n)

    def _coerce(self, other) -> "Fp":
        if isinstance(other, Fp):
            assert self.p == other.p
            return other
        if isinstance(other, int):
            return Fp(self.p, other)
        raise TypeError(f"Unsupported operand type: {type(other)}")

    def __add__(self, other) -> "Fp":
        o = self._coerce(other)
        return Fp(self.p, self.n + o.n)

    def __radd__(self, other) -> "Fp":
        return self.__add__(other)

    def __sub__(self, other) -> "Fp":
        o = self._coerce(other)
        return Fp(self.p, self.n - o.n)

    def __rsub__(self, other) -> "Fp":
        o = self._coerce(other)
        return Fp(self.p, o.n - self.n)

    def __neg__(self) -> "Fp":
        return Fp(self.p, -self.n)

    def __mul__(self, other) -> "Fp":
        o = self._coerce(other)
        return Fp(self.p, self.n * o.n)

    def __rmul__(self, other) -> "Fp":
        return self.__mul__(other)

    def inv(self) -> "Fp":
        return Fp(self.p, inv_mod(self.n, self.p))

    def __truediv__(self, other) -> "Fp":
        o = self._coerce(other)
        return self * o.inv()

    def __rtruediv__(self, other) -> "Fp":
        o = self._coerce(other)
        return o * self.inv()

    def __pow__(self, e: int) -> "Fp":
        if e < 0:
            return (self ** (-e)).inv()
        return Fp(self.p, pow(self.n, e, self.p))

    def is_zero(self) -> bool:
        return self.n == 0

    def legendre(self) -> int:
        """Return 0,1,-1 for quadratic residue in F_p."""
        if self.n == 0:
            return 0
        ls = pow(self.n, (self.p - 1) // 2, self.p)
        return -1 if ls == self.p - 1 else ls

# ---------- Extension field F_{p^k} = F_p[x]/(modulus) ----------

@dataclass(frozen=True)
class Fpk:
    p: int
    modulus: Poly  # monic irreducible of degree k
    coeffs: Tuple[int, ...]  # reduced poly rep deg<k

    def __post_init__(self):
        assert self.modulus.p == self.p
        c = [a % self.p for a in self.coeffs]
        poly = Poly(tuple(c), self.p).mod(self.modulus)
        object.__setattr__(self, "coeffs", poly.coeffs)

    @property
    def k(self) -> int:
        return self.modulus.degree()

    @property
    def q(self) -> int:
        return self.p ** self.k

    @staticmethod
    def zero(p: int, modulus: Poly) -> "Fpk":
        return Fpk(p, modulus, ())

    @staticmethod
    def one(p: int, modulus: Poly) -> "Fpk":
        return Fpk(p, modulus, (1,))

    @staticmethod
    def from_int(p: int, modulus: Poly, n: int) -> "Fpk":
        return Fpk(p, modulus, (n % p,))

    def is_zero(self) -> bool:
        return len(self.coeffs) == 0

    def __repr__(self) -> str:
        if not self.coeffs:
            return "0"
        terms = []
        for i,c in enumerate(self.coeffs):
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                terms.append(f"{c}*alpha")
            else:
                terms.append(f"{c}*alpha^{i}")
        return " + ".join(terms) if terms else "0"

    def _poly(self) -> Poly:
        return Poly(self.coeffs, self.p)

    def _coerce(self, other) -> "Fpk":
        if isinstance(other, Fpk):
            assert self.p == other.p and self.modulus == other.modulus
            return other
        if isinstance(other, int):
            return Fpk.from_int(self.p, self.modulus, other)
        raise TypeError(f"Unsupported operand type: {type(other)}")

    def __add__(self, other) -> "Fpk":
        o = self._coerce(other)
        return Fpk(self.p, self.modulus, (self._poly() + o._poly()).coeffs)

    def __radd__(self, other) -> "Fpk":
        return self.__add__(other)

    def __sub__(self, other) -> "Fpk":
        o = self._coerce(other)
        return Fpk(self.p, self.modulus, (self._poly() - o._poly()).coeffs)

    def __rsub__(self, other) -> "Fpk":
        o = self._coerce(other)
        return Fpk(self.p, self.modulus, (o._poly() - self._poly()).coeffs)

    def __neg__(self) -> "Fpk":
        return Fpk(self.p, self.modulus, (-self._poly()).coeffs)

    def __mul__(self, other) -> "Fpk":
        o = self._coerce(other)
        return Fpk(self.p, self.modulus, (self._poly() * o._poly()).mod(self.modulus).coeffs)

    def __rmul__(self, other) -> "Fpk":
        return self.__mul__(other)

    def inv(self) -> "Fpk":
        if self.is_zero():
            raise ZeroDivisionError("inverse of 0")
        g, s, _ = poly_ext_gcd(self._poly(), self.modulus)
        if g.coeffs != (1,):
            raise ZeroDivisionError("element not invertible (should not happen in a field)")
        return Fpk(self.p, self.modulus, s.coeffs)

    def __truediv__(self, other) -> "Fpk":
        o = self._coerce(other)
        return self * o.inv()

    def __rtruediv__(self, other) -> "Fpk":
        o = self._coerce(other)
        return o * self.inv()

    def __pow__(self, e: int) -> "Fpk":
        if e < 0:
            return (self ** (-e)).inv()
        result = Fpk.one(self.p, self.modulus)
        b = self
        exp = e
        while exp > 0:
            if exp & 1:
                result = result * b
            exp >>= 1
            if exp:
                b = b * b
        return result

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Fpk):
            return False
        return self.p == other.p and self.modulus == other.modulus and self.coeffs == other.coeffs

# ---------- Generic sqrt for finite fields (Tonelli–Shanks) ----------

def field_order(x) -> int:
    if isinstance(x, Fp):
        return x.p
    if isinstance(x, Fpk):
        return x.q
    raise TypeError("unknown field element")

def field_one_like(x):
    if isinstance(x, Fp):
        return Fp.one(x.p)
    if isinstance(x, Fpk):
        return Fpk.one(x.p, x.modulus)
    raise TypeError("unknown field element")

def field_zero_like(x):
    if isinstance(x, Fp):
        return Fp.zero(x.p)
    if isinstance(x, Fpk):
        return Fpk.zero(x.p, x.modulus)
    raise TypeError("unknown field element")

def is_square(a) -> bool:
    q = field_order(a)
    if a == field_zero_like(a):
        return True
    ls = a ** ((q - 1) // 2)
    return ls == field_one_like(a)

def sqrt_ts(a):
    """Square root in F_q (q odd), using Tonelli–Shanks in the multiplicative group."""
    q = field_order(a)
    if q % 2 == 0:
        raise ValueError("field order must be odd")
    if a == field_zero_like(a):
        return field_zero_like(a)
    if not is_square(a):
        return None

    # Simple case: q == 3 mod 4
    if q % 4 == 3:
        return a ** ((q + 1) // 4)

    # Factor q-1 = Q * 2^S with Q odd
    Q = q - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    def rand_elem_like(a):
        import secrets
        if isinstance(a, Fp):
            return Fp(a.p, secrets.randbelow(a.p))
        if isinstance(a, Fpk):
            coeffs = tuple(secrets.randbelow(a.p) for _ in range(a.k))
            return Fpk(a.p, a.modulus, coeffs)
        raise TypeError

    z = rand_elem_like(a)
    while is_square(z):
        z = rand_elem_like(a)

    c = z ** Q
    x = a ** ((Q + 1) // 2)
    t = a ** Q
    m = S
    one = field_one_like(a)

    while t != one:
        i = 1
        t2i = t * t
        while i < m and t2i != one:
            t2i = t2i * t2i
            i += 1
        if i == m:
            return None
        b = c ** (1 << (m - i - 1))
        x = x * b
        t = t * (b * b)
        c = b * b
        m = i
    return x
