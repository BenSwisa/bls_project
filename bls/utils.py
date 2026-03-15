from __future__ import annotations

from typing import List, Tuple, Dict
import math

def egcd_int(a: int, b: int) -> Tuple[int, int, int]:
    """Extended gcd: returns (g,x,y) with ax+by=g."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b != 0:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def inv_mod(a: int, p: int) -> int:
    a %= p
    if a == 0:
        raise ZeroDivisionError("inverse of 0")
    g, x, _ = egcd_int(a, p)
    if g != 1:
        raise ZeroDivisionError("no inverse exists")
    return x % p

def is_probable_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for <2^64-ish integers; fine for course sizes."""
    if n < 2:
        return False
    small_primes = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
    # write n-1 = d*2^s
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # bases good for <2^64
    for a in [2, 325, 9375, 28178, 450775, 9780504, 1795265022]:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n-1:
            continue
        for _ in range(s-1):
            x = (x * x) % n
            if x == n-1:
                break
        else:
            return False
    return True

def factorize(n: int) -> Dict[int, int]:
    """Trial division factorization; sufficient for small/medium n."""
    n0 = n
    f: Dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    # sanity
    prod = 1
    for p,e in f.items():
        prod *= p**e
    assert prod == n0
    return f

def largest_prime_factor(n: int) -> int:
    fac = factorize(n)
    return max(fac.keys())

def prime_factors(n: int) -> List[int]:
    return list(factorize(n).keys())

def multiplicative_order(a: int, mod: int) -> int:
    """Smallest k>0 such that a^k == 1 (mod mod). Assumes gcd(a,mod)=1."""
    a %= mod
    if math.gcd(a, mod) != 1:
        raise ValueError("a and mod not coprime")
    k = 1
    x = a % mod
    while x != 1:
        x = (x * a) % mod
        k += 1
        if k > mod:
            raise RuntimeError("order search failed")
    return k
