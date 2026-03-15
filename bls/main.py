from __future__ import annotations

from .find_generators import setup_system, sign, verify

def _prompt_int(msg: str) -> int:
    while True:
        s = input(msg).strip()
        try:
            return int(s)
        except ValueError:
            print("Please enter an integer.")

def main():
    print("=== BLS Signature (Reduced Tate Pairing) ===")
    p = _prompt_int("Enter prime p (must satisfy p≡3 mod4): ")
    A = _prompt_int("Enter curve parameter A: ")
    B = _prompt_int("Enter curve parameter B: ")
    priv = _prompt_int("Enter private key a (will be reduced mod r): ")
    m = input("Enter message m: ")

    params = setup_system(p, A, B)
    print("\n--- System parameters ---")
    print(f"p = {params.p}")
    print(f"Curve: y^2 = x^3 + {params.A} x + {params.B} over F_p")
    print(f"|E(F_p)| = {params.n}")
    print(f"r (largest prime factor of |E(F_p)|) = {params.r}")
    print(f"embedding degree k = ord_r(p) = {params.k}")
    print(f"irreducible modulus f(x) coefficients (low->high) = {params.f_modulus_coeffs}")
    print(f"P ∈ E(F_p) of order r: {params.P}")
    print(f"Q ∈ E(F_(p^k)) of order r: {params.Q}")

    Hm, sigma, pk = sign(params, priv, m)
    ok, left, right = verify(params, Hm, sigma, pk)

    print("\n--- BLS ---")
    print(f"Alice public key aQ: {pk}")
    print(f"H(m): {Hm}")
    print(f"Signature sigma = a*H(m): {sigma}")
    print(f"Pairing left  e(sigma, Q)    = {left}")
    print(f"Pairing right e(H(m), aQ)    = {right}")
    if ok:
        print(f"\nOK: The message '{m}' was received and verified: e(aH(m),Q) = e(H(m),aQ)")
    else:
        print(f"\nVerification failed.")

if __name__ == "__main__":
    main()
