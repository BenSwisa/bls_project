# BLS Signature (Reduced Tate Pairing) — Course Project

This project implements (educational) BLS signatures over an elliptic curve **y^2 = x^3 + A x + B** over a prime field **F_p**
together with the **reduced Tate pairing** over an extension field **F_{p^k}**.

It follows the course handout requirements:
- prime field arithmetic + fast exponentiation + inverse by extended Euclid
- elliptic curve group law + scalar multiplication (double-and-add)
- compute |E(F_p)|
- choose **r** = largest prime factor of |E(F_p)|
- compute embedding degree **k = ord_r(p)**
- build extension field **F_{p^k} = F_p[x]/(f(x))** using a random monic irreducible polynomial (Rabin test)
- hash-to-curve by "increment-and-try", then clear cofactor to land in the r-subgroup
- Miller's algorithm for f_{r,P}(Q)
- Reduced Tate pairing: e_r(P,Q) = f_{r,P}(Q)^((p^k-1)/r)
- BLS sign: sigma = a * H(m), public key: aQ, verify: e_r(sigma,Q) == e_r(H(m), aQ)

## Quick start

```bash
python -m bls.main
```

You will be prompted for:
- prime p (must satisfy p ≡ 3 (mod 4) and p>3)
- A, B (curve parameters mod p)
- private key a (used mod r)
- message m

The program prints the chosen parameters, the public key, the signature, and
a verification line: `e(sigma,Q) == e(H(m), aQ)`.

## Notes

- This is **not production cryptography** (no constant-time, no secure curves, no subgroup checks beyond what the assignment asks).
- Intended for small p (like the example p=103) for reasonable runtime.
