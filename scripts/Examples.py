'''# Finite Field and Polynomial Algebra Engine

## Overview
This is a pure Python, object-oriented mathematical engine designed to compute arithmetic over finite prime fields, polynomial rings, and quotient rings (extension fields). It handles the underlying algebraic structures required for cryptographic applications, such as Elliptic Curve Cryptography (ECC), without relying on heavy external math libraries.

The engine is built on three stacked layers of abstraction, allowing you to seamlessly move from basic modular arithmetic up to complex polynomial inversions.

## Core Components

### 1. `F_p` (Prime Field Elements)
Represents a scalar integer defined over a finite prime field. 
* Automatically reduces all arithmetic modulo the defined prime.
* Supports standard operators (`+`, `-`, `*`, `/`, `**`).
* **Division:** Implements true division by calculating the modular multiplicative inverse using the Extended Euclidean Algorithm.
* **Square Roots:** Includes a `sqrt()` method based on the Tonelli-Shanks algorithm / Euler's criterion.

### 2. `PolyF_p` (Polynomial Ring Elements)
Represents a polynomial where all coefficients are `F_p` objects.
* Initialized using a sparse dictionary mapping the power to its `F_p` coefficient (e.g., `{power: F_p_coefficient}`). Missing powers are mathematically treated as `0`.
* Supports polynomial addition, subtraction, and multiplication.
* **Division:** Implements polynomial long division over a field via `divmod()`, allowing the use of the `//` (quotient) and `%` (remainder) operators.

### 3. `PolyF_nQuotient` (Quotient Ring / Extension Field Elements)
Represents a polynomial evaluated modulo a defining irreducible polynomial.
* Initialized by providing a base `PolyF_p` polynomial and a quotient `PolyF_p` polynomial.
* Automatically reduces degrees after every arithmetic operation to ensure the polynomial remains within the bounds of the quotient ring.
* **Inverses & Division:** Provides an `inv()` method that uses the Extended Euclidean Algorithm to find the multiplicative inverse of the polynomial in the quotient ring. Division (`/`) is natively supported by multiplying by this inverse.
* **Exponentiation:** Uses a highly optimized "Double-and-Add" (exponentiation by squaring) algorithm that applies modulo reduction at every step, preventing memory/degree explosions when calculating cryptographically large powers.

---

## Important Architectural Notes
* **Modulus Safety:** There is a strict mathematical difference between `rmul` and `mul` in terms of the returning modulus. For example, in the elliptic curve equation `y^2 = x^3 + Ax + B`, the scalar coefficients `A` and `B` might belong to a different modulus than the underlying coordinates. Because the curve operates in the module of `x` and `y`, the returning modulus of cross-operations will strictly inherit the modulus of the coordinate field. 
* **Immutability:** Arithmetic operations return entirely new instances of the respective classes rather than mutating the objects in place.

---

## Quick Start Examples

### 1. Basic Field Arithmetic (`F_p`)
```python
# Define elements over F_7
a = F_p(3, 7)
b = F_p(6, 7)

# Addition automatically wraps around the modulus
print(a + b)  # Output: 2 (mod 7)

# Division calculates the inverse: 3 * (6^-1 mod 7) = 3 * 6 = 18 = 4 mod 7
print(a / b)  # Output: 4 (mod 7)

'''



# Polynomial Arithmetic (PolyF_p)
p = 7
A3 = F_p(3, p)

# Define P1(x) = 3x^2 + 3 over F_7
P1 = PolyF_p({0: A3, 2: A3})

# Define P2(x) = 3x^3 + 3 over F_7
P2 = PolyF_p({0: A3, 3: A3})

# Perform Polynomial Long Division
Quotient, Remainder = divmod(P2, P1)

print(f"Quotient: {Quotient}")
print(f"Remainder: {Remainder}")


# Quotient Rings & Extension Fields (PolyF_nQuotient)

p = 7
A1 = F_p(1, p)
A2 = F_p(2, p)

# 1. Define the irreducible quotient polynomial M(x) = x^2 + 1
M = PolyF_p({0: A1, 2: A1})

# 2. Define a working polynomial P(x) = x + 2
P = PolyF_p({0: A2, 1: A1})

# 3. Wrap it in the quotient ring: P(x) mod M(x)
P_ring = PolyF_nQuotient(P, M)

# 4. Calculate the multiplicative inverse using the Extended Euclidean Algorithm
# The math dictates the inverse of (x + 2) mod (x^2 + 1) over F_7 is (4x + 6)
P_inv = P_ring.inv()

print(f"Inverse: {P_inv}")

# Verification: Multiplying them together should result in exactly 1
print(f"Verification: {P_ring * P_inv}")