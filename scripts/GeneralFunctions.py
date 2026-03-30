

def prime_factors(N):
    """
    Returns the prime factorization of N as a dictionary: {prime: exponent}
    Example: 20 returns {2: 2, 5: 1}
    """
    factors = {}
    
    # 1. Handle the factor of 2 separately to optimize the main loop
    while N % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        N //= 2
        
    # 2. Check odd numbers starting from 3 up to the square root of N
    # We use int(N**0.5) + 1 to ensure we cover the exact square root if N is a perfect square
    i = 3
    while i * i <= N:
        while N % i == 0:
            factors[i] = factors.get(i, 0) + 1
            N //= i
        i += 2
        
    # 3. If N is still greater than 2, then N itself is a prime number
    if N > 2:
        factors[N] = factors.get(N, 0) + 1
        
    return factors