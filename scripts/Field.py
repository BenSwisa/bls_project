

class F_p:
    def __init__(self,n: int,prime_modulus=7):
        self.prime_modulus=prime_modulus
        self._value=n

    
    def __add__(self,number):
        return self.__class__((self._value+number)%self.prime_modulus)
    
    def __sub__(self,number):
        return self.__class__((self._value-number)%self.prime_modulus)
        
    def __mul__(self,number):
        return self.__class__((self._value*number)%self.prime_modulus)

    @staticmethod
    def _mod_inverse(a, m):
        r0, r1 = m, a
        s0, s1 = 1, 0
        t0, t1 = 0, 1

        while r1 != 0:
            q = r0 // r1
            r0, r1 = r1, r0 - q*r1
            s0, s1 = s1, s0 - q*s1
            t0, t1 = t1, t0 - q*t1

        return t0 % m
    
    def __truediv__(self,number):
        return self.__class__((self._value*self._mod_inverse(number,self.prime_modulus))%self.prime_modulos)

    def __repr__(self):
        return f"{self._value}"
