'''

NOTE: there is a difference between rmul and mul in terms of return mod 


'''


class F_p:
    def __init__(self,n,prime_modulus):
        assert isinstance(n,int) , "Unsopported input type for field"
        assert prime_modulus>3 , "Unsopported prime modulus<=3"
        self.prime_modulus=prime_modulus
        self._value=n%self.prime_modulus
    
    def __add__(self,number):
        if isinstance(number,int):
            return self.__class__((self._value+number),self.prime_modulus)
        if isinstance(number,self.__class__):
            assert(number.prime_modulus==self.prime_modulus),"Must add with same field type"
            return self.__class__((self._value+number._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for add with F_p")


    def __radd__(self,number):
        return self.__add__(number)
    

    def __sub__(self,number):
        if isinstance(number,int):
            return self.__class__((self._value-number),self.prime_modulus)
        if isinstance(number,self.__class__):
            assert(number.prime_modulus==self.prime_modulus),"Must sub with same field type"
            return self.__class__((self._value-number._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for sub with F_p")


    def __rsub__(self,number):
        if isinstance(number,int):
            return self.__class__((number-self._value),self.prime_modulus)
        if isinstance(number,self.__class__):
            assert(number.prime_modulus==self.prime_modulus),"Must rsub with same field type"
            return self.__class__((number._value-self._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for rsub with F_p")

    def __neg__(self):
        return self.__class__((-1*self._value),self.prime_modulus)

    def __mul__(self,number):
        if isinstance(number,int):        
            return self.__class__((self._value*number),self.prime_modulus)
        if isinstance(number,self.__class__):
            return self.__class__((self._value*number._value),number.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for mul with F_p")


    def __rmul__(self,number):
        if isinstance(number,int):        
            return self.__class__((self._value*number),self.prime_modulus)
        if isinstance(number,self.__class__):
            return self.__class__((self._value*number._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for rmul with F_p")


    def __pow__(self,number):
        if isinstance(number,int):
            assert not(abs(number)<1 and abs(number)>0) , "Square root is not supported"         
            if number<0: 
                self._value=self._mod_inverse(self._value,self.prime_modulus)
                number*=-1
            return self.__class__((self._value**number),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for pow with F_p")


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
        if isinstance(number,int):  
            number=self._mod_inverse(number,self.prime_modulus)      
        elif isinstance(number,self.__class__):
            number=self._mod_inverse(number._value,self.prime_modulus)      
        else : raise ValueError(f"type {type(number) } unsuported for div with F_p")

        return number*self


    def sqrt(self):
        result= self._sqrt_fp(self._value,self.prime_modulus)
        assert result is not None , "No sqrt result over this field" 
        return self.__class__(result,self.prime_modulus)

    @staticmethod
    def _sqrt_fp(a, p):
        """Return x such that x^2 = a mod p (p must be an odd prime).
        Returns None if no square root exists."""
        
        a %= p
        if a == 0:
            return 0

        # Euler criterion
        if pow(a, (p - 1) // 2, p) != 1:
            return None

        # Fast case: p ≡ 3 mod 4
        if p % 4 == 3:
            return pow(a, (p + 1) // 4, p)

        # Factor p-1 = q * 2^s
        q = p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1

        # Find quadratic non-residue
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1

        c = pow(z, q, p)
        x = pow(a, (q + 1) // 2, p)
        t = pow(a, q, p)
        m = s

        while t != 1:
            i = 1
            temp = pow(t, 2, p)
            while temp != 1:
                temp = pow(temp, 2, p)
                i += 1

            b = pow(c, 2 ** (m - i - 1), p)

            x = (x * b) % p
            t = (t * b * b) % p
            c = (b * b) % p
            m = i

        return x        





    def __eq__(self, number):
        if isinstance(number,int):  
            return self._value==number
        elif isinstance(number,self.__class__):
            assert(number.prime_modulus==self.prime_modulus),"Must check equality with same field type"
            return self._value==number._value
        else : raise ValueError(f"type {type(number) } unsuported for checking equality with F_p")




    def __repr__(self):
        return f"{self._value} (mod {self.prime_modulus})"




A=F_p(3,7)
B=F_p(3,5)


a=5