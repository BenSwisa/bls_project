'''

NOTE: there is a difference between rmul and mul in terms of return mod 

        in the curve y^2=x^3+A*x+B

        
        the 'scalars' A,B could be diferent modulus then the underlying modules , 
        but the curve lives in the module of x,y therfore the ruturning modolus is of x,y 

        

        similar notic for aub and add
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
            return self.__class__((self._value+number._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for add with F_p")


    def __radd__(self,number):
        if isinstance(number,int):
            return self.__class__((self._value+number),self.prime_modulus)
        if isinstance(number,self.__class__):
            return self.__class__((self._value+number._value),number.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for add with F_p")
    

    def __sub__(self,number):
        if isinstance(number,int):
            return self.__class__((self._value-number),self.prime_modulus)
        if isinstance(number,self.__class__):
            return self.__class__((self._value-number._value),self.prime_modulus)
        else : raise ValueError(f"type {type(number) } unsuported for sub with F_p")


    def __rsub__(self,number):
        if isinstance(number,int):
            return self.__class__((number-self._value),self.prime_modulus)
        if isinstance(number,self.__class__):
            return self.__class__((number._value-self._value),number.prime_modulus)
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


    def sqrt(self): #TODO: can i get multiple sqrt values?
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


class PolyF_p:
    def __init__(self,values : dict):

        self.as_dict=values
        for scalar in values.values() : assert isinstance(scalar,F_p), "scalars must be defined over a field of type F_p"
        if len(values.items())>1:
            v0=list(values.values())[0]
            for power,scalar in values.items() : assert v0.prime_modulus==scalar.prime_modulus , "scalars must be over the same field"
    
    def copy(self):
        return self.__class__(self.as_dict)


    
    def __add__(self,poly):
        if isinstance(poly,self.__class__):
            d={}
            for power1,scalar1 in self.as_dict.items():
                flag=0
                for power2,scalar2 in poly.as_dict.items():
                    if power1==power2:
                        d[power2]=scalar1+scalar2 
                        flag=1
                if not flag: d[power1]=scalar1

            for power1,scalar1 in poly.as_dict.items():
                flag=0
                for power2,scalar2 in d.items():
                    if power1==power2:                        
                        flag=1
                if not flag: d[power1]=scalar1

            temp={}
            for power2,scalar2 in d.items(): 
                if scalar2._value!=0: temp[power2]=d[power2]
            return self.__class__(temp)
        else : raise ValueError(f"type {type(poly) } unsuported for add with polyF_p")


    

    def __sub__(self,poly):
        if isinstance(poly,self.__class__):
            d={}
            for power1,scalar1 in self.as_dict.items():
                flag=0
                for power2,scalar2 in poly.as_dict.items():
                    if power1==power2:
                        d[power2]=scalar1-scalar2
                        flag=1
                if not flag: d[power1]=scalar1

            for power1,scalar1 in poly.as_dict.items():
                flag=0
                for power2,scalar2 in d.items():
                    if power1==power2:                        
                        flag=1
                if not flag: d[power1]=-scalar1

            temp={}
            for power2,scalar2 in d.items(): 
                if scalar2._value!=0: temp[power2]=d[power2]
            return self.__class__(temp)
        else : raise ValueError(f"type {type(poly) } unsuported for add with polyF_p")



    def __mul__(self,poly):
        if isinstance(poly,self.__class__):
            d={}
            for power1,scalar1 in self.as_dict.items():
                for power2,scalar2 in poly.as_dict.items():
                    pow= power1+power2
                    if pow not in d : d[pow]=scalar1*scalar2
                    else : d[pow]=d[pow]+scalar1*scalar2

            temp={}
            for power2,scalar2 in d.items(): 
                if scalar2._value!=0: temp[power2]=d[power2]
            return self.__class__(temp)
        else : raise ValueError(f"type {type(poly) } unsuported for add with polyF_p")


                    
    def __divmod__(self, poly):
            """
            Performs Euclidean division. Returns (Quotient, Remainder).
            """
            if not isinstance(poly, self.__class__):
                raise ValueError(f"type {type(poly)} unsupported for divmod with PolyF_p")
            
            if not poly.as_dict:
                raise ZeroDivisionError("Polynomial division by zero")

            # Find the degree and leading coefficient of the divisor B(x)
            deg_divisor = max(poly.as_dict.keys())
            lead_divisor = poly.as_dict[deg_divisor]

            Q_dict = {}
            R = self # The remainder starts as the dividend A(x)

            while R.as_dict:
                # Find the degree and leading coefficient of the current Remainder
                deg_R = max(R.as_dict.keys())
                
                # If the remainder's degree is less than the divisor's, we are done
                if deg_R < deg_divisor:
                    break
                
                lead_R = R.as_dict[deg_R]
                
                # Calculate the degree and coefficient for the new term in the quotient
                term_deg = deg_R - deg_divisor
                
                # This triggers your F_p.__truediv__ to divide the coefficients mod p!
                term_coeff = lead_R / lead_divisor 
                
                # Add this new term to our quotient dictionary
                Q_dict[term_deg] = term_coeff
                
                # Create a single-term polynomial to multiply and subtract
                term_poly = self.__class__({term_deg: term_coeff})
                
                # R = R - (term * B)
                # This requires your __sub__ and __mul__ to be working perfectly
                R = R - (term_poly * poly)

            return self.__class__(Q_dict), R

    def __floordiv__(self, poly):
        """ Maps the // operator to return just the Quotient """
        return divmod(self, poly)[0]

    def __mod__(self, poly):
        """ Maps the % operator to return just the Remainder """
        return divmod(self, poly)[1]

    def __pow__(self,number):
        if isinstance(number,int):
            assert (abs(number)>=1 ) , "only positive powers are suported" 
            temp=self.copy()
            for i in range(number-1):
                temp=temp*self.copy()
            return temp
        else : raise ValueError(f"type {type(number) } unsuported for pow with F_p")

    def inv():
        raise Exception("A general polynomial does not have an inverese, iverse is supported over quetients only")

    # def mod_inverse(self, mod_poly):
    #         """
    #         Calculates the multiplicative inverse of this polynomial modulo mod_poly
    #         using the Extended Euclidean Algorithm.
    #         """
    #         if not isinstance(mod_poly, self.__class__):
    #             raise ValueError(f"type {type(mod_poly)} unsupported for mod_inverse")
                
    #         if not mod_poly.as_dict:
    #             raise ValueError("Cannot modulo by the zero polynomial.")
    #         if not self.as_dict:
    #             raise ValueError("The zero polynomial has no inverse.")
                
    #         # Get a sample scalar to extract the prime field modulus and class type
    #         sample_scalar = list(mod_poly.as_dict.values())[0]
    #         p = sample_scalar.prime_modulus
    #         F_p_class = sample_scalar.__class__
            
    #         # Define the '0' and '1' polynomials for our base cases
    #         zero = self.__class__({})
    #         one = self.__class__({0: F_p_class(1, p)})
            
    #         # Setup variables for the Extended Euclidean Algorithm
    #         # r0 and r1 track the remainders (GCD calculation)
    #         # t0 and t1 track the Bézout coefficients (the inverse)
    #         r0 = mod_poly
    #         r1 = self
    #         t0 = zero
    #         t1 = one
            
    #         # Loop while the remainder polynomial is not zero
    #         while r1.as_dict:
    #             # Get quotient and remainder
    #             q, r_temp = divmod(r0, r1)
                
    #             # Shift the remainders
    #             r0 = r1
    #             r1 = r_temp
                
    #             # Shift the Bézout coefficients: t_new = t0 - q * t1
    #             t_temp = t0 - (q * t1)
    #             t0 = t1
    #             t1 = t_temp
                
    #         # r0 now holds the Greatest Common Divisor (GCD)
    #         deg_gcd = max(r0.as_dict.keys()) if r0.as_dict else -1
            
    #         # If the degree of the GCD is > 0, they share a polynomial factor.
    #         # This means they are not coprime, so no inverse exists.
    #         if deg_gcd > 0 or deg_gcd == -1:
    #             raise ValueError("Polynomials are not coprime; no inverse exists.")
                
    #         # The GCD is a scalar (e.g., 5). We must multiply our result by 
    #         # the inverse of this scalar to normalize the equation to 1.
    #         gcd_scalar = r0.as_dict[0]
            
    #         # Calculate the modular inverse of the scalar (uses F_p.__pow__)
    #         gcd_scalar_inv = gcd_scalar ** -1 
            
    #         # Create a 0-degree polynomial from the scalar inverse
    #         inv_poly = self.__class__({0: gcd_scalar_inv})
            
    #         # Return the normalized coefficient
    #         return t0 * inv_poly

    def __repr__(self):

        values=self.as_dict
        if len(values.items())==0: return "0"
        s=""
        for ii,(power,scalar) in enumerate(list(values.items())) :
            if ii>0: s+="+ "
            if power==0: s+= f"{scalar._value} "
            else: s+=f"{scalar._value}*x^{power} "

        return s+f"(mod {scalar.prime_modulus})"


class PolyF_nQuotient(PolyF_p):
    def __init__(self, values, quotient_poly):
        super().__init__(values)

        #FIXME : defince all operations 


A=F_p(3,7)
B=F_p(3,5)
C=PolyF_p({0: A , 2: A})
D=PolyF_p({0: A , 3: A})

print(C**2)

a=5