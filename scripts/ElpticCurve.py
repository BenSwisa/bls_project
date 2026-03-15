from Verifiers import IsTheElipticCurveSmooth
from Field import F_p

def int_to_bits(n):
    return [int(b) for b in bin(n)[2:]]

class ElipticCurve:
    def __init__(self,A:F_p,B:F_p):
        self.A=A
        self.B=B
        assert self.A.prime_modulus==self.B.prime_modulus, " A,B not from the same field"
        
        self.prime_modulus=self.A.prime_modulus
        
        assert(IsTheElipticCurveSmooth(A,B)),"Curve is not Smooth"


        self.order = self._compute_group_order_naive(A,B,self.prime_modulus) if self.prime_modulus!=1 else None 


    def __repr__(self):
        return f"y^2=x^3+{self.A}*x+{self.B}"
    

    def is_on_curve(self, P, epsilon=10e-10) -> bool:
        if P.is_infinity():
            return True
        x, y = P.x, P.y
        return ((y * y) - (x **3 + self.A * x + self.B))**2<epsilon

    def inf(self):
        return ElipticCurvePoint( None, None,self)

    def getPoint(self,x):
        return (ElipticCurvePoint(x,(x**3+self.A*x+self.B)**0.5,self),
                ElipticCurvePoint(x,-(x**3+self.A*x+self.B)**0.5,self))
    

    def __eq__(self, another_curve):
        return self.A==another_curve.A and self.B==another_curve.B and self.prime_modulus==another_curve.prime_modulus


    @staticmethod
    def _compute_group_order_naive(A, B, p): #TODO: more effiicient way?
        """
        Computes the order |G| of the elliptic curve y^2 = x^3 + Ax + B (mod p)
        using the naive counting algorithm.
        
        Note: This is only efficient for relatively small primes (e.g., p < 10^5).
        """
        # Start with 1 to account for the Point at Infinity (O)
        order = 1 
        
        # Iterate x from 0 to p - 1
        for x in range(p):
            # Calculate z = x^3 + Ax + B (mod p)
            z = (pow(x, 3, p) + A * x + B) % p
            
            if z == 0:
                # If z == 0, then y^2 = 0, which means y = 0. 
                # There is exactly one point: (x, 0)
                order += 1
            else:
                # Calculate the Legendre symbol using Euler's criterion: z^((p-1)/2) mod p
                legendre = pow(z, (p - 1) // 2, p)
                
                if legendre == 1:
                    # z is a quadratic residue. There are two solutions for y: (x, y) and (x, -y)
                    order += 2
                # If legendre == p - 1 (which acts as -1 mod p), there are no solutions for this x.
                
        return order


class ElipticCurvePoint:
    def __init__(self,x,y,curve : ElipticCurve):
        self.x=x
        self.y=y
        self.curve=curve
        self.p=curve.prime_modulus
        if x is not None and y is not None:
            assert(curve.is_on_curve(self)),("Point is not on the Curve")

    def is_infinity(self) -> bool:
        return self.x is None and self.y is None

    def copy(self):
        return self.__class__(self.x,self.y,self.curve)

    def inv(self):
        if self.is_infinity() : return self
        return self.__class__(self.x,-self.y % self.p,self.curve)


    def __add__(self,Q):
        if Q.is_infinity():
            return self
        if self.is_infinity():
            return Q
        
        x1, y1 = self.x, self.y
        x2, y2 = Q.x, Q.y

        #case 1
        if x1 != x2 and y1 != y2:
            lam = (y2 - y1) / (x2 - x1)

        #case 2
        if x1 == x2 and y1 == y2:
            lam = (x1 * x1 * 3 + self.curve.A) / (y1 * 2)
        
        #case 3
        if x1 == x2 and y1 == -y2:
            return self.curve.inf()

        x3 = lam * lam - x1 - x2
        y3 = lam * (x1 - x3) - y1
        return self.__class__(x3,y3,self.curve)




    def __mul__(self, k: int) :
        if k < 0:
            if k==-1 : return self.inv()
            k*=-1
            summand = self.inv().copy()
        else : summand = self.copy()

        for bit in int_to_bits(k)[1:]:
            summand=summand+summand
            if bit: summand= summand+self.copy()

        return summand


    def __rmul__(self, k: int) :
        return self.__mul__(k)

    def __eq__(self, Q,epsilon=10e10):
        return (self.x-Q.x)<epsilon and (self.y-Q.y)<epsilon and self.curve==Q.curve



    # def __sub__(self,number):
    #     return self.__class__((self._value-number)%self.prime_modulos)
        
    # def __mul__(self,number):
    #     return self.__class__((self._value*number)%self.prime_modulos)

    # @staticmethod
    # def _mod_inverse(a, m):
    #     r0, r1 = m, a
    #     s0, s1 = 1, 0
    #     t0, t1 = 0, 1

    #     while r1 != 0:
    #         q = r0 // r1
    #         r0, r1 = r1, r0 - q*r1
    #         s0, s1 = s1, s0 - q*s1
    #         t0, t1 = t1, t0 - q*t1

    #     return t0 % m
    
    # def __truediv__(self,number):
    #     return self.__class__((self._value*self._mod_inverse(number,self.prime_modulos))%self.prime_modulos)

    def __repr__(self) -> str:
        if self.is_infinity():
            return "O"
        return f"({self.x}, {self.y})"
        

a=ElipticCurve(F_p(5,1),F_p(5,1))
P=a.getPoint(5)[0]
Q=a.getPoint(1)[0]
print(P+P)
b=5
