from Verifiers import IsTheElipticCurveSmooth
from scripts.Rings import F_p
from GeneralFunctions import *


def int_to_bits(n):
    return [int(b) for b in bin(n)[2:]]





'''
    the general curve is defined using the 'scalars' A,B :

    E={(x,y): y^2=x^3+Ax_b} union {0} s.t A,B in F_p

    while E(F_p)={(x,y)     in F_p^2 : y^2=x^3+Ax_b} union {0}

    **A,B could belong to different fields than the points**

'''
class ElipticCurve:
    def __init__(self,A:F_p,B:F_p):
        self.A=A
        self.B=B
        assert self.A.prime_modulus==self.B.prime_modulus, " A,B not from the same field"
        
        
        assert(IsTheElipticCurveSmooth(A,B)),"Curve is not Smooth"




    def __repr__(self):
        return f"y^2=x^3+[{self.A}]*x+[{self.B}]"
    

    def is_on_curve(self, P) -> bool:
        if P.is_infinity():
            return True
        x, y = P.x, P.y
        return ((y * y) == (x **3 + self.A * x + self.B))

    def inf(self):
        return ElipticCurvePoint( None, None,self)

    def getPoint(self,x):
        return ElipticCurvePoint(x,(x**3+self.A*x+self.B).sqrt(),self)
    

    def __eq__(self, another_curve):
        return self.A==another_curve.A and self.B==another_curve.B 


class RationalPointsGroupOfElipticCurve:
    def __init__(self,curve :ElipticCurve , prime_modulus):
        self.curve=curve
        self.prime_modulus=prime_modulus

        self.group_order=self._calculate_group_order_naive()

    def __repr__(self):
        return f"E(F_{self.prime_modulus})"
    

    #TODO: this implementation is inefficient
    def _calculate_group_order_naive(self):
        """
        Calculates the order N of the elliptic curve E(F_p): y^2 = x^3 + A*x + B
        using naive counting and Euler's criterion.
        
        A, B: Curve parameters (integers modulo p)
        p: The prime field characteristic
        """
        # Start with 1 to account for the point at infinity (O)
        N = 1
        p=self.prime_modulus
        A=self.curve.A
        B=self.curve.B

        # Exponent for Euler's criterion: (p-1)/2
        euler_exp = (p - 1) // 2

        for x in range(p):
            # Evaluate the right-hand side of the curve equation: z = x^3 + Ax + B mod p
            x=F_p(x,p)
            z = x**3 + A * x + B
            
            if z == 0:
                # If z is 0, y must be 0. There is exactly one point: (x, 0)
                N += 1
            else:
                # Use Euler's criterion to check if z is a quadratic residue (a perfect square)
                legendre = z**euler_exp
                
                if legendre == 1:
                    # z is a quadratic residue, meaning there are two distinct y values (y and -y)
                    N += 2
                # If legendre == p - 1 (which is -1 mod p), z is a non-residue. 
                # There are 0 points for this x, so we add nothing.

        return N

# Example usage with a small prime to verify:
# Curve: y^2 = x^3 + 2x + 2 over F_17
# p = 17, A = 2, B = 2
# Expected order N = 19
# print(calculate_group_order_naive(2, 2, 17))


    # def is_on_curve(self, P) -> bool:
    #     if P.is_infinity():
    #         return True
    #     x, y = P.x, P.y
    #     return ((y * y) == (x **3 + self.A * x + self.B))

    # def inf(self):
    #     return ElipticCurvePoint( None, None,self)

    # def getPoint(self,x):
    #     return ElipticCurvePoint(x,(x**3+self.A*x+self.B).sqrt(),self)
    

    # def __eq__(self, another_curve):
    #     return self.A==another_curve.A and self.B==another_curve.B 



class ElipticCurvePoint:
    def __init__(self,x: F_p,y: F_p,curve : ElipticCurve):
        self.x=x
        self.y=y
        self.curve=curve
        # self.p=curve.prime_modulus #TODO
        if x is not None and y is not None:
            assert isinstance(x,F_p) and isinstance(y,F_p) , "Points must be defined over a field"
            assert x.prime_modulus==y.prime_modulus , "x,y must be defined over the same field"
            # assert(curve.is_on_curve(self)),("Point is not on the Curve") #TODO

        self._order=None

    @property
    def order(self):
        if self._order is None:
            self._order=self._get_point_order()
        return self._order

    def is_infinity(self) -> bool:
        return self.x is None or self.y is None

    def copy(self):
        return self.__class__(self.x,self.y,self.curve)

    def inv(self):
        if self.is_infinity() : return self
        return self.__class__(self.x,-self.y,self.curve)


    def __add__(self,Q):
        assert isinstance(Q,self.__class__) , "Addition is suported over two points only"
        if Q.is_infinity():
            return self
        if self.is_infinity():
            return Q
        
        assert Q.x.prime_modulus==self.x.prime_modulus , "P,Q must be defined over the same field"


        x1, y1 = self.x, self.y
        x2, y2 = Q.x, Q.y

        #case 1
        if x1 != x2 :
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
        assert isinstance(k,int), "multiplication is supported over an integer only"
        if k==0: return self.__class__(None,None,self.curve)
        if k < 0:
            if k==-1 : return self.inv()
            k*=-1
            summand = self.inv().copy()
        else : summand = self.copy()

        for bit in int_to_bits(k)[1:]:
            summand=summand+summand
            if bit: summand= summand+self.copy()

        return summand


    def __rmul__(self, k: int) : #TODO: make sure to check where multiplication is supported from left and from right
        return self.__mul__(k)

    def __eq__(self, Q):
        return (self.x==Q.x) and (self.y==Q.y) and self.curve==Q.curve



    def _get_point_order(self):
        if self.is_infinity(): return 1
        group=RationalPointsGroupOfElipticCurve(self.curve,self.x.prime_modulus)
        prime_fact=prime_factors(group.group_order)

        for fact in prime_fact:
            temp=fact*self
            if temp.is_infinity() : return fact 

        if (group.group_order*self.copy()).is_infinity(): return group.group_order

        return -1

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
        return f"({self.x._value}, {self.y._value}) mod {self.x.prime_modulus}"
        
mod=5
curve=ElipticCurve(F_p(1,mod),F_p(1,mod))
P=ElipticCurvePoint(F_p(0,mod),F_p(1,mod),curve)
Q=ElipticCurvePoint(F_p(2,mod),F_p(1,mod),curve)

# Q=a.getPoint(1)[0]
print(P+Q)



E541=RationalPointsGroupOfElipticCurve(curve,541)


b=5
