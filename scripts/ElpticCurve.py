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


    #     self.order = self._compute_group_order_naive(A,B,self.prime_modulus) if self.prime_modulus!=1 else None 


    def __repr__(self):
        return f"y^2=x^3+[{self.A}]*x+[{self.B}]"
    

    def is_on_curve(self, P) -> bool:
        if P.is_infinity():
            return True
        x, y = P.x, P.y
        return ((y * y) == (x **3 + self.A * x + self.B))

    def inf(self):
        return ElipticCurvePoint( None, None,self)

    # def getPoint(self,x):
    #     return (ElipticCurvePoint(x,(x**3+self.A*x+self.B)**0.5,self),
    #             ElipticCurvePoint(x,-(x**3+self.A*x+self.B)**0.5,self))
    

    def __eq__(self, another_curve):
        return self.A==another_curve.A and self.B==another_curve.B 



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
b=5
