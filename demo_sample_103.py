"""
Demo with the handout's sample input:
  p=103, A=1, B=0, a=7, m="שלום"
Exact Q / pairing value may differ (we choose random f(x) and random Q),
but verification should pass.
"""
from bls.find_generators import setup_system, sign, verify

def run():
    params = setup_system(103, 1, 0)
    Hm, sigma, pk = sign(params, 7, "שלום")
    ok, left, right = verify(params, Hm, sigma, pk)
    print("ok =", ok)
    print("left =", left)
    print("right =", right)

if __name__ == "__main__":
    run()
