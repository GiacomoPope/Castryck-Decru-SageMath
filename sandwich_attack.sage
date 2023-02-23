"""
This file demonstrates how to run an isogeny diamond attack
in the case when 2^(a-alpha) - k*3^b = u^2+v^2
where k is either a sum of squares, or a small smooth or prime integer.

This specific form is used by Damien Robert in the dimension 4 attack
(https://eprint.iacr.org/2022/1038.pdf)

This is notably possible in the SIKEp217 challenge and in the
(withdrawn) SIKEp964 parameter set.

There are 2 possible diagrams for this attack:
- the second auxiliary isogeny can be constructed from EB (always
  possible), in which case we sandwich the secret isogeny
  between the auxiliary isogenies.
  Preferably k is a smooth integer here

         phi       aux2 (deg k)
   Estart ---> EB ------> EB'
     |          |          |
     |          |          |
     |aux=u+iv  |          |
     |          |          |
     v          v          v
     C ------> CB ------> CB'
   
- or it can be realized as an endomorphism of E_start
  (preferably k is a sum of squares)

         aux2 (deg k)
   Estart ---> Estart ------> EB'
     |           |            |
     |           |            |
     |aux=u+iv   |            |
     |           |            |
     v           v            v
   Estart ---> Estart ------> CB'
"""

# Python imports
import time
from itertools import product

# Local Imports
from helpers import supersingular_gens, fast_log3
from richelot_aux import Does22ChainSplit, Pushing3Chain

# Load Sage Files
load('speedup.sage')

# TODO: implement the first strategy (when k is not a sum of squares).

def SandwichAttack(E_start, P2, Q2, EB, PB, QB, two_i, k, alp):
    "Implementation of the second strategy (endomorphism of degree k)"

    tim = time.time()

    # FIXME: a, b are magically in global scope
    # Might be precomputed if *really* required
    v, u = two_squares(2^(a-alp) - k*3^b)
    vk, uk = two_squares(k)

    # aux1 = u + i*v
    # aux2 = uk+ i*vk
    # aux2inv = (uk - i * vk) / k

    # Need to twist 2-torsion and 3-torsion
    kinv2 = pow(k, -1, 2^a)
    kinv3 = pow(k, -1, 3^b)

    def aux(P, kinv):
        # Beware, we are given 2i, so swap coordinates depending on
        # parity.
        if vk % 2 == 0:
            x = uk*P - (vk//2)*two_i(P)
        else:
            x = vk*P - (uk//2)*two_i(P)
        x *= kinv
        if v % 2 == 0:
            return u*x + (v//2)*two_i(x)
        else:
            return v*x + (u//2)*two_i(x)

    P_c = aux(P2, kinv2)
    Q_c = aux(Q2, kinv2)
    # FIXME: P3 and Q3 are magically in scope
    P3_c = aux(P3, kinv3)
    Q3_c = aux(Q3, kinv3)

    chain, (E1, E2) = Does22ChainSplit(E_start, EB,
        2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, a-alp)

    # Evaluate quotient map
    if E1.j_invariant() == E_start.j_invariant():
        index = 1
        CB = E2
    else:
        index = 0
        CB = E1
    def C_to_CB(x):
        pt = (x, None)
        for c in chain:
            pt = c(pt)
        return pt[index]

    P3_CB = C_to_CB(P3_c)
    Q3_CB = C_to_CB(Q3_c)

    print("Computed image of 3-adic torsion in split factor C_B")
    Z3 = Zmod(3^b)
    G1_CB, G2_CB = supersingular_gens(CB)
    G1_CB3 = ((p+1) / 3^b) * G1_CB
    G2_CB3 = ((p+1) / 3^b) * G2_CB
    w = G1_CB3.weil_pairing(G2_CB3, 3^b)

    sk = None
    for G in (G1_CB3, G2_CB3):
        xP = fast_log3(P3_CB.weil_pairing(G, 3^b), w)
        xQ = fast_log3(Q3_CB.weil_pairing(G, 3^b), w)
        if xQ % 3 != 0:
            sk = int(-Z3(xP) / Z3(xQ))
            break

    if sk is not None:
        # Sanity check
        bobscurve, _ = Pushing3Chain(E_start, P3 + sk*Q3, b)
        found = bobscurve.j_invariant() == EB.j_invariant()

        print(f"Bob's secret key revealed as: {sk}")
        print(f"In ternary, this is: {Integer(sk).digits(base=3)}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return sk
    else:
        print("Something went wrong.")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return None
