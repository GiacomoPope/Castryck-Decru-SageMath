"""
Demonstrate how to run the attack when c = small_prime*(u^2+v^2)
This is a mixture of Castryck-Decru and Maino-Martindale

The following variants are implemented:
Case 7: 2^301-3^188 == 7*(u*u+v*v) (3^4 guesses)

        tau
    Estart --> Eguess ---> EB
    | u+iv     |          |
    v          |          |
    Estart     |aux       |
    | phi7     |          |
    v          v          v
    E7 ------> C -------> CB
        tau7

(note that E(p^3) has (p^3+1) points, and a 7-torsion point)

Case 11: 2^305-19*3^189 = 11*(u*u+v*v) (3^3 guesses)

        tau
    Estart --> Eguess ---> EB ---> EB19
    | u+iv       |           phi19  |
    v            |                  |
    Estart       |aux               |
    | phi11      |                  |
    v            v                  v
    E11 -------> C ---------------> CB
        tau11

(note that E(p^4) has (p^2-1)^2 points, and a 11-torsion point)
(note that E(p^9) has p^9+1 points, and a 19-torsion point)

"""

import time
import argparse

# Local imports
from helpers import possibly_parallel, supersingular_gens, fast_log3
import public_values_aux
from public_values_aux import *
from richelot_aux import Pushing3Chain, Does22ChainSplit

# Load Sage Files
load('speedup.sage')

set_verbose(-1)

argp = argparse.ArgumentParser()
argp.add_argument("--parallel", action="store_true")
argp.add_argument("--mode", type=int, choices=(7, 11), default=11,
    help="mode 11 needs 27 guesses, mode 7 needs 81 guesses")
opts = argp.parse_args()

print("Instantiate E_start...")
a, b = 305, 192
p = 2^a*3^b - 1
public_values_aux.p = p

Fp2.<i> = GF(p^2, modulus=x^2+1)
R.<x> = Fp2[]

E_start = EllipticCurve(Fp2, [0,6,0,1,0])
E_start.set_order((p+1)^2, num_checks=0) # Speeds things up in Sage

# Generation of the endomorphism 2i
two_i = generate_distortion_map(E_start)

# Choose an outgoing degree 7 isogeny.
# We don't need to compute a torsion point in GF(p^6),
# we can factor the division polynomial and use Kohel formulas.
if opts.mode == 7:
    print("Using 2^301-3^188 = 7*(u²+v²) (max guesses = 81)")
    print("Precompute an isogeny of degree 7 on E_start")
    P7 = E_start.division_polynomial(7)
    ker7, _ = P7.factor()[0]
    phi_left = E_start.isogeny(ker7)
    print(phi_left)

    u = 714020003029005719823753224880815399155339403
    v = 28031663375683401880549715251102056676622848
    assert 2^301-3^188 == 7*(u*u+v*v)
else:
    print("Using 2^305-19*3^189 = 11*(u²+v²) (max guesses = 27)")
    print("Precompute an isogeny of degree 11 on E_start")
    # Sage is confused by degree 1 factors
    #P11 = E_start.division_polynomial(11)
    #ker11, _ = P11.factor()[0]
    phi_left = E_start.isogenies_prime_degree(11)[0]
    print(phi_left)

    u = 1550193735342211609960431880234414865472178337
    v = 965894233082293540389053428058397267855495894
    assert 2^305-19*3^189 == 11*(u*u+v*v)

# Generate public torsion points, for SIKE implementations
# these are fixed but to save loading in constants we can
# just generate them on the fly
P2, Q2, P3, Q3 = generate_torsion_points(E_start, a, b)
check_torsion_points(E_start, a, b, P2, Q2, P3, Q3)

# Generate Bob's key pair
bob_private_key, EB, PB, QB = gen_bob_keypair(E_start, b, P2, Q2, P3, Q3)
solution = Integer(bob_private_key).digits(3, padto=b)

print(f"If all goes well then the following digits should be found: {solution}")

# Build the following diagram for each guess
#
#       tau
# Estart --> Eguess ---> EB
#   | u+iv     |          |
#   v          |          |
# Estart       |aux       |
#   | phi7     |          |
#   v          v          v
#  E7 -------> C ------> CB
#      tau7
#
# The isogeny Eguess->EB is secret, isogeny diamond is (Eguess, EB, C, CB)
# We don't compute aux, nor CB

if opts.parallel:
    # Set number of cores for parallel computation
    num_cores = os.cpu_count()
    print(f"Performing the attack in parallel using {num_cores} cores")
else:
    num_cores = 1

# Attack starts here

tim = time.time()

if opts.mode == 7:
    phiB = EB.identity_morphism()
    # SAGE identity morphisms are not functions??
    phiB._call_ = lambda x: x
    beta = 4 # 192-188
    alp = 4 # 305-301
else:
    print("Compute (once) an isogeny of degree 19 on E_B")
    P19 = EB.division_polynomial(19)
    ker19, _ = P19.factor()[0]
    phiB = EB.isogeny(ker19)
    print(phiB)
    beta = 3 # 192-189
    alp = 0 # 305-305
    print(f"... done in {time.time()-tim:.3f} seconds")

@possibly_parallel(num_cores)
def CheckGuess(guess):
    guess = Integer(guess)
    first_digits = guess.digits(3, padto=beta)
    print(f"Testing digits: {first_digits}")

    guessker = 3^(b-beta) * (P3 + guess*Q3)
    guessker_left = phi_left(u * guessker + (v//2) * two_i(guessker))
    C, tau_C = Pushing3Chain(phi_left.codomain(), guessker_left, beta)

    # Now compute the image of 2-torsion in C
    def Estart_to_C(x):
        x = u * x + (v//2) * two_i(x) # endomorphism
        x = phi_left(x)
        for c in tau_C:
            x = c(x)
        return x

    P2_C = Estart_to_C(P2)
    Q2_C = Estart_to_C(Q2)

    # Replace EB by the codomain of phiB
    split = Does22ChainSplit(C, phiB.codomain(),
        2^alp*P2_C, 2^alp*Q2_C,
        2^alp*phiB(PB), 2^alp*phiB(QB), a-alp)

    if split:
        Eguess, _ = Pushing3Chain(E_start, guessker, beta)
        chain, (E1, E2) = split
        if E1.j_invariant() == Eguess.j_invariant():
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
        P3_C = Estart_to_C(P3)
        Q3_C = Estart_to_C(Q3)
        P3_CB = C_to_CB(P3_C)
        Q3_CB = C_to_CB(Q3_C)

        print("Computed image of 3-adic torsion in split factor C_B")
        Z3 = Zmod(3^b)
        G1_CB, G2_CB = supersingular_gens(CB)
        G1_CB3 = ((p+1) / 3^b) * G1_CB
        G2_CB3 = ((p+1) / 3^b) * G2_CB
        w = G1_CB3.weil_pairing(G2_CB3, 3^b)

        xP = fast_log3(P3_CB.weil_pairing(G1_CB3, 3^b), w)
        xQ = fast_log3(Q3_CB.weil_pairing(G1_CB3, 3^b), w)
        if xQ % 3 != 0:
            sk = int(-Z3(xP) / Z3(xQ))
            return sk
        xP = fast_log3(P3_CB.weil_pairing(G2_CB3, 3^b), w)
        xQ = fast_log3(Q3_CB.weil_pairing(G2_CB3, 3^b), w)
        if xQ % 3 != 0:
            sk = int(-Z3(xP) / Z3(xQ))
            return sk

        raise Exception("fail?!")

    return None

for result in CheckGuess(list(range(3^beta))):
    ((guess,), _), sk = result
    if sk is not None:
        print("Glue-and-split! These are most likely the secret digits.")
        bobskey = sk
        break

# Sanity check
bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
found = bobscurve.j_invariant() == EB.j_invariant()

if found:
    print(f"Bob's secret key revealed as: {bobskey}")
    print(f"In ternary, this is: {Integer(bobskey).digits(base=3, padto=b)}")
    print(f"Altogether this took {time.time() - tim:.3f} seconds.")
else:
    print("Something went wrong.")
    print(f"Altogether this took {time.time() - tim:.3f} seconds.")
