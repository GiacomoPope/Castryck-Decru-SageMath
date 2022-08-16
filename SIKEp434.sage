import public_values_aux
from public_values_aux import *

load('castryck_decru_shortcut.sage')
load('sandwich_attack.sage')

SIKE_parameters = {
    "SIKEp434" : (216, 137),
    "SIKEp503" : (250, 159),
    "SIKEp610" : (305, 192),
    "SIKEp751" : (372, 239),
    "SIKEp964" : (486, 301), # removed after NIST round 1
}

# Change me to attack different parameter sets
NIST_submission = "SIKEp434"
a, b = SIKE_parameters[NIST_submission]

print(f"Running the attack against {NIST_submission} parameters, which has a prime: 2^{a}*3^{b} - 1")

print(f"Generating public data for the attack...")
# Set the prime, finite fields and starting curve
# with known endomorphism
p = 2^a*3^b - 1
public_values_aux.p = p
Fp2.<i> = GF(p^2, modulus=x^2+1)
R.<x> = PolynomialRing(Fp2)

E_start = EllipticCurve(Fp2, [0,6,0,1,0])
E_start.set_order((p+1)^2, num_checks=0) # Speeds things up in Sage

# Generation of the endomorphism 2i
two_i = generate_distortion_map(E_start)

# Generate public torsion points, for SIKE implementations
# these are fixed but to save loading in constants we can
# just generate them on the fly
P2, Q2, P3, Q3 = generate_torsion_points(E_start, a, b)
check_torsion_points(E_start, a, b, P2, Q2, P3, Q3)

# Generate Bob's key pair
bob_private_key, EB, PB, QB = gen_bob_keypair(E_start, b, P2, Q2, P3, Q3)
solution = Integer(bob_private_key).digits(base=3)

print(f"If all goes well then the following digits should be found: {solution}")

# ===================================
# =====  ATTACK  ====================
# ===================================

def RunAttack(num_cores):
    return CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=num_cores)

if __name__ == '__main__' and '__file__' in globals():
    if '--parallel' in sys.argv:
        # Set number of cores for parallel computation
        num_cores = os.cpu_count()
        print(f"Performing the attack in parallel using {num_cores} cores")
    else:
        num_cores = 1

    if '--sandwich' in sys.argv:
        # Use the fact that 2^a - 5*3^b is a sum of squares
        assert NIST_submission == "SIKEp964"
        assert two_squares(2^a - 5*3^b)
        recovered_key = SandwichAttack(E_start, P2, Q2, EB, PB, QB, two_i, k=5, alp=0)
    else:
        recovered_key = RunAttack(num_cores)

