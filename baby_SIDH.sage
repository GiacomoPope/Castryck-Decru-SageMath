# Local imports
import public_values_aux
from public_values_aux import *

# Load Sage files
load('castryck_decru_shortcut.sage')

# Baby SIKEp64 parameters
a = 33
b = 19

# Set the prime, finite fields and starting curve
# with known endomorphism
p = 2^a*3^b - 1
public_values_aux.p = p

Fp2.<i> = GF(p^2, modulus=x^2+1)
R.<x> = PolynomialRing(Fp2)

E_start = EllipticCurve(Fp2, [0,6,0,1,0])
E_start.set_order((p+1)^2) # Speeds things up in Sage

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

print(f"Running the attack against Baby SIDHp64 parameters, which has a prime: 2^{a}*3^{b} - 1")
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
    recovered_key = RunAttack(num_cores)

