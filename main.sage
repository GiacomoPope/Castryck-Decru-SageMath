import sys
import argparse

pXXX = ['p182', 'p217', 'p377', 'p546', 'p697', 'p434', 'p503', 'p610', 'p751']
# -------------------------------
def arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("--prime", type=str, help="prime field characteristic,", required=True, choices=pXXX)
    parser.add_argument("--strategies", help="optimization using strategies,", action='store_true')
    parser.add_argument("--shortcut", help="shortcut optimization,", action='store_true')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    options = parser.parse_args(args)
    return options

import public_values_aux
from public_values_aux import *

# SIKEpXXX parameters
pxxx = arguments(sys.argv[1:]).prime
strategies = arguments(sys.argv[1:]).strategies
shortcut = arguments(sys.argv[1:]).shortcut

if pxxx == 'p182':
    a, b = 91, 57
elif pxxx == 'p217':
    a, b = 110, 67
elif pxxx == 'p377':
    a, b = 191, 117
elif pxxx == 'p546':
    a, b = 273, 172
elif pxxx == 'p697':
    a, b = 356, 215
elif pxxx == 'p434':
    a, b = 216, 137
elif pxxx == 'p503':
    a, b = 250, 159
elif pxxx == 'p610':
    a, b = 305, 192
elif pxxx == 'p751':
    a, b = 372, 239
else:
    exit(-1)

if not shortcut:
    load('castryck_decru_attack.sage')
else:
    load('castryck_decru_shortcut.sage')


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

print(f"Running the attack against SIDHp{p.bit_length()} parameters, which has a prime: 2^{a}*3^{b} - 1")
print(f"If all goes well then the following digits should be found: {solution}")


def RunAttack(num_cores):
        return CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=num_cores, strategies=strategies)


if __name__ == '__main__' and '__file__' in globals():
    if '--parallel' in sys.argv:
        # Set number of cores for parallel computation
        num_cores = os.cpu_count()
        print(f"Performing the attack in parallel using {num_cores} cores")
    else:
        num_cores = 1
    recovered_key = RunAttack(num_cores)

