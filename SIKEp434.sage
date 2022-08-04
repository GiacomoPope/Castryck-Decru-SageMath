import time
from itertools import product

load('richelot_aux.sage')
load('uvtable.sage')
load('speedup.sage')

# Remove annoying messages about slow Gröbner
set_verbose(-1)

SIKE_parameters = {
    "SIKEp434" : (216, 137),
    "SIKEp503" : (250, 159),
    "SIKEp610" : (305, 192),
    "SIKEp751" : (372, 239)
}

# Change me to attack different parameter sets
NIST_submission = "SIKEp434"
a, b = SIKE_parameters[NIST_submission]

print(f"Running the attack against {NIST_submission} parameters, which has a prime: 2^{a}*3^{b} - 1")

p = 2^a*3^b - 1
Fp2.<i> = GF(p^2, modulus=x^2+1)
assert i^2 == -1

R.<x> = PolynomialRing(Fp2)

E_start = EllipticCurve(Fp2, [0,6,0,1,0])
# Speeds things up in Sage
E_start.set_order((p+1)^2)

phi = EllipticCurveIsogeny(E_start, x)
E1728 = phi.codomain()
# Speeds things up in Sage
E1728.set_order((p+1)^2)

for iota in E1728.automorphisms():
    P = E1728.random_point()
    if iota(iota(P)) == -P:
        two_i = -phi.dual()*iota*phi
        break

infty = E_start(0)

def get_l_torsion_basis(E, l):
    n = (p+1) // l
    return (n*G for G in E.gens())

P2, Q2 = get_l_torsion_basis(E_start, 2^a)
P3, Q3 = get_l_torsion_basis(E_start, 3^b)

# Make sure Torsion points are
# generated correctly
assert 2^(a-1)*P2 != infty
assert 3^(b-1)*P3 != infty
assert P2.weil_pairing(Q2, 2^a)^(2^(a-1)) != 1
assert P3.weil_pairing(Q3, 3^b)^(3^(b-1)) != 1

# generate challenge key
Bobskey = randint(0,3^b)

EB, chain = Pushing3Chain(E_start, P3 + Bobskey*Q3, b)
# Speeds things up in Sage
EB.set_order((p+1)^2)
PB = P2
for c in chain:
    PB = c(PB)
QB = Q2 
for c in chain:
    QB = c(QB)

print(f"If all goes well then the following digits should be found: {Integer(Bobskey).digits(base=3)}")

# ===================================
# =====  ATTACK  ====================
# ===================================
tim = time.time()

skB = [] # TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY

# gathering the alpha_i, u, v from table

expdata = [[0, 0, 0] for _ in range(b-3)]
# for i in [1..b-3] do
for i in range(1,b-2):
    # if IsOdd(b-i) then
    if (b-i)%2 == 1:
        index = (b - i + 1) // 2
        exp = uvtable[index-1][1]
        if exp <= a:
            u = uvtable[index-1][2]
            v = uvtable[index-1][3]
            expdata[i-1] = [exp, u, v]

# gather digits until beta_1
bet1 = 0
while expdata[bet1][0] == 0:
    bet1 += 1
bet1 += 1

ai = expdata[bet1-1][0]
u  = expdata[bet1-1][1]
v  = expdata[bet1-1][2]

print(f"Determination of first {bet1} ternary digits. We are working with 2^{ai}-torsion.")

bi = b - bet1
alp = a - ai

# for j in CartesianPower([0,1,2], bet1) do
for j in product([0,1,2], repeat=int(bet1)):
    print(f"Testing digits: {[j[k] for k in range(bet1)]}")

    # tauhatkernel = 3^bi*P3 + sum([3^(k-1)*j[k-1] for k in range(1,beta+1)])*3^bi*Q3
    tauhatkernel = 3^bi*P3 
    for k in range(1, bet1+1):
        tauhatkernel += (3^(k-1)*j[k-1])*3^bi*Q3

    tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)

    C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, bet1)

    P_c = u*P2 + v*two_i(P2) 
    for taut in tau_tilde:
        P_c = taut(P_c)
    Q_c = u*Q2 + v*two_i(Q2)
    for taut in tau_tilde:
        Q_c = taut(Q_c)

    # if j eq <2 : k in [1..bet1]> or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai) then
    if j == (2,)*bet1 or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai):
        print("Glue-and-split! These are most likely the secret digits.")
        for k in j:
            skB.append(k)
        break

print(skB)


# now compute longest prolongation of Bob's isogeny that may be needed
length = 1
max_length = 0
# for i in [bet1+1..b-3] do
for i in range(bet1 + 1, b - 2):
  if expdata[i-1][0] != 0:
    max_length = max(length, max_length)
    length = 0
  else:
    length += 1

while True:
    K = 2^a*3^(b - max_length)*EB.random_point()
    if K.order() == 3^max_length:
        break

while True:
    alternativeK = 2^a*3^(b - max_length)*EB.random_point()
    if K.weil_pairing(alternativeK, 3^max_length)^(3^(max_length - 1)) != 1:
        break

_, EBprolong = Pushing3Chain(EB, K, max_length)

# gather next digit and change prolongation if needed
i = bet1 + 1
bi = b - i

print(f"Determination of the {i}th ternary digit. We are working with 2^{ai}-torsion.")
prolong = 1
print("Prolonging with 1 steps.")
endPB = EBprolong[0](PB)
endQB = EBprolong[0](QB)
endEB = EBprolong[0].codomain()
# Speeds things up in Sage
endEB.set_order((p+1)^2)

positives = []

for j in range(0,2+1):
    print(f"Testing digit: {j}")
    # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
    tauhatkernel = 3^bi*P3 
    for k in range(1, i):
        tauhatkernel += (3^(k-1)*skB[k-1])*3^bi*Q3
    tauhatkernel += 3^(i-1)*j*3^bi*Q3

    tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)

    C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, i)
    P_c = u*P2 + v*two_i(P2)
    for taut in tau_tilde:
        P_c = taut(P_c)
    Q_c = u*Q2 + v*two_i(Q2)
    for taut in tau_tilde:
        Q_c = taut(Q_c)

    if Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai):
        print("Glue-and-split!")
        positives.append(j)
        break

if len(positives) == 1:
    print(f"Most likely good prolongation and the secret digit is: {positives[0]}")
    skB.append(positives[0])
    next_i = i + 1
else:
    print("All glue-and-splits, must be bad prolongation: changing it and redoing this digit.")
    _, EBprolong = Pushing3Chain(EB, alternativeK, max_length)
    next_i = i
    prolong = 0


# now gather all remaining digits, except for last three (we close that gap by trial and error)

# for i in [next_i..b-3] do
for i in range(next_i, b-2):
    bi = b - i
    if expdata[i-1][0] != 0:
        ai = expdata[i-1][0]
        alp = a - ai
        u = expdata[i-1][1]
        v = expdata[i-1][2]
        prolong = 0
    else:
        prolong += 1

    print(f"Determination of the {i}th ternary digit. We are working with 2^{ai}-torsion.")
    print(f"Prolonging with {prolong} steps.")
    endPB = PB
    endQB = QB
    endEB = EB
    # for j in [1..prolong] do
    for j in range(1,prolong+1):
        endPB = EBprolong[j-1](endPB)
        endQB = EBprolong[j-1](endQB)
    if j == prolong:
        endEB = EBprolong[j-1].codomain()
        # Speeds things up in Sage
        endEB.set_order((p+1)^2)


    for j in range(0,3):
        print(f"Testing digit: {j}")
        # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
        tauhatkernel = 3^bi*P3 
        for k in range(1, i):
            tauhatkernel += (3^(k-1)*skB[k-1])*3^bi*Q3
        tauhatkernel += 3^(i-1)*j*3^bi*Q3

        tauhatkernel_distort =u*tauhatkernel + v*two_i(tauhatkernel)
        
        C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, i)
        P_c = u*P2 + v*two_i(P2)
        for taut in tau_tilde:
            P_c = taut(P_c)
        Q_c = u*Q2 + v*two_i(Q2)
        for taut in tau_tilde:
            Q_c = taut(Q_c)
        if j == 2 or Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai):
            print("Glue-and-split! This is most likely the secret digit.")
            skB.append(j)
            break

key = sum([skB[i-1]*3^(i-1) for i in range(1,b-2)])

# bridge last safety gap
tim2 = time.time()

found = false
for i in range(3^5+1):
    bobskey = key + i*3^(b-3)
    bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
    if bobscurve.j_invariant() == EB.j_invariant():
        found = true
        break

print(f"Bridging last gap took: {time.time() - tim2}")

if found:
    print(f"Bob's secret key revealed as: {bobskey}")
    print(f"In ternary, this is: {Integer(bobskey).digits(base=3)}")
else:
    print("Something went wrong.")

print(f"Altogether this took {time.time() - tim} seconds.")







