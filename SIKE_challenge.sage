import time
from itertools import product

load('richelot_aux.sage')
load('uvtable.sage')
load('speedup.sage')

# Remove annoying messages about slow Gröbner
set_verbose(-1)

# Stop slow primality checks GF(p^k) construction
proof.arithmetic(False)

# $IKEp217 parameters
a = 110
b = 67
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
        two_i = phi.post_compose(iota).post_compose(phi.dual())
        break

# $IKEp217 public parameters

xQ2 = 86445170414599058485968715662346922796016566052081315743781191908 + i*40452781127453632342062231901634407992350274206219311799642330230
yQ2 = 65598021239445364766945421641383356785387864046503211101901536990 + i*4185188043442820950612067662642469583838730096408235972690714750
Q2 = E_start(xQ2, yQ2)

xP2 = 37927755131506746411537091732732279589710442710570639382902295344 + i*37017495259074475516858568039563512752905548684528955375066616976
yP2 = 47930070613074499206421889519713097155043539387599009924307417901 + i*7176302698603683442486948733484190513242470057553561741081582098
P2 = E_start(xP2, yP2)

# CONSISTENCY WITH R (NOT NEEDED BUT OK)
xR2 = 46724562430016373759602265681716628459894400021233033976074511127 + i*62869247191516167619032311426583862201277153531817575978673280654
assert (P2 - Q2)[0] == xR2

xQ3 = 47793138785202808493310870472269548499693686682045542672866243781
yQ3 = i*59233392251697228025115695338640088714243311845504632846426147327
Q3 = E_start(xQ3, yQ3)

xP3 = 92537981321677359984282752681526321709848257167914581295182029482
yP3 = 64890706732687703644418730466418873818100958793359354306462487533
P3 = E_start(xP3, yP3)

# CONSISTENCY WITH R (NOT NEEDED BUT OK)
# Typo in magma, should be P3 and Q3
xR3 = 100985822528123790926639173045977043567721100360555352880158477546 + i*53651151307208479216007649001480928514654306829312202244997088803
assert (P3 - Q3)[0] == xR3

# ALICE'S PUBLIC KEY:
xPA = 78851325149093883127710876413676714831509071567295995722742563459 + i*21481241277029422048465721240261620296276964367410557936597294375
xQA = 109443172641179151776707590487317622581952970379528972130069645642 + i*3867221956981918915643365445875236474767408154492041549977730573
xRA = 24146826348939386714009375843953121070061474437444339669850030464 + i*9794094030587590044286487045494277248551007280921263840310791516

A = (1 - xPA*xQA - xPA*xRA - xQA*xRA)^2/(4*xPA*xQA*xRA) - xPA - xQA - xRA

yPA = sqrt(xPA^3 + A*xPA^2 + xPA)
yQA = sqrt(xQA^3 + A*xQA^2 + xQA)

# SMALL ERROR IN SIDH-spec.pdf, CORRECTED HERE
if xRA + xQA + xPA + A != (yQA + yPA)^2 / (xQA - xPA)^2:
    yQA = -yQA

# let's check:
EA = EllipticCurve(Fp2, [0,A,0,1,0])
EA.set_order((p+1)^2)
PA = EA(xPA, yPA)
QA = EA(xQA, yQA)

# BOB'S PUBLIC KEY:
xPB = 52037618715847826453371077000320105687598036562145407135988121710 + i*62945436285055860346151337655131657491042243534644871894809196747
xQB = 94057161062674597281795314311864564004565620907834550169224722966 + i*91420731496759657779126063859508682663377955903334296321639551249
xRB = 43790287819500432145214110821932450371863522319238208485657321972 + i*98694640376206066779482191725776091085259044935342665789389325446

B = (1 - xPB*xQB - xPB*xRB - xQB*xRB)^2/(4*xPB*xQB*xRB) - xPB - xQB - xRB

yPB = sqrt(xPB^3 + B*xPB^2 + xPB)
yQB = sqrt(xQB^3 + B*xQB^2 + xQB)

# SMALL ERROR IN SIDH-spec.pdf, CORRECTED HERE
if xRB + xQB + xPB + B != (yQB + yPB)^2 / (xQB - xPB)^2:
    yQB = -yQB

# let's check:
EB = EllipticCurve(Fp2, [0,B,0,1,0])
EB.set_order((p+1)^2)
PB = EB(xPB, yPB)
QB = EB(xQB, yQB)

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
    print(f"Prolonging with {prolong} steps. Current: {skB}")
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









