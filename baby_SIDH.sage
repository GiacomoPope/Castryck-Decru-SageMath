load('richelot_aux.sage')
load('uvtable.sage')

a = 33
b = 19

p = 2^a*3^b - 1
Fp2.<i> = GF(p^2, modulus=x^2+1)
assert i^2 == -1
R.<x> = PolynomialRing(Fp2)

E_start = EllipticCurve(Fp2, [0,6,0,1,0])

phi = EllipticCurveIsogeny(E_start, x)
E1728 = phi.codomain()

for iota in E1728.automorphisms():
    P = E1728.random_point()
    if iota(iota(P)) == -P:
        # two_i = phi.post_compose(iota).post_compose(phi.dual())
        two_i = -phi.dual()*iota*phi
        break


infty = E_start(0)

P2x = 3001978773937227303*i + 844826078775277193
P2y = 7016480731488920513*i + 3336734539702591894
P2 = E_start(P2x, P2y)

Q2x = 9714335965426121848*i + 3634229810067794918
Q2y = 6843179620215919990*i + 283325142855416748
Q2 = E_start(Q2x, Q2y)

P3x = 856716429329350834*i + 3275086448477871897
P3y = 789318291276569898*i + 5513314229840294380
P3 = E_start(P3x, P3y)

Q3x = 9698695942032140636*i + 4611223447894666827
Q3y = 4739425428244553293*i + 434003219133716652
Q3 = E_start(Q3x, Q3y)

assert 2^(a-1)*P2 != infty
assert 3^(b-1)*P3 != infty
assert P2.weil_pairing(Q2, 2^a)^(2^(a-1)) != 1
assert P3.weil_pairing(Q3, 3^b)^(3^(b-1)) != 1

# generate challenge key
# Bobskey = randint(0,3^b)
Bobskey = 15002860

EB, chain = Pushing3Chain(E_start, P3 + Bobskey*Q3, b);
PB = P2
for c in chain:
    PB = c(PB)
QB = Q2 
for c in chain:
    QB = c(QB)

skB = [] # DIGITS IN EXPANSION OF BOB'S SECRET KEY

print(f"If all goes well then the following digits should be found: {Integer(Bobskey).digits(base=3)      }")



# ===================================
# =====  ATTACK  ====================
# ===================================
import time
tim = time.time()

skB = [] # TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY

# gathering the alpha_i, u, v from table

expdata = [[0, 0, 0] for _ in range(b-3)]
# for i in [1..b-3] do
for i in range(1,b-2):
    # if IsOdd(b-i) then
    if (b-i)%2 == 1:
        index = (b - i + 1) // 2
        exp = uvtable[index-1][1];
        if exp <= a:
            u = uvtable[index-1][2];
            v = uvtable[index-1][3];
            expdata[i-1] = [exp, u, v];

# gather digits until beta_1
bet1 = 0
while expdata[bet1][0] == 0:
    bet1 += 1
bet1 += 1

ai = expdata[bet1-1][0];
u  = expdata[bet1-1][1];
v  = expdata[bet1-1][2];

print(f"Determination of first {bet1} ternary digits. We are working with 2^{ai}-torsion.")

bi = b - bet1
alp = a - ai

# for j in CartesianPower([0,1,2], bet1) do
from itertools import product
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

print(skB)





"""
# ===================================
# === Remaining magma to convert! ===
# ===================================

if #positives eq 1 then
  print "Most likely good prolongation and the secret digit is", positives[1];
  skB cat:= [positives[1]];
  next_i := i + 1;
else
  print "All glue-and-splits, must be bad prolongation: changing it and redoing this digit.";
  _, EBprolong := Pushing3Chain(EB, alternativeK, max_length);
  next_i := i;
  prolong := 0;
end if;

// now gather all remaining digits, except for last three (we close that gap by trial and error)

for i in [next_i..b-3] do
  bi := b - i;
  if expdata[i][1] ne 0 then
    ai := expdata[i][1];
    alp := a - ai;
    u := expdata[i][2];
    v := expdata[i][3];
    prolong := 0;
  else
    prolong +:= 1;
  end if;
  print "Determination of the "*IntegerToString(i)*"th ternary digit. We are working with 2^"*IntegerToString(ai)*"-torsion.";
  print "Prolonging with", prolong, "steps.";
  endPB := PB;
  endQB := QB;
  endEB := EB;
  for j in [1..prolong] do
    endPB := EBprolong[j](endPB);
    endQB := EBprolong[j](endQB);
    if j eq prolong then
      endEB := Codomain(EBprolong[j]);
    end if;
  end for;

  for j in [0..2] do
    print "Testing digit", j;
    tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
    tauhatkernel_distort := u*tauhatkernel + v*two_i(tauhatkernel);
    C, tau_tilde := Pushing3Chain(E_start, tauhatkernel_distort, i);
    P_c := u*P2 + v*two_i(P2); for taut in tau_tilde do P_c := taut(P_c); end for;
    Q_c := u*Q2 + v*two_i(Q2); for taut in tau_tilde do Q_c := taut(Q_c); end for;
    if j eq 2 or Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai) then
      print "Glue-and-split! This is most likely the secret digit.";
      skB cat:= [j];
      break;
    end if;
  end for;

end for;


key := &+[skB[i]*3^(i-1) : i in [1..b-3]];

// bridge last safety gap

tim2 := Cputime();

found := false;
for i in [0..3^5] do
  bobskey := key + i*3^(b-3);
  bobscurve := Pushing3Chain(E_start, P3 + bobskey*Q3, b);
  if jInvariant(bobscurve) eq jInvariant(EB) then
    found := true;
    break;
  end if;
end for;

print "Bridging last gap took", Cputime(tim2);

if found then
  print "Bob's secret key revealed as", bobskey;
else
  print "Something went wrong.";
end if;

print "Altogether this took", Cputime(tim), "seconds.";
"""







