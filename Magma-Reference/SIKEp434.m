/*
This code was downloaded from https://homes.esat.kuleuven.be/~wcastryc/ on 3 August 2022
*/

load "richelot_aux.m";
load "uvtable.m";

a := 216;
b := 137;

p := 2^a*3^b - 1;
Fp2<I> := GF(p, 2);
assert I^2 eq -1;
R<x> := PolynomialRing(Fp2);

E_start := EllipticCurve(x^3 + 6*x^2 + x);

// NAIVE GENERATION OF AUTOMORPHISM 2i
E1728, phi := IsogenyFromKernel(E_start, x);
for iota in Automorphisms(E1728) do
  P := Random(E1728);
  if iota(iota(P)) eq -P then
    two_i := phi*iota*DualIsogeny(phi);
    break;
  end if;
end for;

infty := E_start ! 0;

repeat
  P2 := 3^b*Random(E_start);
until 2^(a-1)*P2 ne infty;
repeat
  Q2 := 3^b*Random(E_start);
until WeilPairing(P2, Q2, 2^a)^(2^(a-1)) ne 1;

repeat
  P3 := 2^a*Random(E_start);
until 3^(b-1)*P3 ne infty;
repeat
  Q3 := 2^a*Random(E_start);
until WeilPairing(P3, Q3, 3^b)^(3^(b-1)) ne 1;

// generate challenge key

Bobskey := Random(3^b);

EB, chain := Pushing3Chain(E_start, P3 + Bobskey*Q3, b);
PB := P2; for c in chain do PB := c(PB); end for;
QB := Q2; for c in chain do QB := c(QB); end for;

skB := []; // DIGITS IN EXPANSION OF BOB'S SECRET KEY
// kappas := []; // CHAIN OF ISOGENIES // is eigenlijk niet nodig, merkwaardig genoeg!

print "If all goes well then the following digits should be found", Intseq(Bobskey, 3);

tim := Cputime();

skB := []; // TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY

// gathering the alpha_i, u, v from table

expdata := [[0, 0, 0] : i in [1..b-3]];
for i in [1..b-3] do
  if IsOdd(b-i) then
    index := (b - i + 1) div 2;
    exp := uvtable[index][2];
    if exp le a then
      u := uvtable[index][3];
      v := uvtable[index][4];
      expdata[i] := [exp, u, v];
    end if;
  end if;
end for;

// gather digits until beta_1

bet1 := 0;
repeat
  bet1 +:= 1;
until expdata[bet1][1] ne 0;

ai := expdata[bet1][1];
u := expdata[bet1][2];
v := expdata[bet1][3];

print "Determination of first", bet1, "ternary digits. We are working with 2^"*IntegerToString(ai)*"-torsion.";

bi := b - bet1;
alp := a - ai;

for j in CartesianPower([0,1,2], bet1) do
  print "Testing digits", &*[IntegerToString(j[k])*" " : k in [1..bet1]];
  tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*j[k] : k in [1..bet1]])*3^bi*Q3;
  tauhatkernel_distort := u*tauhatkernel + v*two_i(tauhatkernel);
  C, tau_tilde := Pushing3Chain(E_start, tauhatkernel_distort, bet1);
  P_c := u*P2 + v*two_i(P2); for taut in tau_tilde do P_c := taut(P_c); end for;
  Q_c := u*Q2 + v*two_i(Q2); for taut in tau_tilde do Q_c := taut(Q_c); end for;
  if j eq <2 : k in [1..bet1]> or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai) then
    print "Glue-and-split! These are most likely the secret digits.";
    skB cat:= [j[k] : k in [1..bet1]];
    break;
  end if;
end for;

// now compute longest prolongation of Bob's isogeny that may be needed

length := 1;
max_length := 0;
for i in [bet1+1..b-3] do
  if expdata[i][1] ne 0 then
    max_length := Max(length, max_length);
    length := 0;
  else
    length +:= 1;
  end if;
end for;

repeat
  K := 2^a*3^(b - max_length)*Random(EB);
until Order(K) eq 3^max_length;

repeat
  alternativeK := 2^a*3^(b - max_length)*Random(EB);
until WeilPairing(K, alternativeK, 3^max_length)^(3^(max_length - 1)) ne 1;

_, EBprolong := Pushing3Chain(EB, K, max_length);

// gather next digit and change prolongation if needed

i := bet1 + 1;
bi := b - i;

print "Determination of the "*IntegerToString(i)*"th ternary digit. We are working with 2^"*IntegerToString(ai)*"-torsion.";
prolong := 1;
print "Prolonging with 1 steps.";
endPB := EBprolong[1](PB);
endQB := EBprolong[1](QB);
endEB := Codomain(EBprolong[1]);

positives := [];

for j in [0..2] do
  print "Testing digit", j;
  tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
  tauhatkernel_distort := u*tauhatkernel + v*two_i(tauhatkernel);
  C, tau_tilde := Pushing3Chain(E_start, tauhatkernel_distort, i);
  P_c := u*P2 + v*two_i(P2); for taut in tau_tilde do P_c := taut(P_c); end for;
  Q_c := u*Q2 + v*two_i(Q2); for taut in tau_tilde do Q_c := taut(Q_c); end for;
  if Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai) then
    print "Glue-and-split!";
    positives cat:=[j];
    if #positives gt 1 then
      break;
    end if;
  end if;
end for;

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

// now gather all remaining digits, except for last three

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
for i in [0..3^3] do
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
