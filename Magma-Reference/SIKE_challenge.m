/*
This code was downloaded from https://homes.esat.kuleuven.be/~wcastryc/ on 3 August 2022
*/

load "richelot_aux.m";
load "uvtable.m";

a := 110;
b := 67;

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

xQ2 := 86445170414599058485968715662346922796016566052081315743781191908 + I*40452781127453632342062231901634407992350274206219311799642330230;
yQ2 := 65598021239445364766945421641383356785387864046503211101901536990 + I*4185188043442820950612067662642469583838730096408235972690714750;
Q2 := E_start ! [xQ2, yQ2];

xP2 := 37927755131506746411537091732732279589710442710570639382902295344 + I*37017495259074475516858568039563512752905548684528955375066616976;
yP2 := 47930070613074499206421889519713097155043539387599009924307417901 + I*7176302698603683442486948733484190513242470057553561741081582098;
P2 := E_start ! [xP2, yP2];

// CONSISTENCY WITH R (NOT NEEDED BUT OK)
// xR2 := 46724562430016373759602265681716628459894400021233033976074511127 + I*62869247191516167619032311426583862201277153531817575978673280654;
// (P2 - Q2)[1] eq xR2;

xQ3 := 47793138785202808493310870472269548499693686682045542672866243781;
yQ3 := I*59233392251697228025115695338640088714243311845504632846426147327;
Q3 := E_start ! [xQ3, yQ3];

xP3 := 92537981321677359984282752681526321709848257167914581295182029482;
yP3 := 64890706732687703644418730466418873818100958793359354306462487533;
P3 := E_start ! [xP3, yP3];

// CONSISTENCY WITH R (NOT NEEDED BUT OK)
// xRB := 100985822528123790926639173045977043567721100360555352880158477546 + I*53651151307208479216007649001480928514654306829312202244997088803;
// (PB - QB)[1] eq xRB;

// ALICE'S PUBLIC KEY:

xPA := 78851325149093883127710876413676714831509071567295995722742563459 + I*21481241277029422048465721240261620296276964367410557936597294375;
xQA := 109443172641179151776707590487317622581952970379528972130069645642 + I*3867221956981918915643365445875236474767408154492041549977730573;
xRA := 24146826348939386714009375843953121070061474437444339669850030464 + I*9794094030587590044286487045494277248551007280921263840310791516;

A := (1 - xPA*xQA - xPA*xRA - xQA*xRA)^2/(4*xPA*xQA*xRA) - xPA - xQA - xRA;

yPA := Sqrt(xPA^3 + A*xPA^2 + xPA);
yQA := Sqrt(xQA^3 + A*xQA^2 + xQA);

if xRA + xQA + xPA + A ne (yQA + yPA)^2 / (xQA - xPA)^2 then yQA := -yQA; end if; // SMALL ERROR IN SIDH-spec.pdf, CORRECTED HERE

// let's check:

EA := EllipticCurve(x^3 + A*x^2 + x);
PA := EA ! [xPA, yPA];
QA := EA ! [xQA, yQA];

// BOB'S PUBLIC KEY:

xPB := 52037618715847826453371077000320105687598036562145407135988121710 + I*62945436285055860346151337655131657491042243534644871894809196747;
xQB := 94057161062674597281795314311864564004565620907834550169224722966 + I*91420731496759657779126063859508682663377955903334296321639551249;
xRB := 43790287819500432145214110821932450371863522319238208485657321972 + I*98694640376206066779482191725776091085259044935342665789389325446;

B := (1 - xPB*xQB - xPB*xRB - xQB*xRB)^2/(4*xPB*xQB*xRB) - xPB - xQB - xRB;

yPB := Sqrt(xPB^3 + B*xPB^2 + xPB);
yQB := Sqrt(xQB^3 + B*xQB^2 + xQB);

if xRB + xQB + xPB + B ne (yQB + yPB)^2 / (xQB - xPB)^2 then yQB := -yQB; end if; // SMALL ERROR IN SIDH-spec.pdf, CORRECTED HERE

// let's check:

EB := EllipticCurve(x^3 + B*x^2 + x);
PB := EB ! [xPB, yPB];
QB := EB ! [xQB, yQB];

// ===================================
// =====  ATTACK  ====================
// ===================================

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
  tauhatkernel := 3^bi*P3 + (j[1] + 3*j[2])*3^bi*Q3;
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
