/*
This code was downloaded from https://homes.esat.kuleuven.be/~wcastryc/ on 3 August 2022
*/

function FromProdToJac(C, E, P_c, Q_c, P, Q, a)
  Fp2 := BaseField(C);
  R<x> := PolynomialRing(Fp2);

  P_c2 := 2^(a-1)*P_c;
  Q_c2 := 2^(a-1)*Q_c;
  P2 := 2^(a-1)*P;
  Q2 := 2^(a-1)*Q;

  alp1 := P_c2[1];
  alp2 := Q_c2[1];
  alp3 := (P_c2 + Q_c2)[1];
  bet1 := P2[1];
  bet2 := Q2[1];
  bet3 := (P2 + Q2)[1];
  a1 := (alp3 - alp2)^2/(bet3 - bet2) + (alp2 - alp1)^2/(bet2 - bet1) + (alp1 - alp3)^2/(bet1 - bet3);
  b1 := (bet3 - bet2)^2/(alp3 - alp2) + (bet2 - bet1)^2/(alp2 - alp1) + (bet1 - bet3)^2/(alp1 - alp3);
  a2 := alp1*(bet3 - bet2) + alp2*(bet1 - bet3) + alp3*(bet2 - bet1);
  b2 := bet1*(alp3 - alp2) + bet2*(alp1 - alp3) + bet3*(alp2 - alp1);
  Deltalp := (alp1 - alp2)^2*(alp1 - alp3)^2*(alp2 - alp3)^2;
  Deltbet := (bet1 - bet2)^2*(bet1 - bet3)^2*(bet2 - bet3)^2;
  A := Deltbet*a1/a2;
  B := Deltalp*b1/b2;

  h := - (A*(alp2 - alp1)*(alp1 - alp3)*x^2 + B*(bet2 - bet1)*(bet1 - bet3)) *
         (A*(alp3 - alp2)*(alp2 - alp1)*x^2 + B*(bet3 - bet2)*(bet2 - bet1)) *
         (A*(alp1 - alp3)*(alp3 - alp2)*x^2 + B*(bet1 - bet3)*(bet3 - bet2));

  t1 := -(A/B)*b2/b1;
  t2 := (bet1*(bet3 - bet2)^2/(alp3 - alp2) + bet2*(bet1 - bet3)^2/(alp1 - alp3) + bet3*(bet2 - bet1)^2/(alp2 - alp1))/b1;
  s1 := -(B/A)*a2/a1;
  s2 := (alp1*(alp3 - alp2)^2/(bet3 - bet2) + alp2*(alp1 - alp3)^2/(bet1 - bet3) + alp3*(alp2 - alp1)^2/(bet2 - bet1))/a1;

  Uff<u0, u1, v0, v1> := FunctionField(Fp2, 4);
  A4<U0, U1, V0, V1> := AffineSpace(Fp2, 4); U := Parent(U0);

  u0tilde := 1/u0;
  u1tilde := u1/u0;
  v0tilde := (u1*v0 - u0*v1)/u0^2;
  v1tilde := (u1^2*v0 - u0*v0 - u0*u1*v1)/u0^2;

  lamb1 := - (Deltbet/A^3)*v1tilde/(s1*u1tilde);
  lamb2 := - (Deltalp/B^3)*v1/(t1*u1);

  x1 := lamb1^2 + alp1 + alp2 + alp3 - s1*(u1tilde^2 - 2*u0tilde) - 2*s2;
  y1 := -lamb1*(x1 - s2 + (u0tilde*v1tilde - u1tilde*v0tilde)*s1/v1tilde);

  x2 := lamb2^2 + bet1 + bet2 + bet3 - t1*(u1^2 - 2*u0) - 2*t2;
  y2 := -lamb2*(x2 - t2 + (u0*v1 - u1*v0)*t1/v1);

  eq1 := U ! Numerator(x1 - P_c[1]);
  eq2 := U ! Numerator(y1 - P_c[2]);
  eq3 := U ! Numerator(x2 - P[1]);
  eq4 := U ! Numerator(y2 - P[2]);
  eq5 := 2*V0^2 - 2*V0*V1*U1 + V1^2*(U1^2 - 2*U0)
         - 2*Coefficient(h, 0)
         - (-U1)*Coefficient(h, 1)
         - (U1^2 - 2*U0)*Coefficient(h, 2)
         - (-U1^3 + 3*U0*U1)*Coefficient(h, 3)
         - (U1^4 - 4*U1^2*U0 + 2*U0^2)*Coefficient(h, 4)
         - (-U1^5 + 5*U1^3*U0 - 5*U1*U0^2)*Coefficient(h, 5)
         - (U1^6 - 6*U1^4*U0 + 9*U1^2*U0^2 - 2*U0^3)*Coefficient(h, 6);

  V := Scheme(A4, [eq1, eq2, eq3, eq4, eq5]);

  // point with zero coordinates probably correspond to "extra" solutions, we should be left with 4 sols
  // (code may fail over small fields)

  realsols := [];
  for D in Points(V) do
    Dseq := Eltseq(D);
    if not 0 in Dseq then
      realsols cat:= [Dseq];
    end if;
  end for;

  // print "Number of inverse images found:", #realsols, "(hopefully 4)";

  J := Jacobian(HyperellipticCurve(h));
  sol := Random(realsols);
  D := elt<J | x^2 + sol[2]*x + sol[1], sol[4]*x + sol[3]>;
  imPcP := 2*D;

  // now for (Q_c, Q)

  eq1 := U ! Numerator(x1 - Q_c[1]);
  eq2 := U ! Numerator(y1 - Q_c[2]);
  eq3 := U ! Numerator(x2 - Q[1]);
  eq4 := U ! Numerator(y2 - Q[2]);
  V := Scheme(A4, [eq1, eq2, eq3, eq4, eq5]);
  realsols := [];
  for D in Points(V) do
    Dseq := Eltseq(D);
    if not 0 in Dseq then
      realsols cat:= [Dseq];
    end if;
  end for;
  // print "Number of inverse images found:", #realsols, "(hopefully 4)";
  sol := Random(realsols);
  D := elt<J | x^2 + sol[2]*x + sol[1], sol[4]*x + sol[3]>;
  imQcQ := 2*D;

  return h, imPcP[1], imPcP[2], imQcQ[1], imQcQ[2];
end function;

function FromJacToJac(h, D11, D12, D21, D22, a)
  R<x> := Parent(h);
  Fp2 := BaseRing(R);

  J := Jacobian(HyperellipticCurve(h));
  D1 := elt<J | D11, D12>;
  D2 := elt<J | D21, D22>;

  G1 := (2^(a-1)*D1)[1];
  G2 := (2^(a-1)*D2)[1];
  G3 := h div (G1*G2);

  delta := Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                              Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                              Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)]);
  delta := Determinant(delta)^(-1);

  H1 := delta*(Derivative(G2)*G3 - G2*Derivative(G3));
  H2 := delta*(Derivative(G3)*G1 - G3*Derivative(G1));
  H3 := delta*(Derivative(G1)*G2 - G1*Derivative(G2));

  hnew := H1*H2*H3;
  Jnew := Jacobian(HyperellipticCurve(hnew));

  // now compute image points
  // first the point [D11, D12]:

  u0 := Coefficient(D11, 0);
  u1 := Coefficient(D11, 1);
  v0 := Coefficient(D12, 0);
  v1 := Coefficient(D12, 1);
  S<x1,y1,y2,x2> := PolynomialRing(Fp2, 4);
  pr := hom<S -> R | 0, 0, 0, x>;

  eq1 := x1^2 + u1*x1 + u0;
  eq2 := v1*x1 + v0 - y1;
  eq3 := Evaluate(G1, x1)*Evaluate(H1, x2) + Evaluate(G2, x1)*Evaluate(H2, x2);
  eq4 := y1*y2 - Evaluate(G1, x1)*Evaluate(H1, x2)*(x1 - x2);
  eq5 := y1^2 - Evaluate(h, x1);
  I := Ideal([eq1, eq2, eq3, eq4, eq5]);
  G := GroebnerBasis(I); // last two are in non-reduced Mumford form: y2 + cubic(x2), quartic(x2)
  unew := pr(G[#G]);
  vnew := -pr(G[#G-1]);
  // sanity check: (vnew^2 - hnew) mod unew;
  imD1 := elt<Jnew | pr(G[#G]), -pr(G[#G-1])>;

  // now same for the point [D21, D22]:

  u0 := Coefficient(D21, 0);
  u1 := Coefficient(D21, 1);
  v0 := Coefficient(D22, 0);
  v1 := Coefficient(D22, 1);
  eq1 := x1^2 + u1*x1 + u0;
  eq2 := v1*x1 + v0 - y1;
  I := Ideal([eq1, eq2, eq3, eq4, eq5]);
  G := GroebnerBasis(I);
  unew := pr(G[#G]);
  vnew := -pr(G[#G-1]);
  imD2 := elt<Jnew | pr(G[#G]), -pr(G[#G-1])>;

  return hnew, imD1[1], imD1[2], imD2[1], imD2[2];
end function;

function Does22ChainSplit(C, E, P_c, Q_c, P, Q, a)
  Fp2 := BaseField(C);
  // gluing step
  h, D11, D12, D21, D22 := FromProdToJac(C, E, P_c, Q_c, P, Q, a);
  // print "order 2^", a-1, "on hyp curve", h;
  for i in [1..a-2] do
    h, D11, D12, D21, D22 := FromJacToJac(h, D11, D12, D21, D22, a-i);
    // print "order 2^", a - i - 1, "on hyp curve", h;
  end for;
  // now we are left with a quadratic splitting: is it singular?
  G1 := D11;
  G2 := D21;
  G3 := h div (G1*G2);
  // print G1, G2, G3;
  delta := Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                              Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                              Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)]);
  delta := Determinant(delta);
  return delta eq 0;
end function;

function OddCyclicSumOfSquares(n, factexpl, provide_own_fac)
  if provide_own_fac then
    fac := factexpl;
  else
    fac := Factorization(n);
  end if;
  if {f[1] mod 4 : f in fac} ne {1} then
    return false, 0, 0;
  else
    Z<I> := GaussianIntegers();
    prod := 1;
    for f in fac do
      p := f[1];
      repeat
        z := Random(GF(p));
      until z^((p-1) div 2) eq -1;
      z := Integers() ! (z^((p-1) div 4)); // square root of -1
      prod *:= GCD(Z ! p, Z ! z + I)^f[2];
    end for;
    u := [Integers() ! e : e in Eltseq(prod)];
    if IsOdd(u[1]) then u1 := u[1]; u2 := u[2] div 2; else u1 := u[2]; u2 := u[1] div 2; end if;
    return true, Abs(u1), Abs(u2);
  end if;
end function;


function Pushing3Chain(E, P, i) // compute chain of isogenies quotienting out a point P of order 3^i
  Fp2 := BaseField(E);
  R<x> := PolynomialRing(Fp2);
  chain := [];
  C := E;
  remainingker := P;
  for j in [1..i] do
    kerpol := x - (3^(i-j)*remainingker)[1];
    C, comp := IsogenyFromKernel(C, kerpol);
    remainingker := comp(remainingker);
    chain cat:=[comp];
  end for;
  return C, chain;
end function;


function Pushing9Chain(E, P, i) // compute chain of isogenies quotienting out a point P of order 9^i (obsolete)
  Fp2 := BaseField(E);
  R<x> := PolynomialRing(Fp2);
  chain := [];
  C := E;
  remainingker := P;
  for j in [1..i] do
    kerpol := &*[x - (k*9^(i-j)*remainingker)[1] : k in [1..4]];
    C, comp := IsogenyFromKernel(C, kerpol);
    remainingker := comp(remainingker);
    chain cat:=[comp];
  end for;
  return C, chain;
end function;
