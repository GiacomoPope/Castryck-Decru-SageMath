def Coefficient(h, n):
    """
    Helper function to make things look similar!
    """
    return h[n]

def FromProdToJac(C, E, P_c, Q_c, P, Q, a):
    Fp2 = C.base()
    R.<x> = PolynomialRing(Fp2)

    P_c2 = 2^(a-1)*P_c
    Q_c2 = 2^(a-1)*Q_c
    P2 = 2^(a-1)*P
    Q2 = 2^(a-1)*Q

    alp1 = P_c2[0]
    alp2 = Q_c2[0]
    alp3 = (P_c2 + Q_c2)[0]
    bet1 = P2[0]
    bet2 = Q2[0]
    bet3 = (P2 + Q2)[0]
    a1 = (alp3 - alp2)^2/(bet3 - bet2) + (alp2 - alp1)^2/(bet2 - bet1) + (alp1 - alp3)^2/(bet1 - bet3)
    b1 = (bet3 - bet2)^2/(alp3 - alp2) + (bet2 - bet1)^2/(alp2 - alp1) + (bet1 - bet3)^2/(alp1 - alp3)
    a2 = alp1*(bet3 - bet2) + alp2*(bet1 - bet3) + alp3*(bet2 - bet1)
    b2 = bet1*(alp3 - alp2) + bet2*(alp1 - alp3) + bet3*(alp2 - alp1)
    Deltalp = (alp1 - alp2)^2*(alp1 - alp3)^2*(alp2 - alp3)^2
    Deltbet = (bet1 - bet2)^2*(bet1 - bet3)^2*(bet2 - bet3)^2

    A = Deltbet*a1/a2
    B = Deltalp*b1/b2

    h  = - (A*(alp2 - alp1)*(alp1 - alp3)*x^2 + B*(bet2 - bet1)*(bet1 - bet3)) 
    h *=   (A*(alp3 - alp2)*(alp2 - alp1)*x^2 + B*(bet3 - bet2)*(bet2 - bet1)) 
    h *=   (A*(alp1 - alp3)*(alp3 - alp2)*x^2 + B*(bet1 - bet3)*(bet3 - bet2))

    t1 = -(A/B)*b2/b1
    t2 = (bet1*(bet3 - bet2)^2/(alp3 - alp2) + bet2*(bet1 - bet3)^2/(alp1 - alp3) + bet3*(bet2 - bet1)^2/(alp2 - alp1))/b1
    s1 = -(B/A)*a2/a1
    s2 = (alp1*(alp3 - alp2)^2/(bet3 - bet2) + alp2*(alp1 - alp3)^2/(bet1 - bet3) + alp3*(alp2 - alp1)^2/(bet2 - bet1))/a1

    """
    We cannot do multivariate function fields, 
    so I followed 

    https://ask.sagemath.org/question/37584/multivariate-rational-function-field-of-rank-3-over-rational-field/
    
    To create a rational function field from a polynomial ring
    """
    Uff_poly.<u0, u1, v0, v1> = PolynomialRing(Fp2, 4)
    Uff = Uff_poly.fraction_field()
    Uff.inject_variables()

    A4.<U0, U1, V0, V1> = AffineSpace(Fp2, 4)
    # U.<U0, U1, V0, V1> = PolynomialRing(Fp2, 4)
    U = U0.parent()

    u0tilde = 1/u0
    u1tilde = u1/u0
    v0tilde = (u1*v0 - u0*v1)/u0^2
    v1tilde = (u1^2*v0 - u0*v0 - u0*u1*v1)/u0^2

    lamb1 = - ((Deltbet/A^3)*v1tilde)/(s1*u1tilde)
    lamb2 = - ((Deltalp/B^3)*v1)/(t1*u1)

    x1 = lamb1^2 + alp1 + alp2 + alp3 - s1*(u1tilde^2 - 2*u0tilde) - 2*s2
    y1 = -lamb1*(x1 - s2 + (u0tilde*v1tilde - u1tilde*v0tilde)*s1/v1tilde)

    x2 = lamb2^2 + bet1 + bet2 + bet3 - t1*(u1^2 - 2*u0) - 2*t2
    y2 = -lamb2*(x2 - t2 + (u0*v1 - u1*v0)*t1/v1)

    eq1 = U((x1 - P_c[0]).numerator())
    eq2 = U((y1 - P_c[1]).numerator())
    eq3 = U((x2 - P[0]).numerator())
    eq4 = U((y2 - P[1]).numerator())

    eq5  = 2*V0^2 - 2*V0*V1*U1 + V1^2*(U1^2 - 2*U0) 
    eq5 -= 2*Coefficient(h, 0)
    eq5 -= (-U1)*Coefficient(h, 1)
    eq5 -= (U1^2 - 2*U0)*Coefficient(h, 2)
    eq5 -= (-U1^3 + 3*U0*U1)*Coefficient(h, 3)
    eq5 -= (U1^4 - 4*U1^2*U0 + 2*U0^2)*Coefficient(h, 4)
    eq5 -= (-U1^5 + 5*U1^3*U0 - 5*U1*U0^2)*Coefficient(h, 5)
    eq5 -= (U1^6 - 6*U1^4*U0 + 9*U1^2*U0^2 - 2*U0^3)*Coefficient(h, 6)

    # UP TO HERE
    # Find a way to quickly find solutions!!

    V = A4.subscheme([eq1, eq2, eq3, eq4, eq5])

    # # point with zero coordinates probably correspond to "extra" solutions, we should be left with 4 sols
    # # (code may fail over small fields)

    # from sage.matrix.matrix2 import Matrix 
    # def resultant(f1, f2, var):
    #     return Matrix.determinant(f1.sylvester_matrix(f2, var))

    # tmp1 = resultant(eq3,  eq4, U0)
    # print(tmp1)
    # tmp2 = resultant(tmp1, eq5, V0)
    # print(tmp2)  
    # tmp3 = resultant(tmp2, eq1, V0)
    # print(tmp3)  
    # tmp4 = resultant(tmp3, eq2, V1)
    # print(tmp4)

    realsols = []
    for D in V.rational_points():
        print(D)
        Dseq = list(D)
        if not 0 in Dseq:
            realsols.append(Dseq)

    print(f"Number of inverse images found: {len(realsols)} (hopefully 4)")
    
    return False, False, False, False, False
    # TODO below (similar to above)

    J = Jacobian(HyperellipticCurve(h))
    sol = choice(realsols)
    D = J(x^2 + sol[1]*x + sol[0], sol[3]*x + sol[2])
    imPcP = 2*D

    # now for (Q_c, Q)

    eq1 = U(Numerator(x1 - Q_c[0]))
    eq2 = U(Numerator(y1 - Q_c[1]))
    eq3 = U(Numerator(x2 - Q[0]))
    eq4 = U(Numerator(y2 - Q[1]))

    V = A4.subscheme([eq1, eq2, eq3, eq4, eq5])
    realsols = []
    for D in V.rational_points():
        print(D)
        Dseq = list(D)
        if not 0 in Dseq:
            realsols.append(Dseq)

    print(f"Number of inverse images found: {len(realsols)} (hopefully 4)")
    sol = choice(realsols)
    D = J(x^2 + sol[1]*x + sol[0], sol[3]*x + sol[2])
    imQcQ = 2*D

    return h, imPcP[0], imPcP[1], imQcQ[0], imQcQ[1]

def FromJacToJac(h, D11, D12, D21, D22, a):
    R.<x> = h.parent()
    Fp2 = R.base()

    J = Jacobian(HyperellipticCurve(h))
    D1 = J(D11, D12)
    D2 = J(D21, D22)

    G1 = (2^(a-1)*D1)[0]
    G2 = (2^(a-1)*D2)[0]
    G3 = h / (G1*G2)

    delta = Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                              Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                              Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)])
    delta = delta.determinant()^(-1)

    H1 = delta*(derivative(G2)*G3 - G2*derivative(G3))
    H2 = delta*(derivative(G3)*G1 - G3*derivative(G1))
    H3 = delta*(derivative(G1)*G2 - G1*derivative(G2))

    hnew = H1*H2*H3
    Jnew = Jacobian(HyperellipticCurve(hnew))

    # now compute image points
    # first the point [D11, D12]:

    u0 = Coefficient(D11, 0)
    u1 = Coefficient(D11, 1)
    v0 = Coefficient(D12, 0)
    v1 = Coefficient(D12, 1)
    S.<x1,y1,y2,x2> = PolynomialRing(Fp2, 4)
    # TODO: How is this computed in sage?
    # pr = hom<S -> R | 0, 0, 0, x>

    eq1 = x1^2 + u1*x1 + u0
    eq2 = v1*x1 + v0 - y1
    # TODO: What is evaluate doing?
    # eq3 = Evaluate(G1, x1)*Evaluate(H1, x2) + Evaluate(G2, x1)*Evaluate(H2, x2)
    # eq4 = y1*y2 - Evaluate(G1, x1)*Evaluate(H1, x2)*(x1 - x2)
    # eq5 = y1^2 - Evaluate(h, x1)

    I = Ideal([eq1, eq2, eq3, eq4, eq5])
    G = GroebnerBasis(I) # last two are in non-reduced Mumford form: y2 + cubic(x2), quartic(x2)
    # TODO what is pr?
    # unew = pr(G[len(G)])
    # vnew = -pr(G[len(G)-1])
    # sanity check: (vnew^2 - hnew) mod unew
    # TODO what is pr?
    # imD1 = Jnew(pr(G[len(G)]), -pr(G[len(G)-1]))

    # now same for the point [D21, D22]:

    u0 = Coefficient(D21, 0)
    u1 = Coefficient(D21, 1)
    v0 = Coefficient(D22, 0)
    v1 = Coefficient(D22, 1)
    eq1 = x1^2 + u1*x1 + u0
    eq2 = v1*x1 + v0 - y1
    I = Ideal([eq1, eq2, eq3, eq4, eq5])
    G = GroebnerBasis(I)
    # TODO what is pr?
    # unew = pr(G[len(G)])
    # vnew = -pr(G[len(G)-1])
    # TODO what is pr?
    # imD2 = Jnew(pr(G[len(G)]), -pr(G[len(G)-1]))

    return hnew, imD1[0], imD1[1], imD2[0], imD2[1]

def Does22ChainSplit(C, E, P_c, Q_c, P, Q, a):
    Fp2 = C.base()
    # gluing step
    h, D11, D12, D21, D22 = FromProdToJac(C, E, P_c, Q_c, P, Q, a);
    # TODO
    return False
    # TODO
    print(f"order 2^{a-1} on hyp curve {h}")
    for i in range(1,a-2+1):
        h, D11, D12, D21, D22 = FromJacToJac(h, D11, D12, D21, D22, a-i)
        print(f"order 2^{a - i - 1} on hyp curve {h}")
    # now we are left with a quadratic splitting: is it singular?
    G1 = D11
    G2 = D21
    G3 = h # (G1*G2)
    print(G1, G2, G3)

    delta = Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                               Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                               Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)])
    delta = delta.determinant();
    return delta == 0

def OddCyclicSumOfSquares(n, factexpl, provide_own_fac):
    return NotImplemented

def Pushing3Chain(E, P, i):
    """
    Compute chain of isogenies quotienting 
    out a point P of order 3^i
    """
    Fp2 = E.base()
    R.<x> = PolynomialRing(Fp2)
    chain = []
    C = E
    remainingker = P
    # for j in [1..i] do
    for j in range(1, i+1):
        kerpol = x - (3^(i-j)*remainingker)[0]
        comp = EllipticCurveIsogeny(C, kerpol)
        C = comp.codomain()
        remainingker = comp(remainingker)
        chain.append(comp)
    return C, chain

def Pushing9Chain(E, P, i):
    """
    Compute chain of isogenies quotienting 
    out a point P of order 3^i
    """
    Fp2 = E.base()
    R.<x> = PolynomialRing(Fp2)
    chain = []
    C = E
    remainingker = P
    # for j in [1..i] do
    for j in range(1, i+1):
        # kerpol := &*[x - (k*9^(i-j)*remainingker)[1] : k in [1..4]];
        kerpol = x
        for k in range(1,4+1):
                kerpol -= (k*9^(i-j)*remainingker)[0]
        comp = EllipticCurveIsogeny(C, kerpol)
        C = comp.codomain()
        remainingker = comp(remainingker)
        chain.append(comp)
    return C, chain