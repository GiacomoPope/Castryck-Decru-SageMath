set_verbose(-1)

def Coefficient(h, n):
    """
    Helper function to make things look similar!
    """
    assert h.denominator().is_one()
    h = h.numerator()
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

    H = HyperellipticCurve(h)
    J = H.jacobian()

    # We need the image of (P_c, P) and (Q_c, Q) in J
    # The image of (P_c, P) is the image of P_c as a divisor on H
    # plus the image of P as a divisor on H.
    # This allows for direct computation without solving equations
    # as in Castryck-Decru's paper.

    # The projection maps are:
    # H->C: (xC = s1/x²+s2, yC = (Deltbet/A³)(y/x³))
    # so we compute Mumford coordinates of the divisor f^-1(P_c): a(x), y-b(x)
    xPc, yPc = P_c.xy()
    mumPc = [x^2 - s1 / (xPc - s2), yPc * x^3 * A^3 / Deltbet]
    JPc = J(mumPc)
    # same for Q_c
    xQc, yQc = Q_c.xy()
    mumQc = [x^2 - s1 / (xQc - s2), yQc * x^3 * A^3 / Deltbet]
    JQc = J(mumQc)

    # Same for E
    # H->E: (xE = t1 x² + t2, yE = (Deltalp/B³)y)
    xP, yP = P.xy()
    JP = J([t1* x^2 + t2 - xP, R(yP * B^3 / Deltalp)])
    xQ, yQ = Q.xy()
    JQ = J([t1* x^2 + t2 - xQ, R(yQ * B^3 / Deltalp)])

    imPcP = JP + JPc
    imQcQ = JQ + JQc

    # Validate result
    # For debugging
    # def projC(_x, _y):
    #     return (s1 / _x^2 + s2, Deltbet / A^3 * _y / _x^3)
    # def projE(_x, _y):
    #     return (t1 * _x^2 + t2, Deltalp / B^3 * _y)
    # Fp4 = Fp2.extension(2)
    # E4 = E.change_ring(Fp4)
    # C4 = C.change_ring(Fp4)
    # divP = [(xr, imPcP[1](xr)) for xr, _ in imPcP[0].roots(Fp4)]
    # assert 2*E4(P) == sum(E4(*projE(*pt)) for pt in divP)
    # assert 2*C4(P_c) == sum(C4(*projC(*pt)) for pt in divP)
    # divQ = [(xr, imQcQ[1](xr)) for xr, _ in imQcQ[0].roots(Fp4)]
    # assert 2*E4(Q) == sum(E4(*projE(*pt)) for pt in divQ)
    # assert 2*C4(Q_c) == sum(C4(*projC(*pt)) for pt in divQ)

    return h, imPcP[0], imPcP[1], imQcQ[0], imQcQ[1]

def test_FromProdToJac():
    # Choose some supersingular curves and 2^a torsion generators.
    p = 2**61 - 1
    assert p.is_prime()
    k = GF(p^2)
    E = EllipticCurve(k, [-1, 0])
    assert E.is_supersingular()
    assert E.order() == 2**122
    C = EllipticCurve(k, [1, 0])
    assert C.is_supersingular()
    assert C.order() == 2**122
    a = 61
    Pc, Qc = C.gens()
    P, Q = E.gens()
    wc = Pc.weil_pairing(Qc, 2^a)
    we = P.weil_pairing(Q, 2^a)
    # make it an anti-isometry
    k = we.log(wc)
    Q = -pow(k, -1, 2^a) * Q
    assert P.weil_pairing(Q, 2^a) * Pc.weil_pairing(Qc, 2^a) == 1
    return FromProdToJac(C, E, Pc, Qc, P, Q, a)

#test_FromProdToJac()

def FromJacToJac(h, D11, D12, D21, D22, a, power=None):
    # power is an optional precomputed tuple (l, 2^l D1, 2^l D2)
    # where l < a
    R = h.parent()
    Fp2 = R.base()

    J = HyperellipticCurve(h).jacobian()
    D1 = J(D11, D12)
    D2 = J(D21, D22)

    next_power = None
    if power is None:
        # Precompute some power of D1, D2 to save computations later.
        if a >= 16:
            _D1 = 2^(a-8)*D1
            _D2 = 2^(a-8)*D2
            G1, _ = 2^7 * _D1
            G2, _ = 2^7 * _D2
            next_power = (a-8, _D1, _D2)
        else:
            G1, _ = 2^(a-1) * D1
            G2, _ = 2^(a-1) * D2
    else:
        (l, _D1, _D2) = power
        if l < a-1 and a >= 16:
            next_power = power
        G1, _ = 2^(a-1-l)*J(*_D1)
        G2, _ = 2^(a-1-l)*J(*_D2)

    #assert 2^a*D1 == 0
    #assert 2^a*D2 == 0
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0

    delta = Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                              Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                              Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)])
    delta = delta.determinant()^(-1)

    H1 = delta*(derivative(G2)*G3 - G2*derivative(G3))
    H2 = delta*(derivative(G3)*G1 - G3*derivative(G1))
    H3 = delta*(derivative(G1)*G2 - G1*derivative(G2))

    hnew = H1*H2*H3

    # Now compute image points: Richelot isogeny is defined by the degree 2
    # correspondence:
    # g1(x1) h1(x2) + g2(x1) h2(x2) = 0
    # or y1 y2 = g1(x1) h1(x2) (x1 - x2)

    R2.<y,x,xa,xb> = PolynomialRing(Fp2, order="lex")
    def image(D):
        # Let x^2 + u1 x + u0 = U(x) = 0
        #     y =   v1 x + v0
        # be Mumford cordinates of D.
        # Introduce formal coordinates xa, xb for the roots of U.
        # We should compute the image of (xa, ya) and (xb, yb) through
        # the correspondance.
        U, V = D
        u1, u0 = U[1], U[0]
        assert U[2] == 1
        # Perform manually-assisted elimination:
        # Compute g1 % U and g2 % U
        #    (g11 x1a + g10) h1(x2) + (g21 x1a + g20) h2(x2) = 0
        # OR (g11 x1b + g10) h1(x2) + (g21 x1b + g20) h2(x2) = 0
        # Multiply these equations:
        G1red = G1 % U
        G2red = G2 % U
        Qx = (G1red.subs(x=xa) * H1 + G2red.subs(x=xa) * H2) * \
             (G1red.subs(x=xb) * H1 + G2red.subs(x=xb) * H2)
        # This is a polynomial in:
        # xa*xb = u0
        # xa + xb = -u1
        Px = u0 * Qx.coefficient({xa: 1, xb: 1}) \
           - u1 * Qx.coefficient({xa: 1, xb: 0}) \
           + Qx.coefficient({xa: 0, xb: 0})
        Px = Px.univariate_polynomial()

        # Similarly
        #    (v1 x1a + v0) y2 = (g11 x1a + g10) h1(x2) (x1a - x2)
        # OR (v1 x1b + v0) y2 = (g11 x1b + g10) h1(x2) (x1b - x2)
        # Multiply and reduce:
        Qy = (V(xa) * y - G1red(xa) * H1 * (xa - x)) * \
             (V(xb) * y - G1red(xb) * H1 * (xb - x))
        # Reduce symmetric functions:
        # xa^2 xb^2 = u0^2
        # xa^2 xb + xa xb^2 = -u0*u1
        # xa^2 + xb^2 = u1**2 - 2*u0
        # etc.
        Py = u0^2 * Qy.coefficient({xa: 2, xb: 2}) \
           - u0 * u1 * Qy.coefficient({xa: 2, xb: 1}) \
           + (u1**2 - 2*u0) * Qy.coefficient({xa: 2, xb: 0}) \
           + u0 * Qy.coefficient({xa: 1, xb: 1}) \
           - u1 * Qy.coefficient({xa: 1, xb: 0}) \
           + Qy.coefficient({xa: 0, xb: 0})
        # Py = A y^2 + B y + C
        #    = B y + (A hnew + C) = 0
        A = Py.coefficient({y: 2})
        assert A.is_constant()
        A = A.constant_coefficient()
        B = Py.coefficient({y: 1}).univariate_polynomial()
        C = Py.coefficient({y: 0}).univariate_polynomial()

        _, Binv, _ = B.xgcd(Px)
        Py = (- Binv * (A * hnew + C)) % Px

        # Inlined Cantor reduction
        # Px has degree 4, Py has degree 3
        Dx = ((hnew - Py ** 2) // Px).monic()
        Dy = (-Py) % Dx
        return (Dx, Dy)

    imD1 = image(D1)
    imD2 = image(D2)
    if next_power:
        l, _D1, _D2 = next_power
        next_power = (l, image(_D1), image(_D2))
    return hnew, imD1[0], imD1[1], imD2[0], imD2[1], next_power

def test_FromJacToJac():
    print("test_FromJacToJac")
    h, D11, D12, D21, D22 = test_FromProdToJac()
    FromJacToJac(h, D11, D12, D21, D22, 60)

#test_FromJacToJac()

def Does22ChainSplit(C, E, P_c, Q_c, P, Q, a):
    Fp2 = C.base()
    # gluing step
    h, D11, D12, D21, D22 = FromProdToJac(C, E, P_c, Q_c, P, Q, a);
    next_power = None
    # print(f"order 2^{a-1} on hyp curve {h}")
    for i in range(1,a-2+1):
        h, D11, D12, D21, D22, next_power = FromJacToJac(h, D11, D12, D21, D22, a-i, power=next_power)
        # print(f"order 2^{a - i - 1} on hyp curve {h}")
    # now we are left with a quadratic splitting: is it singular?
    G1 = D11
    G2 = D21
    G3 = h // (G1*G2)
    # print(G1, G2, G3)

    delta = Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                               Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                               Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)])
    delta = delta.determinant();
    return delta == 0

def test_ChainSplit():
    print("test_ChainSplit")
    p = 2**61 - 1
    k = GF(p^2)
    E = EllipticCurve(k, [-1, 0])
    C = EllipticCurve(k, [1, 0])
    a = 61
    Pc, Qc = C.gens()
    P, Q = E.gens()
    wc = Pc.weil_pairing(Qc, 2^a)
    we = P.weil_pairing(Q, 2^a)
    # make it an anti-isometry
    k = we.log(wc)
    Q = -pow(k, -1, 2^a) * Q
    assert P.weil_pairing(Q, 2^a) * Pc.weil_pairing(Qc, 2^a) == 1
    assert not Does22ChainSplit(C, E, Pc, Qc, P, Q, 61)

#test_ChainSplit()

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
