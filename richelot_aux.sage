set_verbose(-1)

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

class RichelotCorr:
    """
    The Richelot correspondance between hyperelliptic
    curves y²=g1*g2*g3 and y²=h1*h2*h3=hnew(x)

    It is defined by equations:
        g1(x1) h1(x2) + g2(x1) h2(x2) = 0
    and y1 y2 = g1(x1) h1(x2) (x1 - x2)

    Given a divisor D in Mumford coordinates:
        U(x) = x^2 + u1 x + u0 = 0
        y = V(x) = v1 x + v0
    Let xa and xb be the symbolic roots of U.
    Let s, p by the sum and product (s=-u1, p=u0)

    Then on x-coordinates, the image of D is defined by equation:
        (g1(xa) h1(x) + g2(xa) h2(x))
      * (g1(xb) h1(x) + g2(xb) h2(x))
    which is a symmetric function of xa and xb.
    This is a non-reduced polynomial of degree 4.

    Write gred = g-U = g1*x + g0
    then gred(xa) gred(xb) = g1^2*p + g1*g0*s + g0^2
    and  g1red(xa) g2red(xb) + g1red(xb) g2red(xa)
       = 2 g11 g21 p + (g11*g20+g10*g21) s + 2 g10*g20

    On y-coordinates, the image of D is defined by equations:
           V(xa) y = Gred1(xa) h1(x) (xa - x)
        OR V(xb) y = Gred1(xb) h1(x) (xb - x)
    If we multiply:
    * y^2 has coefficient V(xa)V(xb)
    * y has coefficient h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
      (x-degree 3)
    * 1 has coefficient Gred1(xa) Gred1(xb) h1(x)^2 (x-xa)(x-xb)
                      = Gred1(xa) Gred1(xb) h1(x)^2 U(x)
      (x-degree 4)
    """
    def __init__(self, G1, G2, H1, H2, hnew):
        assert G1[2].is_one() and G2[2].is_one()
        self.G1 = G1
        self.G2 = G2
        self.H1 = H1
        self.H11 = H1*H1
        self.H12 = H1*H2
        self.H22 = H2*H2
        self.hnew = hnew
        self.x = hnew.parent().gen()

    def map(self, D):
        U, V = D
        assert U[2].is_one()
        # Sum and product of (xa, xb)
        s, p = -U[1], U[0]
        # Compute X coordinates (non reduced, degree 4)
        g1red = self.G1 - U
        g2red = self.G2 - U
        assert g1red[2].is_zero() and g2red[2].is_zero()
        g11, g10 = g1red[1], g1red[0]
        g21, g20 = g2red[1], g2red[0]
        # see above
        Px = (g11*g11*p + g11*g10*s + g10*g10) * self.H11 \
           + (2*g11*g21*p + (g11*g20+g21*g10)*s + 2*g10*g20) * self.H12 \
           + (g21*g21*p + g21*g20*s + g20*g20) * self.H22

        # Compute Y coordinates (non reduced, degree 3)
        assert V[2].is_zero()
        v1, v0 = V[1], V[0]
        # coefficient of y^2 is V(xa)V(xb)
        Py2 = v1*v1*p + v1*v0*s + v0*v0
        # coefficient of y is h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
        # so we need to symmetrize:
        # V(xa) Gred1(xb) (x-xb)
        # = (v1*xa+v0)*(g11*xb+g10)*(x-xb)
        # = (v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)*x
        # - xb*(v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)
        # Symmetrizing xb^2 gives u1^2-2*u0
        Py1 = (2*v1*g11*p + v1*g10*s + v0*g11*s + 2*v0*g10)*self.x \
          - (v1*g11*s*p + 2*v1*g10*p + v0*g11*(s*s-2*p) + v0*g10*s)
        Py1 *= self.H1
        # coefficient of 1 is Gred1(xa) Gred1(xb) h1(x)^2 U(x)
        Py0 = self.H11 * U * (g11*g11*p + g11*g10*s + g10*g10)

        # Now reduce the divisor, and compute Cantor reduction.
        # Py2 * y^2 + Py1 * y + Py0 = 0
        # y = - (Py2 * hnew + Py0) / Py1
        _, Py1inv, _ = Py1.xgcd(Px)
        Py = (- Py1inv * (Py2 * self.hnew + Py0)) % Px
        assert Px.degree() == 4
        assert Py.degree() == 3

        Dx = ((self.hnew - Py ** 2) // Px).monic()
        Dy = (-Py) % Dx
        return (Dx, Dy)

def jacobian_double(h, u, v):
    """
    Computes the double of a jacobian point (u,v)
    given by Mumford coordinates.

    See SAGE cantor_composition() and cantor_reduction
    """
    assert u.degree() == 2
    assert u[2] == 1
    # Replace u by u^2
    # Compute h3 the inverse of 2*v modulo u
    # Replace v by (v + h3 * (h - v^2)) % u
    q, r = u.quo_rem(2*v)
    if r[0] == 0: # gcd(u, v) = v, very improbable
        a = q.monic()**2
        b = (v + (h - v^2) // v) % a
        return a, b
    else: # gcd(u, v) = 1
        h3 = 1 / (-r[0]) * q
        a = u*u
        b = (v + h3 * (h - v**2)) % a
        # Cantor reduction
        Dx = ((h - b**2) // a).monic()
        Dy = (-b) % Dx
        return Dx, Dy

def jacobian_iter_double(h, u, v, n):
    for _ in range(n):
        u, v = jacobian_double(h, u, v)
    return u, v

def FromJacToJac(h, D11, D12, D21, D22, a, powers=None):
    # power is an optional list of precomputed tuples
    # (l, 2^l D1, 2^l D2) where l < a are increasing
    R = h.parent()
    Fp2 = R.base()

    J = HyperellipticCurve(h).jacobian()
    D1 = J(D11, D12)
    D2 = J(D21, D22)

    next_powers = None
    if not powers:
        # Precompute some powers of D1, D2 to save computations later.
        # We are going to perform O(a^1.5) squarings instead of O(a^2)
        if a >= 16:
            gap = Integer(a).isqrt()
            doubles = [(0, D1, D2)]
            _D1, _D2 = D1, D2
            for i in range(a-1):
                _D1 = jacobian_double(h, _D1[0], _D1[1])
                _D2 = jacobian_double(h, _D2[0], _D2[1])
                doubles.append((i+1, _D1, _D2))
            _, (G1, _), (G2, _) = doubles[a-1]
            next_powers = [doubles[a-2*gap], doubles[a-gap]]
        else:
            G1, _ = jacobian_iter_double(h, D1[0], D1[1], a-1)
            G2, _ = jacobian_iter_double(h, D2[0], D2[1], a-1)
    else:
        (l, _D1, _D2) = powers[-1]
        if a >= 16:
            next_powers = powers if l < a-1 else powers[:-1]
        G1, _ = jacobian_iter_double(h, _D1[0], _D1[1], a-1-l)
        G2, _ = jacobian_iter_double(h, _D2[0], _D2[1], a-1-l)

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
    R = RichelotCorr(G1, G2, H1, H2, hnew)

    imD1 = R.map(D1)
    imD2 = R.map(D2)
    if next_powers:
        next_powers = [(l, R.map(_D1), R.map(_D2))
            for l, _D1, _D2 in next_powers]
    return hnew, imD1[0], imD1[1], imD2[0], imD2[1], next_powers

def Does22ChainSplit(C, E, P_c, Q_c, P, Q, a):
    Fp2 = C.base()
    # gluing step
    h, D11, D12, D21, D22 = FromProdToJac(C, E, P_c, Q_c, P, Q, a);
    next_powers = None
    # print(f"order 2^{a-1} on hyp curve {h}")
    for i in range(1,a-2+1):
        h, D11, D12, D21, D22, next_powers = FromJacToJac(
            h, D11, D12, D21, D22, a-i, powers=next_powers)
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
        comp = EllipticCurveIsogeny(C, kerpol, degree=3, check=False)
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
