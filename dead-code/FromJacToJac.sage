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