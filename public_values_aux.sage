def generate_automorphism(E):
    phi = EllipticCurveIsogeny(E, x)
    E1728 = phi.codomain()
    # Speeds things up in Sage
    E1728.set_order((p+1)^2)

    for iota in E1728.automorphisms():
        P = E1728.random_point()
        if iota(iota(P)) == -P:
            two_i = -phi.dual()*iota*phi
            return two_i

def generate_torsion_points(E, a, b):
    def get_l_torsion_basis(E, l):
        n = (p+1) // l
        return (n*G for G in E.gens())

    P2, Q2 = get_l_torsion_basis(E, 2^a)
    P3, Q3 = get_l_torsion_basis(E, 3^b)

    return P2, Q2, P3, Q3

def check_torsion_points(E, a, b, P2, Q2, P3, Q3):
    # Make sure Torsion points are
    # generated correctly
    infty = E(0)
    assert 2^(a-1)*P2 != infty
    assert 3^(b-1)*P3 != infty
    assert P2.weil_pairing(Q2, 2^a)^(2^(a-1)) != 1
    assert P3.weil_pairing(Q3, 3^b)^(3^(b-1)) != 1

def gen_bob_keypair(E_start, P3, Q3):
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
    return Bobskey, EB, PB, QB