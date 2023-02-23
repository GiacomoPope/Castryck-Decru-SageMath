from richelot_aux import *
load('speedup.sage')

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

test_FromProdToJac()


def test_FromJacToJac():
    print("test_FromJacToJac")
    h, D11, D12, D21, D22, _ = test_FromProdToJac()
    for e in FromJacToJac(h, D11, D12, D21, D22, 60):
        print(e)

test_FromJacToJac()

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

test_ChainSplit()
