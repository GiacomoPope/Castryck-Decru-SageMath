import time
from itertools import product

load('richelot_aux.sage')
load('uvtable.sage')
load('speedup.sage')

# ===================================
# =====  ATTACK  ====================
# ===================================

def possibly_parallel(num_cores):
    if num_cores == 1:
        def _wrap(fun):
            def _fun(args):
                for a in args:
                    yield ((a,), None), fun(a)
            return _fun
        return _wrap
    return parallel(num_cores)

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=1):
    tim = time.time()

    skB = [] # TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY

    # gathering the alpha_i, u, v from table

    expdata = [[0, 0, 0] for _ in range(b-3)]
    # for i in [1..b-3] do
    for i in range(1,b-2):
        # if IsOdd(b-i) then
        if (b-i)%2 == 1:
            index = (b - i + 1) // 2
            exp = uvtable[index-1][1]
            if exp <= a:
                u = uvtable[index-1][2]
                v = uvtable[index-1][3]
                expdata[i-1] = [exp, u, v]

    # gather digits until beta_1
    bet1 = 0
    while expdata[bet1][0] == 0:
        bet1 += 1
    bet1 += 1

    ai = expdata[bet1-1][0]
    u  = expdata[bet1-1][1]
    v  = expdata[bet1-1][2]

    print(f"Determination of first {bet1} ternary digits. We are working with 2^{ai}-torsion.")

    bi = b - bet1
    alp = a - ai

    @possibly_parallel(num_cores)
    def CheckGuess(first_digits):
        print(f"Testing digits: {first_digits}")

        # tauhatkernel = 3^bi*P3 + sum([3^(k-1)*j[k-1] for k in range(1,beta+1)])*3^bi*Q3
        scalar = sum(3^k*d for k,d in enumerate(first_digits))
        tauhatkernel = 3^bi * (P3 + scalar*Q3)

        tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)

        C, P_c, Q_c = AuxiliaryIsogeny(bet1, u, v, E_start, P2, Q2, tauhatkernel, two_i)

        return Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)

    guesses = [i.digits(3, padto=bet1) for i in [0..3^bet1-2]]

    for result in CheckGuess(guesses):
        ((first_digits,), _), is_split = result

        # if j eq <2 : k in [1..bet1]> or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai) then
        if is_split:
            print("Glue-and-split! These are most likely the secret digits.")
            skB += first_digits
            break

    else:
        print("All other guesses failed, so first digits must be all 2!")
        skB += [2]*bet1

    print(skB)

    # now compute longest prolongation of Bob's isogeny that may be needed
    length = 1
    max_length = 0
    # for i in [bet1+1..b-3] do
    for i in range(bet1 + 1, b - 2):
      if expdata[i-1][0] != 0:
        max_length = max(length, max_length)
        length = 0
      else:
        length += 1

    while True:
        K = 2^a*3^(b - max_length)*EB.random_point()
        if K.order() == 3^max_length:
            break

    while True:
        alternativeK = 2^a*3^(b - max_length)*EB.random_point()
        if K.weil_pairing(alternativeK, 3^max_length)^(3^(max_length - 1)) != 1:
            break

    _, EBprolong = Pushing3Chain(EB, K, max_length)

    # gather next digit and change prolongation if needed
    i = bet1 + 1
    bi = b - i

    print(f"Determination of the {i}th ternary digit. We are working with 2^{ai}-torsion.")
    prolong = 1
    print("Prolonging with 1 steps.")
    endPB = EBprolong[0](PB)
    endQB = EBprolong[0](QB)
    endEB = EBprolong[0].codomain()
    # Speeds things up in Sage
    endEB.set_order((p+1)^2, num_checks=0)

    positives = []

    @possibly_parallel(num_cores)
    def CheckGuess(j):
        print(f"Testing digit: {j}")

        # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
        scalar = sum(3^k*d for k,d in enumerate(skB)) + 3^(i-1)*j
        tauhatkernel = 3^bi * (P3 + scalar*Q3)

        C, P_c, Q_c = AuxiliaryIsogeny(i, u, v, E_start, P2, Q2, tauhatkernel, two_i)

        return Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai)

    for result in CheckGuess([0,1,2]):
        ((j,), _), is_split = result
        if is_split:
            print("Glue-and-split!")
            positives.append(j)
            # continue testing other digits unless we reached an
            # ambiguity already:
            # By Remark 4 of [Castryck-Decru], there a probability
            # that K cancels the tail of Bob's isogeny, creating false
            # positives, in this case, we have to switch to alternativeK.
            if len(positives) > 1:
                break

    if len(positives) == 1:
        print(f"Most likely good prolongation and the secret digit is: {positives[0]}")
        skB.append(positives[0])
        print(skB)
        next_i = i + 1
    else:
        print("All glue-and-splits, must be bad prolongation: changing it and redoing this digit.")
        _, EBprolong = Pushing3Chain(EB, alternativeK, max_length)
        next_i = i
        prolong = 0


    # now gather all remaining digits, except for last three (we close that gap by trial and error)

    # for i in [next_i..b-3] do
    for i in range(next_i, b-2):
        bi = b - i
        if expdata[i-1][0] != 0:
            ai = expdata[i-1][0]
            alp = a - ai
            u = expdata[i-1][1]
            v = expdata[i-1][2]
            prolong = 0
        else:
            prolong += 1

        print(f"Determination of the {i}th ternary digit. We are working with 2^{ai}-torsion.")
        print(f"Prolonging with {prolong} steps.")
        endPB = PB
        endQB = QB
        endEB = EB
        # for j in [1..prolong] do
        for j in range(1,prolong+1):
            endPB = EBprolong[j-1](endPB)
            endQB = EBprolong[j-1](endQB)
        if j == prolong:
            endEB = EBprolong[j-1].codomain()
            # Speeds things up in Sage
            endEB.set_order((p+1)^2, num_checks=0)

        @possibly_parallel(num_cores)
        def CheckGuess(j):
            print(f"Testing digit: {j}")

            # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
            scalar = sum(3^k*d for k,d in enumerate(skB)) + 3^(i-1)*j
            tauhatkernel = 3^bi * (P3 + scalar * Q3)

            C, P_c, Q_c = AuxiliaryIsogeny(i, u, v, E_start, P2, Q2, tauhatkernel, two_i)

            return Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai)

        for result in CheckGuess([0,1]):
            ((j,), _), is_split = result

            if is_split:
                print("Glue-and-split! This is most likely the secret digit.")
                skB.append(j)
                print(skB)
                break
        else:
            print("All other guesses failed, so the digit must be 2")
            skB.append(2)
            print(skB)

    key = sum([skB[i]*3^(i) for i in range(b-3)])

    # bridge last safety gap
    tim2 = time.time()

    @possibly_parallel(num_cores)
    def CheckGuess(i):
        bobskey = key + i*3^(b-3)
        bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
        return bobscurve.j_invariant() == EB.j_invariant()

    for result in CheckGuess([0..3^3]):
        ((i,), _), found = result
        if found:
            bobskey = key + i*3^(b-3)
            break

    print(f"Bridging last gap took: {time.time() - tim2}")

    if found:
        print(f"Bob's secret key revealed as: {bobskey}")
        print(f"In ternary, this is: {Integer(bobskey).digits(base=3)}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return bobskey
    else:
        print("Something went wrong.")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return 0
