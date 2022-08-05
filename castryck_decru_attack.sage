import time
from itertools import product

load('richelot_aux.sage')
load('uvtable.sage')
load('speedup.sage')

# Remove annoying messages about slow Gröbner
set_verbose(-1)

# ===================================
# =====  ATTACK  ====================
# ===================================

def check_progress(solution, skB, tim):
    """
    There is currently a bug in which sometimes the
    wrong digit is calculated. When this happens, the
    whole algorithm runs till the end before failing

    We include a cheat-test on each digit so we can exit 
    earlier whenever we know the secret (such as in 
    baby SIDH or in the SIKEp434 file)
    """
    if skB != solution[:len(skB)]:
        print(f"The last digited calculated is incorrect.")
        print(f"Expected: {solution[:len(skB)]}")
        print(f"Computed: {skB}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        exit()

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, solution=None):
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

    # for j in CartesianPower([0,1,2], bet1) do
    for first_digits in product([0,1,2], repeat=int(bet1)):
        print(f"Testing digits: {[first_digits[k] for k in range(bet1)]}")

        # tauhatkernel = 3^bi*P3 + sum([3^(k-1)*j[k-1] for k in range(1,beta+1)])*3^bi*Q3
        tauhatkernel = 3^bi*P3 
        for k in range(1, bet1+1):
            tauhatkernel += (3^(k-1)*first_digits[k-1])*3^bi*Q3

        tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)

        C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, bet1)

        P_c = u*P2 + v*two_i(P2) 
        for taut in tau_tilde:
            P_c = taut(P_c)
        Q_c = u*Q2 + v*two_i(Q2)
        for taut in tau_tilde:
            Q_c = taut(Q_c)

        # if j eq <2 : k in [1..bet1]> or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai) then
        if first_digits == (2,)*bet1 or Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai):
            print("Glue-and-split! These are most likely the secret digits.")
            skB += first_digits
            break

    print(skB)
    if solution: check_progress(solution, skB, tim)

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
    endEB.set_order((p+1)^2)

    positives = []

    for j in range(0,3):
        print(f"Testing digit: {j}")
        # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
        tauhatkernel = 3^bi*P3 
        for k in range(1, i):
            tauhatkernel += (3^(k-1)*skB[k-1])*3^bi*Q3
        tauhatkernel += 3^(i-1)*j*3^bi*Q3

        tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)

        C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, i)
        P_c = u*P2 + v*two_i(P2)
        for taut in tau_tilde:
            P_c = taut(P_c)
        Q_c = u*Q2 + v*two_i(Q2)
        for taut in tau_tilde:
            Q_c = taut(Q_c)

        if Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai):
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
        if solution: check_progress(solution, skB, tim)
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
            endEB.set_order((p+1)^2)


        for j in range(0,3):
            print(f"Testing digit: {j}")
            # tauhatkernel := 3^bi*P3 + (&+[3^(k-1)*skB[k] : k in [1..i-1]] + 3^(i-1)*j)*3^bi*Q3;
            tauhatkernel = 3^bi*P3 
            for k in range(1, i):
                tauhatkernel += (3^(k-1)*skB[k-1])*3^bi*Q3
            tauhatkernel += 3^(i-1)*j*3^bi*Q3

            tauhatkernel_distort =u*tauhatkernel + v*two_i(tauhatkernel)
            
            C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, i)
            P_c = u*P2 + v*two_i(P2)
            for taut in tau_tilde:
                P_c = taut(P_c)
            Q_c = u*Q2 + v*two_i(Q2)
            for taut in tau_tilde:
                Q_c = taut(Q_c)
            if j == 2 or Does22ChainSplit(C, endEB, 2^alp*P_c, 2^alp*Q_c, 2^alp*endPB, 2^alp*endQB, ai):
                print("Glue-and-split! This is most likely the secret digit.")
                skB.append(j)
                if solution: check_progress(solution, skB, tim)
                print(skB)
                break

    key = sum([skB[i]*3^(i) for i in range(b-3)])

    # bridge last safety gap
    tim2 = time.time()

    found = false
    for i in range(3^3+1):
        bobskey = key + i*3^(b-3)
        bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
        if bobscurve.j_invariant() == EB.j_invariant():
            found = true
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








