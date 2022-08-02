# Groebner Problems

In the file `richelot_aux.m` in the function `FromProdToJac()` a system of five equations in the multivariate Polynomial Ring $F[U_0, U_1, V_0, U_1]$ over the field $K = F_{p^2}$ are solved using the following Magma code:

```c
A4<U0, U1, V0, V1> := AffineSpace(Fp2, 4);
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
```

A similar piece of code using SageMath would be:

```py
A4.<U0, U1, V0, V1> = AffineSpace(Fp2, 4)
V = A4.subscheme([eq1, eq2, eq3, eq4, eq5])

# point with zero coordinates probably correspond to "extra" solutions, we should be left with 4 sols
# (code may fail over small fields)

realsols = []
for D in V.rational_points():
    print(D)
    Dseq = list(D)
    if not 0 in Dseq:
        realsols.append(Dseq)
```

However, running this code we fall back to the INCREDIBLY slow Groebner basis implementation
```py
verbose 0 (3848: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
verbose 0 (1081: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
```

The current goal is to re-write the above code in a form that either `singular` or `Macaulay2` interfaces to play nicely with our code, but both seem to be really upset at trying to do this in the extension field `GF(p^2, modulus=x^2+1)` and give error messages.

One suggestion has to be Weil-restrict the system to move form $F_{p^2}$ to $F_p$ however, this increases the complexity of what we need to solve.

Another suggestion ahs been to instead use resultants, but I get errors unless I compute resultants with

```py
from sage.matrix.matrix2 import Matrix 
def resultant(f1, f2, var):
    return Matrix.determinant(f1.sylvester_matrix(f2, var))
```

Which is too slow for our current needs.