# Castryck-Decru Key Recovery Attack on SIDH

SageMath implementation of [An efficient key recovery attack on SIDH, preliminary version](https://eprint.iacr.org/2022/975.pdf), based on supplied Magma code from [https://homes.esat.kuleuven.be/~wcastryc/](https://homes.esat.kuleuven.be/~wcastryc/).

**Sage version**: This code was developed using SageMath 9.5. Certain isogeny functions require SageMath 9.5 and above, so if the code does not run, check your current version with `sage --version`.

## A Note on Reimplementing the Castryck-Decru Attack

A note [A Note on Reimplementing the Castryck-Decru Attack and Lessons Learned for SageMath](https://eprint.iacr.org/2022/1283) has been uploaded to eprint. We hope that this gives a good background for the work accomplished in this repo, as well as some more general tips for implementing cryptographic attacks in SageMath.

## Deviation from Castryck-Decru Attack

A recent commit introduces a modification of the Castryck-Decru attack in which only the first $\beta_1$ ternary digits must be guessed. Now, instead of recovering the remaining digits one by one, the secret isogeny is directly calculated from the result of the (2,2)-isogeny chain. 

This modification was introduced in [#12 Implement direct computation of isogeny once the first splitting is found](https://github.com/jack4818/Castryck-Decru-SageMath/pull/19) and was acomplished by [Rémy Oudompheng](https://twitter.com/oudomphe).

A description of the attack is described in: [A note on implementing direct isogeny determination in the Castryck-Decru SIKE attack](https://www.normalesup.org/~oudomphe/textes/202208-castryck-decru-shortcut.pdf).

Additional derivations appear in how the Glue-and-Split computations can be found, which can be seen in the file `richelot_aux.sage` in the functions: 

* `FromProdToJac()`
* `FromJacToJac()`

Thanks again to [Rémy Oudompheng](https://twitter.com/oudomphe) for deriving and implementing these faster algorithms.

## Baby Example

During development of the code, we created a weaker parameter set `SIKEp64` with $p = 2^{33}\cdot 3^{19} - 1$. This has the benefit of helping debug our implementation while also giving instant gratification of seeing an attack in real time.

Running `sage baby_SIDH.sage` on a laptop recovers Bob's private key in less than ten seconds.

## Breaking SIDH on a Laptop

| ~ Running Time                  | `SIKEp64`  | `$IKEp217` | `SIKEp434` | `SIKEp503`  | `SIKEp610`   | `SIKEp751`                     |
|---------------------------------|:----------:|------------|------------|-------------|--------------|--------------------------------|
| Paper Implementation (Magma)    |   -        | 6 minutes  | 62 minutes | 2h19m       | 8h15m        | 20h37m                         |
| Our implementation (SageMath)   | 5 seconds  | 2 minutes  | 10 minutes | 15 minutes  | 25 minutes   | 1-2 hours                      |
| Direct Computation (Oudompheng) | 2 seconds  | 9 seconds  | 22 seconds | 2 minutes   | 15 minutes   | 1 hour                         |
 
**Note**: Especially for the higher NIST levels, a lot of time is spent getting the first digits, and so performance time varies based on whether or not the first few values are `0` (fastest) or `2` (slowest). 

### Parameter choice

* To run the attack on the baby parameters, run `sage baby_SIDH.sage`
* To run the attack on the Microsoft `$IKEp217` challenge, run `sage SIKE_challenge.sage`
* To run the attack on the parameters submitted to the NIST PQ competition:
    * Default: `NIST_submission = "SIKEp434"`. Simply run `sage SIKEp434.sage` for an attack on `SIKEp434`.
    * Modify line 12: `NIST_submission = "SIKEp503"` for an attack against `SIKEp503`
    * Modify line 12: `NIST_submission = "SIKEp610"` for an attack against `SIKEp610`
    * Modify line 12: `NIST_submission = "SIKEp751"` for an attack against `SIKEp751`

### Parallelism

You can now run the attack partly in parallel thanks to work by Lorenz Panny.

* To run the attack on all available cores, simply add `--parallel` to the command line. Example: `sage SIKEp434.sage --parallel`.
* To choose the number of cores to use, manually change the value of `num_cores` in any of the attack scripts.
  * If someone wants to make a nice CLI please feel free to make a pull request.

Essentially, we can guess the first $\beta_1$ digits in parallel (which has a dramatic improvement for higher level NIST parameters) and then guess both $j=0$ and $j=1$ in parallel rather than serially for the remaining digits. This means we expect an approximate 1.6× speedup for `SIKEp434` and even more for higher levels.

Note that this optimization improves latency at the expense of throughput: The overall amount of work is higher, but the attack finishes quicker. Parallelism more fine-grained than simply testing all guesses simultaneously will certainly improve this, but this seems to be much less trivial to implement in SageMath.

## Estimating the running time (Castryck-Decru Attack)

We can estimate an average running time from the expected number of calls to the oracle `Does22ChainSplit()`. 

* For the first $\beta_1$ digits, we have to run at most $3^{\beta_1}$ calls to `Does22ChainSplit()` and half this on average
* For the remaining $b - \beta_1$ digits we will call `Does22ChainSplit()` once when $j = 0$ and twice when $j = 1,2$. We can then expect on average to call the oracle once a third of the time and twice for the rest
* Expressing the cost of a single call as $c$ we can estimate the total cost as:

$$
\textsf{Cost} = c \cdot \frac{3^{\beta_1}}{2} + c \cdot \frac{(b - \beta_1)}{3} + 2c \cdot \frac{2(b - \beta_1)}{3}
$$

$$
\textsf{Cost} = c \left(\frac{3^{\beta_1}}{2} + \frac{5(b - \beta_1)}{3} \right)
$$

|             | $c$   | $\beta_1$ | Cost        |
|-------------|-------|-----------|-------------|
| `SIKEp64`   | 0.2s  | 2         | 6.5 seconds |
| `$IKEp217`  | 1s    | 2         | 2 minutes   |
| `SIKEp434`  | 3.4s  | 2         | 13 minutes  |
| `SIKEp503`  | 4.5s  | 4         | 22 minutes  |
| `SIKEp610`  | 6s    | 5         | 43 minutes  |
| `SIKEp751`  | 8.4s  | 6         | 1.75 hours  |

Where $c$ has been estimated using a MacBook Pro using a Intel Core i7 CPU @ 2.6 GHz. **Note** as $c$ was benchmarked for the *first* oracle calls, these are over-estimates, as the oracle calls are faster as more digits are collected.

## Estimating the running time (Oudompheng Modification)

```
Coming soon...
```

## Speeding SageMath up using a cache

There is a SageMath performance issue with the group law for the Jacobian of a hyperelliptic curve. When testing equality, the code invokes `GF(p^k)(...)` for all coefficients. The constructor of the Finite Field includes a primality test for every call, which for larger primes is incredibly expensive.

Rémy managed to avoid this by patching SageMath itself, modifying `sage.categories.fields` so that the vector space is cached:

```py
from sage.misc.cachefunc import cached_method

        @cached_method
        def vector_space(self, *args, **kwds):
            ...
```

A gentler fix is to use `proof.arithmetic(False)`. This still requires constructing the same vector space again and again, but drops the expensive primality test on every call.

However, the easiest fix for fast performance is thanks to [Robin Jadoul](https://ur4ndom.dev). He found that we can achieve a similar result to Rémy's Sage patch with the following in-line monkey patch to our finite field by including the line

```py
Fp2.<i> = GF(p^2, modulus=x^2+1)
type(Fp2).vector_space = sage.misc.cachefunc.cached_method(type(Fp2).vector_space)
```

This speed up is included by default through loading in the file `speedup.sage` for each of the attack files.

Included below are some recorded times for running the scripts with and without various patches, before the `JacToJac()` optimisations which were implemented in pull requests #6-#9. 

### Additional Monkey patch for fixing the dimension

A slow call to compute the dimension is made when running `HyperElliptic(h).jacobian()` which always returns `1`. Caching should improve performance up to 20% (not yet fully benchmarked).

```
# No use calculating the dimension of HyperElliptic every single time
from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
AlgebraicScheme_subscheme_projective.dimension = lambda self: 1
```

## Conversion Progress

### Magma files to convert

- [x] Convert `uvtable.m`
- [x] Convert `SIKE_challenge.m`
- [x] Convert `SIKEp434.m`
- [x] Convert `richelot_aux.m`

### Functions to rewrite in `richelot_aux.m`:

- [x] `FromProdToJac()` (Deviation from original code) 
- [x] `FromJacToJac()` (Deviation from original code)
- [x] `Does22ChainSplit()`
- [x] `OddCyclicSumOfSquares()` (Code used to generate `uvtable`, not necessary to reimplement)
- [x] `Pushing3Chain()`
- [x] `Pushing9Chain()` (Obselete code, not implemented)

### Main attack to re-write:

- [x] Set small primes and torsion points
- [x] Fill `expdata` data table from `uvtable.sage`
- [x] Compute the first bits using `Glue-and-split`
- [x] Compute longest prolongation of Bob's isogeny
- [x] Compute next digit
- [x] Compute all digits except last 3
- [x] Find last three digits
