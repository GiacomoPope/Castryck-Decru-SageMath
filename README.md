# Castryck-Decru Key Recovery Attack on SIDH

SageMath implementation of [An efficient key recovery attack on SIDH, with Thomas Decru, preliminary version](https://eprint.iacr.org/2022/975.pdf), based on supplied Magma code from [https://homes.esat.kuleuven.be/~wcastryc/](https://homes.esat.kuleuven.be/~wcastryc/).

**Sage version**: This was written using SageMath 9.5, and works on the latest stable version: 9.6. I have been told that it doesn't work for 9.2 and below.

## Baby Example

During development of the code, we created a weaker parameter set `SIKEp64` with $p = 2^{33}\*3^{19} - 1$. This has the benefit of helping debug our implementation while also giving instant gratification of seeing an attack in real time.

Running `sage baby_SIDH.sage` on a laptop recovers Bob's private key in less than one minute.

## Breaking SIDH on a Laptop

|                          | `SIKEp64`  | `$IKEp217` | `SIKEp434` | `SIKEp503` | `SIKEp610` | `SIKEp751`   |
|--------------------------|------------|------------|------------|:----------:|------------|:------------:|
| Approximate Running Time | 30 seconds | 5 minutes  | 30 minutes | 45 minutes | 1.5 hours  | 2.5-9 hours  |

**Note**: Especially for the higher NIST levels, a lot of time is spent getting the first digits, and so performance time varies based on whether or not the first few values are `0` (fastest) or `2` (slowest). For example, attacking `SIKEp751`, similar hardware has been run multiple times with a compute times ranging from 2.5 hours to 9 hours. 

### Parameter choice

* To run the attack on the baby parameters, run `sage baby_SIDH.sage`
* To run the attack on the Microsoft `$IKEp217` challenge, run `sage SIKE_challenge.sage`
* To run the attack on the parameters submitted to the NIST PQ competition:
    * Default: `NIST_submission = "SIKEp434"`. Simply run `sage SIKEp434.sage` for an attack on `SIKEp434`.
    * Modify line 12: `NIST_submission = "SIKEp503"` for an attack against `SIKEp503`
    * Modify line 12: `NIST_submission = "SIKEp610"` for an attack against `SIKEp610`
    * Modify line 12: `NIST_submission = "SIKEp751"` for an attack against `SIKEp751`

## Estimating the running time

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

|             | $c$   | $\beta_1$ | Cost       |
|-------------|-------|-----------|------------|
| `SIKEp64`   | 1s    | 2         | 32 seconds |
| `$IKEp217`  | 4.5s  | 2         | 8 minutes  |
| `SIKEp434`  | 12s   | 2         | 45 minutes |
| `SIKEp503`  | 13s   | 4         | 1 hour     |
| `SIKEp610`  | 19s   | 5         | 2 hours    |
| `SIKEp751`  | 26s   | 6         | 5.5 hours  |

Where $c$ has been estimated using a MacBook Pro using a 6-Core Intel Core i7 CPU @ 2.6 GHz.


## Deviation from Castryck-Decru Attack

Roughly: points are now directly computed rather than derived by solving equations. This means we can avoid the very slow Grobner basis computation which Sage uses.

Deviation can be see in the file `richelot_aux.sage` in the functions: 

* `FromProdToJac()`
* `FromJacToJac()`

Thanks to [Rémy Oudompheng](https://twitter.com/oudomphe) for deriving and implementing these algorithms.

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

### Performance estimates with different patches

|                       | Vanilla :icecream: | No Proof :sleeping: | Monkey Patch :monkey_face: | Sage Patch 🩹 |
|-----------------------|:------------------:|:-------------------:|:--------------------------:|:-------------:|
| Baby SIDH (`SIKEp64`) | 1 minute           | 1 minute            | 1 minute                   | 1 minute      |
| `$IKEp217` Challenge  |          -         | 30 minutes          | 15 minutes                 | 15 minutes    |
| `SIKEp434`            |          -         |          -          | 1.5 hours                  | 1.5 hours     |
| `SIKEp503`            |          -         |          -          | 3.5 hours                  |       -       |
| `SIKEp610`            |          -         |          -          |              -             |       -       |
| `SIKEp751`            |          -         |          -          |              -             |       -       |


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
