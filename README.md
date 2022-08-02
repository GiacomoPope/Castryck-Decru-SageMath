# Castryck-Decru Key Recovery Attack on SIDH

SageMath implementation of [An efficient key recovery attack on SIDH, with Thomas Decru, preliminary version](https://eprint.iacr.org/2022/975.pdf), based on supplied Magma code from [https://homes.esat.kuleuven.be/~wcastryc/](https://homes.esat.kuleuven.be/~wcastryc/).

## Deviation from Castryck-Decru Attack

**TODO: Further explaination**

Roughly: points are now directly computed rather than derived by solving equations. This means we can avoid the very slow Grobner basis computation which Sage uses.

Deviation can be see in the file `richelot_aux.sage` in the functions: 

* `FromProdToJac()`
* `FromJacToJac()`

Thanks to [Rémy Oudompheng](https://twitter.com/oudomphe) for deriving and implementing these algorithms.

## Baby Example

SageMath is significantly slower than Magma for Hyperelliptic computations, and so running the attack on SIKE parameters is still slow. To show the attack working we include `baby_SIDH.sage` which has a smaller prime $p = 2^{33}\*3^{19} - 1$. 

Running `baby_SIDH.sage` on a laptop recovers Bob's private key in less than one minute.

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
- [x] `OddCyclicSumOfSquares()` (Dead code, not implemented)
- [x] `Pushing3Chain()`
- [x] `Pushing9Chain()`

### Main attack to re-write:

- [x] Set small primes and torsion points
- [x] Fill `expdata` data table from `uvtable.sage`
- [x] Compute the first bits using `Glue-and-split`
- [x] Compute longest prolongation of Bob's isogeny
- [x] Compute next digit (Partial progress)
- [x] Compute all digits except last 3
- [x] Find last three digits
