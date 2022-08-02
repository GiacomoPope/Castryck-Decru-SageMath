# Castryck-Decru Key Recovery Attack on SIDH

Attempt at a SageMath implementation of [An efficient key recovery attack on SIDH, with Thomas Decru, preliminary version](https://eprint.iacr.org/2022/975.pdf), based on supplied Magma code from [https://homes.esat.kuleuven.be/~wcastryc/](https://homes.esat.kuleuven.be/~wcastryc/).

## Where we at?

Current goal is to re-implement `richelot_aux.m` and use this to break the key for a small $p = 2^{33}\*3^{19} - 1$. This can be found partially finished in `baby_SIDH.sage`. 

After this is done, it should be easy to make the changes for the harder challenges.

## Magma files to convert

- [x] Convert `uvtable.m`
- [ ] Convert `SIKE_challenge.m`
- [ ] Convert `SIKEp434.m`
- [ ] Convert `richelot_aux.m`

## Functions to rewrite in `richelot_aux.m`:

- [ ] `FromProdToJac()` (Partial Progress) 
- [ ] `FromJacToJac()`
- [x] `Does22ChainSplit()`
- [ ] `OddCyclicSumOfSquares()`
- [x] `Pushing3Chain()`
- [x] `Pushing9Chain()`

## Progress in Baby SIDH

- [x] Set small primes and torsion points
- [x] Fill `expdata` data table from `uvtable.sage`
- [ ] Compute the first bits using `Glue-and-split` (Partial progress)
- [ ] Compute longest prolongation of Bob's isogeny
- [ ] Compute next digit
- [ ] Compute all digits except last 3
- [ ] Find last three digits
