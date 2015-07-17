# hsqc
haskell ab-initio quantum chemistry

![Build Status](https://travis-ci.org/peterspackman/hsqc.svg)

Currently in development, intended as a playground/learning tool
for functional concepts in ab-initio computational chemistry.

Current features:
  - Hartree Fock (SCF)
  - Traditional exact integrals (only for s and p type orbitals)
  - STO-3G basis
  - Simple Schwarz Inequality

TODO:
  - Implement PRISM algorithm for integrals
  - Add basis proper basis sets (i.e. not just minimal basis)
  - Pretty printing of matrices
  - Optional debug flags
  - TESTS
