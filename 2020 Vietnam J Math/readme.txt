READ ME

This repository contains some codes that are useful to reproduce the results in the paper

Scarabel F, Breda D, Diekmann O, Gyllenberg M, Vermiglio R (2020). Numerical Bifurcation Analysis of Physiologically Structured Population Models via Pseudospectral 
Approximation, Vietnam J Mathematics, https://doi.org/10.1007/s10013-020-00421-3

Each example consists of two files:
PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
To perform the Matcont continuation, the system definition file ``PS_example'' must be copied into the subfolder "systems" of the Matcont folder.

MC_example: script for the Matcont continuation of the system defined in "PS_example".

The codes are tested on Matlab 2019b and Matcont version matcont6p6
