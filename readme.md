# PS-method-PDE-structured-models
 Pseudospectral discretization of structured models as PDEs

This repository contains the MATLAB codes for the pseudospectral discretization of structured population models in PDE formulation, and for the numerical bifurcation analysis using the package MatCont for MATLAB.

The main reference papers is
[VJM2020] Scarabel F, Breda D, Diekmann O, Gyllenberg M, Vermiglio R (2020). Numerical Bifurcation Analysis of Physiologically Structured Population Models via Pseudospectral Approximation, Vietnam J Mathematics, https://doi.org/10.1007/s10013-020-00421-3

Each example consists of two files:
PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
MC_example: script for the Matcont continuation of the system defined in "PS_example".

The codes are released under the MIT license (see file LICENSE.txt for details).
