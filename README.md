# Lambda-coalescents May 21 2020
C++ code for sampling from Lambda-coalescents -  only Lambda-coalescents, including Beta-coalescents (Schweinsberg, 2003) and  the Dirac coalescent (Eldon and Wakeley 2006); also eventually time-changed coalescents (see Fabian Freund (2020))
## Assumes the GNU GPL version 3+
## Compilation
requires the Boost library, and Rcpp from CRAN
in R:
> library(Rcpp)
> Rcpp::sourceCpp("clambdabeta.cpp")
