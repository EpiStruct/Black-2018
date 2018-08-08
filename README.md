# Importance sampling for partially observed temporal epidemic models

This is code to support the paper: Black, A. J. (2018) Importance sampling for partially observed temporal epidemic models. Statistics and Computing. [doi:10.1007/s11222-018-9827-1](https://doi.org/10.1007/s11222-018-9827-1).

If you make use of this code, please cite the paper. 

The first parts are coded in MATLAB and illustrate the basic importance sampling technique for a number of increasingly complex epidemic models (SIR, SEIR, SEEIIRp and SEIAR). The second part is a particle marginal Metropolis Hastings algorithm for implementing Bayesian inference using the importance sampling (SEIAR_pmMH). This is coded in C and a different version based on an alive filter is provided for comparison.

