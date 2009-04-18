# Test functions are from Tzitzler E. et al, Comparison of Multiobjective 
# Evolutionary Algorithms: Empirical Results, Evolutionary computation, 2000,
# vol 8(2), pp. 173-195.
#
# Arguments: 
# contenders - those whose fitness is to be evaluated

from numpy import array, sqrt, hstack, atleast_2d

def tau1(contenders):
    contenders = array(contenders, copy=False, ndmin=2)
    fit1 = contenders[:,0]
    g = 1 + 9./(contenders.shape[1] - 1) * contenders[:,1:].sum(axis=1)
    h = 1 - sqrt(fit1 / g)
    return hstack((fit1[:,None], (g*h)[:,None]))
