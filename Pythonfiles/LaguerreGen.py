import numpy as np
from scipy.special import genlaguerre
from numpy.polynomial import Polynomial

def LaguerreGen(*args):
    if len(args) == 1:
        n = args[0]
        alpha = 0
    elif len(args) == 2:
        n = args[0]
        alpha = args[1]
    else:
        raise ValueError('n must be integer, and (optional) alpha >= -1')

    if not isinstance(n, int) or n < 0 or alpha < -1:
        raise ValueError('n must be integer, and (optional) alpha >= -1')

    L = genlaguerre(n, alpha)
    p = Polynomial(L)
    return p.coef[::-1]

# Example usage:
# y = LaguerreGen(3, 2)
# print(y)


