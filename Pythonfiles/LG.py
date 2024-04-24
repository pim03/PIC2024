import numpy as np
from scipy.special import factorial, genlaguerre

def LG(R, Phi, P, L, weights, w0):
    """
    This function computes a superposition of LG mode at the plane z=0
    R, Phi are coordinate matrices
    "weights" is a weight vector for the coefficients in the superposition
    """
    U = np.zeros(R.shape)  # initialise field
    for i in range(len(weights)):
        coeff = (np.sqrt(2 * factorial(P[i]) / (np.pi * factorial(P[i] + abs(L[i])))) *
                 (1 / w0) * (np.sqrt(2) * R / w0) ** abs(L[i]) *
                 np.exp(-R ** 2 / w0 ** 2) *
                 genlaguerre(P[i], abs(L[i]))(2 * R ** 2 / w0 ** 2) *
                 np.exp(1j * L[i] * Phi))
        U += weights[i] * coeff
    return U

# Example usage:
# x = np.linspace(-1, 1, 200)
# y = np.linspace(-1, 1, 200)
# X, Y = np.meshgrid(x, y)
# Phi, R = np.arctan2(Y, X), np.sqrt(X**2 + Y**2)
# P = [0, 0, 1, 3]
# L = [-5, 3, 4, 4]
# weights = [5, 1j, -1, 2]
# w0 = 1
# U = LG(R, Phi, P, L, weights, w0)


