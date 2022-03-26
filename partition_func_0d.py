"""Compute partition function in zero dimensions exactly and compare to asymptotic series."""

from sympy import Symbol, integrate, exp, oo, factorial


hbar = 1
coupling = 0.1
m = 1


x = Symbol('x')

# action
S = 0.5 * m**2 * x**2 + coupling * x**4 / factorial(4)

# partition function
Z = integrate(exp(-S/hbar), (x, -oo, oo))
Z = float(Z)


print("Partition function is", Z)
