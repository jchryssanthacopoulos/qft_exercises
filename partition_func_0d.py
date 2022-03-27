"""Compute partition function in zero dimensions exactly and compare to asymptotic series."""

import matplotlib.pyplot as plt
import numpy as np
from sympy import Symbol, integrate, exp, oo, factorial


hbar = 1
coupling = 0.1
m = 1


def asymptotic_term(n):
    # get nth term in asymptotic series
    term = float(factorial(4*n) / (factorial(4)**n * factorial(n) * 4**n * factorial(2*n)))
    term *= (-hbar * coupling / m**4)**n
    return np.sqrt(2*np.pi*hbar)/m * term


x = Symbol('x')

# action
S = 0.5 * m**2 * x**2 + coupling * x**4 / factorial(4)

# partition function
Z = integrate(exp(-S/hbar), (x, -oo, oo))
Z = float(Z)


nterms = 50
aterms = []
asums = []
asum = 0

min_error = 0
n_min_error = 0
n_2 = 0
is_next_term_greater = False

for n in range(nterms):
    aterm = asymptotic_term(n)
    aterms.append(aterm)

    asum += aterm
    asums.append(asum)

    if n == 0 or np.abs(Z - asum) < min_error:
        min_error = np.abs(Z - asum)
        n_min_error = n

    if not is_next_term_greater and n > 0:
        if np.abs(aterm) > np.abs(aterms[n - 1]):
            n_2 = n
            is_next_term_greater = True


print("Partition function is", Z)
print("Minimum error occurs at n =", n_min_error)
print("Order n where terms start to grow =", n_2)

# plot
plt.figure()
plt.plot(range(nterms), [np.log(np.abs(Z - a)) for a in asums], 'rx-', label="asymptoic series")
plt.title("$Z-Z_n$")
plt.xlabel("Perturbation order $n$")
plt.xlim([-1, nterms])
plt.show()
