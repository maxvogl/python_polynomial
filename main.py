import numpy as np
import matplotlib.pyplot as plt

from polynomial import Polynomial
from pylab import mpl


plt.style.use("ggplot")
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"] = 10


def myplot(x, p, n, title):
    fig, ax=plt.subplots()
    ax.plot(x, p, linewidth=2, alpha=0.25, label="Polynomial class")
    ax.plot(x, n, linewidth=1, linestyle="dashed", label="Numpy")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    coefficients = [1, 2, 3]
    p1 = Polynomial(*coefficients)
    n1 = np.poly1d(coefficients)

    X = np.linspace(-5, 5, 101, True)
    myplot(X, p1(X), n1(X), "Compare Polynomial Class with Numpy")

    coeff = [0, 1, 2, 3, 4]
    p2 = Polynomial(*coeff)
    n2 = np.poly1d(coeff)
    
    p_sum = p1 + p2
    print("Sum:", str(p_sum))
    p_diff = p1 - p2
    print("Subtract:", str(p_diff))

    p_prime = p2.derivative()
    n_prime = n2.deriv()
    print("First derivative:", str(p_prime))
    p_prime_prime = p_prime.derivative()
    n_prime_prime = n2.deriv(2)
    print("Second derivative:", str(p_prime_prime))

    fig, ax=plt.subplots(nrows=3, sharex=True)
    ax[0].plot(X, p2(X), linewidth=2, alpha=0.25, label="Polynomial class")
    ax[0].plot(X, n2(X), linewidth=1, linestyle="dashed", label="Numpy")
    ax[0].set_ylabel("y")
    ax[0].set_title("Polynomials")

    ax[1].plot(X, p_prime(X), linewidth=2, alpha=0.25, label="Polynomial class")
    ax[1].plot(X, n_prime(X), linewidth=1, linestyle="dashed", label="Numpy")
    ax[1].set_ylabel("y")
    ax[1].set_title("First derivative")

    ax[2].plot(X, p_prime_prime(X), linewidth=2, alpha=0.25, label="Polynomial class")
    ax[2].plot(X, n_prime_prime(X), linewidth=1, linestyle="dashed", label="Numpy")
    ax[2].set_xlabel("x")
    ax[2].set_ylabel("y")
    ax[2].set_title("Second derivative")
    ax[2].legend()
    
    plt.show()