"""Check the new closed-disk proof's moment recurrence at a complex point."""

from math import prod

from sympy import I, QQ_I
from sympy.polys.rings import ring


def main():
    size = 6
    zero, one = QQ_I.zero, QQ_I.one
    t = QQ_I.from_sympy(7 * I / 6)  # |t| = 1+1/K
    polynomial_ring, x = ring("x", QQ_I)
    nodes = [t**a for a in range(1, size + 1)]
    denominator = prod((1 - x / r for r in nodes), start=polynomial_ring.one)
    derivative = denominator.diff(x)
    weights = [one / prod((one - r / s for s in nodes if r != s), start=one)
               for r in nodes]
    tails = []
    tail = zero
    for r in nodes:
        tail += r * (r + one) / (r - one)**3
        tails.append(tail)
    mu = [zero] + [j**2 / (t**j - one) for j in range(1, 2 * size)]

    # Two distinct choices determine the full affine dependence on Y.
    for y in (zero, QQ_I.from_sympy(1 + I)):
        simple_moments = [y - s for s in tails]
        quotient_used = []
        for k in range(2 * size - 1):
            by_recurrence = sum((w * u for w, u in zip(weights, simple_moments)), zero)
            quotient, remainder = divmod(x**k, denominator)
            polynomial_part = sum((c * mu[j] for (j,), c in quotient.items()), zero)
            if polynomial_part != zero:
                quotient_used.append(k)
            # The coefficient of 1/(1-x/r) is -remainder(r)/(r*D'(r)).
            proper_part = sum((-remainder(r) * (y - s) / (r * derivative(r))
                               for r, s in zip(nodes, tails)), zero)
            assert by_recurrence == polynomial_part + proper_part, (y, k)
            simple_moments = [r * (u - mu[k]) for r, u in zip(nodes, simple_moments)]
        assert quotient_used == [7, 8, 9, 10]
    print("All 11 affine moments agree over Q(i) at t=7i/6, K=6. "
          "The polynomial quotient contributes at k=7,8,9,10.")


if __name__ == "__main__":
    main()
