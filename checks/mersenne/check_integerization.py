"""Exact checks of the two moment normalizations and the cleared polynomial."""

from fractions import Fraction as F
from itertools import permutations
from math import ceil, prod
import json


def root_product(roots):
    poly = [1]
    for root in roots:
        new = [0] * (len(poly) + 1)
        for j, c in enumerate(poly):
            new[j] -= root * c
            new[j + 1] += c
        poly = new
    return poly


def monic_quotient(k, divisor):
    degree = len(divisor) - 1
    remainder = [0] * k + [1]
    quotient = [0] * max(0, k - degree + 1)
    for j in range(k, degree - 1, -1):
        c = remainder[j]
        quotient[j - degree] = c
        for v, coeff in enumerate(divisor):
            remainder[j - degree + v] -= c * coeff
    return quotient


def moments(base, size, mixture):
    a, b, c = mixture
    roots = [F(base**j) for j in range(1, size + 1)]
    monic = root_product(roots)
    scale = (-1)**size * prod(roots)
    max_k = 2 * size - 2
    mu = [0] + [F(a + b*j + c*j*j, base**j - 1) for j in range(1, max_k + 1)]
    tails = []
    tail = F(0)
    for r in roots:
        term = a/(r-1) + b*r/(r-1)**2 + c*r*(r+1)/(r-1)**3
        if mixture == (2, -3, 1):
            assert term == 2/(r-1)**3
        tail += term
        tails.append(tail)

    divided = []
    recurrence = [[F(0), F(0)] for _ in range(max_k + 1)]
    for k in range(max_k + 1):
        qpoly = monic_quotient(k, monic)
        const = sum(v * m for v, m in zip(qpoly, mu))
        slope = F(0)
        for r, tail in zip(roots, tails):
            coefficient = -r**(k-1) / prod(r-s for s in roots if s != r)
            const -= coefficient * tail
            slope += coefficient
        divided.append([scale * const, scale * slope])

    for r, tail in zip(roots, tails):
        weight = 1 / prod((1-r/s for s in roots if s != r), start=F(1))
        const, slope = -tail, F(1)
        for k in range(max_k + 1):
            recurrence[k][0] += weight * const
            recurrence[k][1] += weight * slope
            const, slope = r*(const-mu[k]), r*slope
    assert divided == recurrence
    return divided


def determinant_polynomial(size, sequence):
    coeffs = [F(0) for _ in range(size + 1)]
    for perm in permutations(range(size)):
        inversions = sum(v > w for j, v in enumerate(perm) for w in perm[j+1:])
        term = [F((-1)**inversions)]
        for i, j in enumerate(perm):
            const, slope = sequence[i+j]
            new = [F(0) for _ in range(len(term) + 1)]
            for k, val in enumerate(term):
                new[k] += const * val
                new[k+1] += slope * val
            term = new
        for j, val in enumerate(term):
            coeffs[j] += val
    return coeffs


def normalization(base, size):
    phi_values = {}
    clearer = 1
    for d in range(1, size + 1):
        previous = prod(v for h, v in phi_values.items() if d % h == 0)
        phi_values[d], rem = divmod(base**d - 1, previous)
        assert rem == 0
        if 2*d <= size:
            height = ceil(F(size*size, 4*d) + 2*size + F(17*d, 16))
        else:
            height = 5*(size-d) + 3
        clearer *= phi_values[d]**height
    loss = 0
    for i in range(size):
        h = max(0, size - 2*i)
        loss += h * (h+1) // 2
    exponent = size**2 * (size+1) // 2 - loss
    return F(clearer, base**exponent)


def main():
    mixtures = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (2, -3, 1)]
    cases = [(2, size, mix) for size in range(1, 5) for mix in mixtures]
    cases.append((6, 4, (2, -3, 1)))
    moment_comparisons = coefficient_checks = 0
    for base, size, mixture in cases:
        sequence = moments(base, size, mixture)
        moment_comparisons += len(sequence)
        coeffs = determinant_polynomial(size, sequence)
        leading = (-1)**(size*(size-1)//2) * base**(size*(size*size-1)//2)
        assert coeffs[-1] == leading
        factor = normalization(base, size)
        assert all((factor * c).denominator == 1 for c in coeffs)
        coefficient_checks += len(coeffs)
    assert F(19, 36) > F(321, 64) / F(31, 10)**2
    print(json.dumps({
        "status": "pass",
        "cases": len(cases),
        "base_2": {"sizes": [1, 2, 3, 4], "mixtures": mixtures},
        "composite_base": {"base": 6, "size": 4, "mixture": [2, -3, 1]},
        "affine_moment_comparisons": moment_comparisons,
        "cleared_polynomial_coefficients": coefficient_checks,
        "leading_coefficients": len(cases),
        "strict_rational_margin_check": True
    }, indent=2))


if __name__ == "__main__":
    main()
