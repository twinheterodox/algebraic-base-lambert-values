"""Exact mixed determinants for K <= 3 and the first large-order quotient of positive rank."""

from math import comb

import sympy as sp


t, x, A, B, C, X = sp.symbols("t x A B C X")
variables = (A, B, C, X)
polynomial_variables = (t, *variables)


def rational_moments(size):
    """Use ordinary polynomial division and the simple-fraction values at each pole."""
    roots = [t**a for a in range(1, size + 1)]
    base_field = sp.QQ.frac_field(t)
    denominator = sp.Poly(sp.prod(1 - x / r for r in roots), x, domain=base_field)

    tails = [sp.S.Zero] * 3
    poles = []
    for root in roots:
        tails[0] += 1 / (root - 1)
        tails[1] += root / (root - 1)**2
        tails[2] += root * (root + 1) / (root - 1)**3
        other_factors = sp.prod(1 - root / s for s in roots if s != root)
        poles.append((root, other_factors, (-tails[0], -tails[1], -tails[2], 1)))

    moments = []
    for k in range(2 * size - 1):
        quotient, _ = sp.Poly(x**k, x, domain=base_field).div(denominator)
        coefficients = [sp.S.Zero] * 4
        for (j,), coefficient in quotient.terms():
            if j > 0:
                for power in range(3):
                    coefficients[power] += coefficient * j**power / (t**j - 1)
        for root, other_factors, tail_moment in poles:
            residue = root**k / other_factors
            for index, coefficient in enumerate(tail_moment):
                coefficients[index] += residue * coefficient
        moments.append([sp.cancel(coefficient) for coefficient in coefficients])
    return moments


def clearer(size):
    factors = []
    for order in range(1, size + 1):
        if 2 * order <= size:
            numerator = 4 * size**2 + 32 * size * order + 17 * order**2
            exponent = (numerator + 16 * order - 1) // (16 * order)
        else:
            exponent = 5 * (size - order) + 3
        factors.append(sp.cyclotomic_poly(order, t)**exponent)
    return sp.prod(factors)


def check(size):
    moments = rational_moments(size)
    common = sp.S.One
    for moment in moments:
        for coefficient in moment:
            common = sp.lcm(common, sp.fraction(coefficient)[1])
    entries = [
        sum(s * sp.cancel(common * c) for s, c in zip(variables, moment))
        for moment in moments
    ]
    matrix = sp.Matrix(size, size, lambda i, j: entries[i + j])
    numerator = sp.Poly(matrix.det(method="domain-ge"), *polynomial_variables, domain=sp.ZZ)
    multiplier = sp.Poly(clearer(size), *polynomial_variables, domain=sp.ZZ)

    # Check the order at zero in the same normalization.
    loss = sum(comb(max(size - 2 * i, 0) + 1, 2) for i in range(size))
    exponent = size**2 * (size + 1) // 2 - loss
    divisor = sp.Poly(common**size * t**exponent, *polynomial_variables, domain=sp.ZZ)
    cleared = (multiplier * numerator).exquo(divisor)
    assert all(sum(powers[1:]) == size for powers in cleared.monoms())
    if size == 1:
        concrete = cleared.as_expr().subs({A: 2, B: -3, C: 1})
        assert sp.expand(concrete - ((t - 1)**3 * X - 2)) == 0
    return exponent, cleared.degree(t), len(cleared.terms())


def check_polynomial_part():
    # K=5 and d=3 is the first double-cluster case with a singular polynomial moment.
    size, order = 5, 3
    base_field = sp.QQ.frac_field(t)
    denominator = sp.Poly(sp.prod(1 - x / t**a for a in range(1, size + 1)), x,
                          domain=base_field)
    entries = []
    for k in range(2 * size - 1):
        quotient, _ = sp.Poly(x**k, x, domain=base_field).div(denominator)
        entries.append(quotient.nth(order))
    matrix = sp.Matrix(size, size, lambda i, j: entries[i + j])
    assert matrix == sp.diag(0, 0, 0, 0, -t**15)


if __name__ == "__main__":
    for size in (1, 2, 3):
        exponent, degree, terms = check(size)
        print(f"K={size}  e_K={exponent}  degree in t={degree}  nonzero monomials={terms}  PASS")
    check_polynomial_part()
    print("K=5, d=3  polynomial coefficient form has rank one  PASS")
