"""Independent small-size checks for the 18 November Lambert manuscript.

Run with Python 3, SymPy and mpmath. No code from the manuscript or its
integer-base predecessor is imported. Exact and numerical outputs are labeled.
"""

from fractions import Fraction
from math import factorial, prod

import mpmath as mp
from sympy import QQ
from sympy.polys.fields import field


def moment_coefficients(t, size):
    # 1/D(x) = sum_a w_a/(1-x/r_a). For this simple fraction,
    # L(x^(k+1)/(1-x/r_a)) = r_a*(L(x^k/(1-x/r_a)) - L(x^k)).
    one = t**0
    zero = one - one
    nodes = [t**a for a in range(1, size + 1)]
    count = 2 * size - 1
    moments = [zero] + [j**2 / (t**j - one) for j in range(1, count)]
    slope = [zero] * count
    offset = [zero] * count
    pole_tail = zero
    for a, node in enumerate(nodes):
        weight = one / prod((one - node / other for i, other in enumerate(nodes)
                             if i != a), start=one)
        pole_tail += node * (node + one) / (node - one)**3
        linear, constant = one, -pole_tail
        for k in range(count):
            slope[k] += weight * linear
            offset[k] += weight * constant
            linear = node * linear
            constant = node * (constant - moments[k])
    return slope, offset


def hankel_determinant(sequence, size):
    matrix = [[sequence[i + j] for j in range(size)] for i in range(size)]
    sign, last = 1, 1
    for j in range(size - 1):
        pivot_row = next((i for i in range(j, size) if matrix[i][j] != 0), None)
        if pivot_row is None:
            return sequence[0] * 0
        if j != pivot_row:
            matrix[j], matrix[pivot_row] = matrix[pivot_row], matrix[j]
            sign *= -1
        pivot = matrix[j][j]
        for i in range(j + 1, size):
            for k in range(j + 1, size):
                matrix[i][k] = (pivot * matrix[i][k] - matrix[i][j] * matrix[j][k]) / last
        last = pivot
    return sign * matrix[-1][-1]


def power_coefficients(values):
    # Newton's forward polynomial, converted to the ordinary power basis.
    differences = list(values)
    result = [values[0] * 0 for _ in values]
    basis = [1]
    for j in range(len(values)):
        scale = differences[0] / factorial(j)
        for k, coeff in enumerate(basis):
            result[k] += scale * coeff
        basis = [left - j * right for left, right in zip([0] + basis, basis + [0])]
        differences = [b - a for a, b in zip(differences, differences[1:])]
    return result


def exponent(size, d):
    if 2 * d > size:
        return 5 * (size - d) + 3
    numerator = 4 * size**2 + 32 * size * d + 17 * d**2
    return (numerator + 16 * d - 1) // (16 * d)


def clearer(t, size):
    polys = {}
    for d in range(1, size + 1):
        polys[d] = (t**d - 1) / prod(polys[j] for j in range(1, d) if d % j == 0)
    return prod(p**exponent(size, d) for d, p in polys.items())


def normalized_coefficients(t, size):
    slopes, offsets = moment_coefficients(t, size)
    c = clearer(t, size)
    prefactor = c / t**((size**3 - size) // 2)
    data = []
    for x in range(size + 1):
        entries = [a * x + b for a, b in zip(slopes, offsets)]
        data.append(prefactor * hankel_determinant(entries, size))
    coefficients = power_coefficients(data)
    assert coefficients[-1] == (-1)**(size * (size - 1) // 2) * c
    return coefficients


def exact_checks():
    field_t, t = field('t', QQ)
    for size in range(1, 5):
        coefficients = normalized_coefficients(t, size)
        assert all(p.denom == field_t.ring.one for p in coefficients)
        assert all(c.denominator == 1 for p in coefficients for _, c in p.numer.terms())
        if size == 1:
            assert coefficients == [-t * (t + 1), (t - 1)**3]
        print(f'Exact over Q(t), K={size}: integer coefficients in t,X; '
              f'degree in t = {max(p.numer.degree() for p in coefficients)}.', flush=True)
    counter = 0
    for base in (-3, -2, 2, 3, 4, 6):
        for size in range(1, 9):
            coefficients = normalized_coefficients(Fraction(base), size)
            assert all(p.denominator == 1 for p in coefficients)
            counter += 1
    print(f'Exact at integer bases: all {counter} full polynomials in X passed.', flush=True)


def relative_error(actual, expected):
    return abs(actual - expected) / max(1, abs(actual), abs(expected))


def expanding_sample(t, size, terms):
    # Evaluate the original Euler derivatives of x^k/D at the geometric
    # nodes. Use a monomial Lambert series for the independent value of X.
    q = 1 / t
    q_powers = [q**a for a in range(1, size + 1)]
    z_value = mp.fsum(j**2 * q**j / (1 - q**j) for j in range(1, terms + 1))
    numeric = [mp.mpc(0) for _ in range(2 * size - 1)]
    for m in range(1, terms + 1):
        x = q**m
        monomials = [p * x for p in q_powers]
        inv_d = 1 / mp.fprod(1 - u for u in monomials)
        log_derivative = sum(u / (1 - u) for u in monomials)
        derivative_again = sum(u / (1 - u)**2 for u in monomials)
        factor = inv_d
        for k in range(len(numeric)):
            numeric[k] += factor * ((k + log_derivative)**2 + derivative_again)
            factor *= x
    a, b = moment_coefficients(t, size)
    return max(relative_error(value, linear * z_value + constant)
               for value, linear, constant in zip(numeric, a, b))


def contracting_sample(t, y, size, terms):
    # Compare the convergent h-series formula for x^k/F with the separate
    # partial-fraction recurrence above, including the unique k=0 correction.
    n = terms + 2 * size
    hs = [mp.mpc(1)] + [mp.mpc(0) for _ in range(n)]
    for a in range(1, size + 1):
        node = t**a
        for j in range(1, n + 1):
            hs[j] += node * hs[j - 1]
    cs = [j**2 * t**j / (1 - t**j) for j in range(1, terms + 1)]
    u = y + mp.fsum(cs)
    a_size = size * (size + 1) // 2
    scale = (-1)**size * t**a_size
    slopes, offsets = moment_coefficients(t, size)
    errors = []
    for k in range(2 * size - 1):
        quotient = sum(j**2 / (t**j - 1) * hs[k - size - j]
                       for j in range(1, k - size + 1))
        series = mp.fsum(c * hs[k + j - size] for j, c in enumerate(cs, start=1)
                         if k + j >= size)
        constant = u / scale if k == 0 else (-u * hs[k - size] if k >= size else 0)
        comparison = scale * (constant + series + quotient)
        errors.append(relative_error(comparison, slopes[k] * y + offsets[k]))
    return max(errors)


if __name__ == '__main__':
    exact_checks()
    mp.mp.dps = 82
    phi = (1 + mp.sqrt(5)) / 2
    worst = mp.mpf(0)
    for t in (phi, mp.mpc(1, 1), mp.mpf(-2)):
        for size in range(1, 6):
            worst = max(worst, expanding_sample(t, size, 540))
    assert worst < mp.mpf('1e-60')
    print(f'Expanding moment comparison, K=1..5, three bases: {mp.nstr(worst, 5)} '
          'maximum scaled discrepancy.', flush=True)
    worst = mp.mpf(0)
    for t in (-1 / phi, mp.mpc(1, 2) / 5):
        for size in range(1, 7):
            worst = max(worst, contracting_sample(t, mp.mpc(2, -3), size, 430))
    assert worst < mp.mpf('1e-60')
    print(f'Contracting moment comparison, K=1..6, two bases: {mp.nstr(worst, 5)} '
          'maximum scaled discrepancy.', flush=True)
    gamma = 321 / (64 * mp.pi**2)
    print(f'gamma = {mp.nstr(gamma, 20)}; epsilon = '
          f'{mp.nstr(mp.mpf(11) / 18 - gamma, 20)}.')
