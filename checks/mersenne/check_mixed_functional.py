"""Exact check of the finite-node version of the mixed Lambert functional.

Run with Python 3. The comparison uses ordinary derivatives on one side and
shifted finite tails plus a partial-fraction recurrence on the other.
"""

from fractions import Fraction
import json
from math import prod


def basis_product(q, first_degree, second_degree):
    """Multiply the two monic triple-root basis polynomials in ascending order."""
    coefficients = [Fraction(1)]
    for degree in (first_degree, second_degree):
        for index in range(1, degree + 1):
            root = q ** ((index + 2) // 3)
            expanded = [Fraction(0)] * (len(coefficients) + 1)
            for power, coefficient in enumerate(coefficients):
                expanded[power] -= root * coefficient
                expanded[power + 1] += coefficient
            coefficients = expanded
    return coefficients


def values_at_node(coefficients, x, inverse_poles):
    """Return f, theta f, theta^2 f for f=P/prod(1-x/r_a)."""
    value = first = second = Fraction(0)
    for coefficient in reversed(coefficients):
        second = second * x + 2 * first
        first = first * x + value
        value = value * x + coefficient

    weight = 1 / prod(1 - c * x for c in inverse_poles)
    slopes = [c / (1 - c * x) for c in inverse_poles]
    weight_first = weight * sum(slopes)
    weight_second = weight * (sum(slopes) ** 2 + sum(s * s for s in slopes))
    derivative = weight * first + weight_first * value
    second_derivative = weight * second + 2 * weight_first * first + weight_second * value
    return weight * value, x * derivative, x * derivative + x * x * second_derivative


def finite_functional(coefficients, q, size, nodes, mixture, early_count=0):
    a, b, c = mixture
    inverse_poles = [q**j for j in range(1, size + 1)]
    result = Fraction(0)
    for m in range(1, nodes + 1):
        value, theta, theta_squared = values_at_node(coefficients, q**m, inverse_poles)
        summand = a * (value - coefficients[0]) + b * theta + c * theta_squared
        if m <= early_count:
            assert (value, theta, theta_squared) == (0, 0, 0)
            assert summand == -a * coefficients[0]
        result += summand
    return result


def finite_partial_fraction_moments(base, size, nodes, mixture):
    a, b, c = mixture
    q = Fraction(1, base)
    degrees = 2 * size - 1
    polynomial_moments = [Fraction(0)]
    for j in range(1, degrees):
        geometric_sum = (1 - q ** (j * nodes)) / (base**j - 1)
        polynomial_moments.append((a + b * j + c * j**2) * geometric_sum)

    finite_tails = [Fraction(0)]
    for m in range(1, nodes + size + 1):
        y = Fraction(base**m)
        h = a / (y - 1) + b * y / (y - 1)**2 + c * y * (y + 1) / (y - 1)**3
        if mixture == (2, -3, 1):
            assert h == 2 / (y - 1)**3
        finite_tails.append(finite_tails[-1] + h)

    moments = [Fraction(0)] * degrees
    roots = [Fraction(base**a_index) for a_index in range(1, size + 1)]
    for a_index, root in enumerate(roots, start=1):
        weight = 1 / prod((1 - root / other for other in roots if other != root),
                          start=Fraction(1))
        shifted = finite_tails[a_index + nodes] - finite_tails[a_index]
        for j in range(degrees):
            moments[j] += weight * shifted
            shifted = root * (shifted - polynomial_moments[j])
    return moments


def main():
    monomial_count = basis_count = early_count = 0
    cases = [(2, 1), (2, 3), (2, 5), (3, 5)]
    mixtures = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (2, -3, 1)]
    nodes = 6
    for base, size in cases:
        q = Fraction(1, base)
        for mixture in mixtures:
            expected = finite_partial_fraction_moments(base, size, nodes, mixture)
            for j, value in enumerate(expected):
                monomial = [Fraction(0)] * j + [Fraction(1)]
                assert finite_functional(monomial, q, size, nodes, mixture) == value
                monomial_count += 1
            for i in range(size):
                for j in range(i, size):
                    coefficients = basis_product(q, i, j)
                    early = max(0, (max(i, j) + 2) // 3 - 1)
                    from_nodes = finite_functional(coefficients, q, size, nodes, mixture, early)
                    assert from_nodes == sum(a * b for a, b in zip(coefficients, expected))
                    basis_count += 1
                    early_count += early
    print(json.dumps({
        "status": "pass",
        "base_and_size": cases,
        "mixtures": mixtures,
        "nodes_per_functional": nodes,
        "exact_monomial_comparisons": monomial_count,
        "exact_triple_root_basis_comparisons": basis_count,
        "early_node_subtractions_checked": early_count,
    }, indent=2))


if __name__ == "__main__":
    main()
