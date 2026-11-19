"""Reproduce the height example and the first four cases of Lemma 4.1."""

from fractions import Fraction
import json

import mpmath as mp
from sympy import QQ, cancel, diff, symbols
from sympy.polys.fields import field

from lambert_checks import clearer, normalized_coefficients


def check_example():
    # The proof uses sqrt(3) > 12/7 and the classical lower bound pi > 31/10.
    assert 12**2 < 3 * 7**2
    assert 9**6 < 33 * 7**5
    assert Fraction(7, 6) < Fraction(352, 2889) * Fraction(31, 10)**2

    with mp.workdps(60):
        root = mp.sqrt(3)
        cutoff = 352 * mp.pi**2 / 2889
        plus = mp.log(6) / mp.log(3 + root)
        minus = mp.log(6) / mp.log(3 - root)
        return {
            "ratio_at_larger_root": mp.nstr(plus, 15),
            "ratio_at_smaller_root": mp.nstr(minus, 15),
            "cutoff": mp.nstr(cutoff, 15),
        }


def check_pole_moment():
    x, r = symbols("x r")
    pole = 1 / (1 - x / r)
    euler_twice = x * diff(x * diff(pole, x), x)
    y = x / r
    assert cancel(euler_twice - y * (1 + y) / (1 - y)**3) == 0


def check_leading_terms():
    base, t = field("t", QQ)
    results = []
    for size in range(1, 5):
        coefficients = normalized_coefficients(t, size)
        assert len(coefficients) == size + 1
        assert all(p.denom == base.ring.one for p in coefficients)
        cyclotomic = clearer(t, size)
        degree = cyclotomic.numer.degree()
        sign = (-1) ** (size * (size - 1) // 2)
        assert cyclotomic.numer.LC == 1
        assert coefficients[-1] == sign * cyclotomic
        assert all(p == 0 or p.numer.degree() < degree for p in coefficients[:-1])
        results.append({"K": size, "degree_in_t": degree, "leading_sign": sign})
    return results


if __name__ == "__main__":
    check_pole_moment()
    result = {
        "exact_pole_moment": "passed",
        "rational_inequalities_for_the_larger_root": "passed",
        "rounded_values": check_example(),
        "exact_determinant_coefficients": check_leading_terms(),
    }
    print(json.dumps(result, indent=2))
