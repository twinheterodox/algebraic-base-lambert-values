"""Independent exact checks for the 18 November 2026 reciprocal Mersenne cubes PDF.

Run with Python 3 and SymPy. All inputs are below. The program reads no other files.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations, product
from math import comb, prod
import json

from sympy import QQ, QQ_I, Symbol, cyclotomic_poly
from sympy.polys.fields import field
from sympy.polys.matrices import DomainMatrix


def polymul(left, right):
    coefficients = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            coefficients[i + j] += a * b
    return coefficients


def evaluate(poly, x):
    value = 0
    for coefficient in reversed(poly):
        value = value * x + coefficient
    return value


def divrem(poly, denominator):
    remainder = list(poly)
    quotient = [0] * max(0, len(poly) - len(denominator) + 1)
    for index in reversed(range(len(quotient))):
        value = remainder[index + len(denominator) - 1] / denominator[-1]
        quotient[index] = value
        for j, coefficient in enumerate(denominator):
            remainder[index + j] -= value * coefficient
    assert not any(remainder[len(denominator) - 1:])
    return quotient, remainder[:len(denominator) - 1]


def prepare_poles(base, size):
    roots = [base**a for a in range(1, size + 1)]
    denominator = [1]
    for root in roots:
        denominator = polymul(denominator, [1, -1/root])
    derivative = [j * denominator[j] for j in range(1, size + 1)]
    factors = [-1/(root * evaluate(derivative, root)) for root in roots]
    return roots, denominator, factors


def decompose(poly, poles):
    roots, denominator, factors = poles
    quotient, remainder = divrem(poly, denominator)
    residues = [factor * evaluate(remainder, root)
                for factor, root in zip(factors, roots)]
    return quotient, residues


def translated_jet(poly, u):
    # Coefficients of 1, h, h^2 in poly(u+h).
    return [sum(comb(n, j) * poly[n] * u**(n - j)
                for n in range(j, len(poly))) for j in range(3)]


def quotient_jet(poly, denominator, u):
    p0, p1, p2 = translated_jet(poly, u)
    d0, d1, d2 = translated_jet(denominator, u)
    z0 = p0 / d0
    z1 = (p1 - d1*z0) / d0
    z2 = (p2 - d1*z1 - d2*z0) / d0
    return z0, u*z1, u*z1 + 2*u*u*z2


def h_parts(v):
    return 1/(v-1), v/(v-1)**2, v*(v+1)/(v-1)**3


def finite_cases(base, size, nodes, mixtures):
    poles = prepare_poles(base, size)
    denominator = poles[1]
    basis = [[1]]
    for k in range(1, size):
        basis.append(polymul(basis[-1], [-base**(-((k+2)//3)), 1]))
    tests = [([0]*k + [1], 0, 0) for k in range(2*size-1)]
    tests += [(polymul(basis[i], basis[j]), max(0, (j+2)//3 - 1), 1)
              for i in range(size) for j in range(i, size)]
    counts = [0, 0, 0]
    for poly, early, kind in tests:
        quotient, residues = decompose(poly, poles)
        jets = [quotient_jet(poly, denominator, base**(-m))
                for m in range(1, nodes+1)]
        for A, B, C in mixtures:
            direct = sum(A*(value-poly[0]) + B*first + C*second
                         for value, first, second in jets)
            # Polynomial quotient plus exact residues, with both finite endpoints shifted.
            split = sum(value * (A+B*j+C*j*j) * (1-base**(-nodes*j))/(base**j-1)
                        for j, value in enumerate(quotient) if j)
            for a, residue in enumerate(residues, 1):
                tail = sum(sum(c*v for c, v in zip((A, B, C), h_parts(base**m)))
                           for m in range(a+1, a+nodes+1))
                split += residue * tail
            assert direct == split
            for value, first, second in jets[:early]:
                assert not any((value, first, second))
                assert A*(value-poly[0])+B*first+C*second == -A*poly[0]
            counts[kind] += 1
            counts[2] += early
    return counts


def base_exponent(k):
    return k*k*(k+1)//2 - sum(comb(max(0, k-2*i)+1, 2) for i in range(k))


def cyclotomic_exponent(k, d):
    if 2*d > k:
        return 5*(k-d)+3
    numerator = 4*k*k + 32*k*d + 17*d*d
    return (numerator + 16*d-1)//(16*d)


def symbolic_entries(size, F, t):
    # Each entry is a linear form in (A,B,C,X) over QQ(t).
    poles = prepare_poles(t, size)
    tail = [F.zero, F.zero, F.zero]
    pole_values = []
    for root in poles[0]:
        tail = [u+v for u, v in zip(tail, h_parts(root))]
        pole_values.append([-v for v in tail] + [F.one])
    entries = []
    for n in range(2*size-1):
        quotient, residues = decompose([0]*n + [1], poles)
        entry = [sum((r*p[j] for r, p in zip(residues, pole_values)), F.zero)
                 for j in range(4)]
        for j, coefficient in enumerate(quotient):
            if j:
                for k in range(3):
                    entry[k] += coefficient * j**k/(t**j-1)
        entries.append(entry)
    return entries


def zero_order(poly):
    return min(m[0] for m in poly.to_dict())


def multiplicity(poly, divisor):
    count = 0
    while True:
        quotient, remainder = poly.div(divisor)
        if remainder:
            return count
        poly = quotient
        count += 1


def symbolic_certificate(size):
    F, t = field("t", QQ)
    T = Symbol("t")
    entries = symbolic_entries(size, F, t)

    determinant = defaultdict(lambda: F.zero)
    for perm in permutations(range(size)):
        inversions = sum(perm[i]>perm[j] for i in range(size)
                         for j in range(i+1, size))
        for selected in product(range(4), repeat=size):
            powers = tuple(selected.count(j) for j in range(4))
            term = prod(entries[i+perm[i]][k] for i, k in enumerate(selected))
            determinant[powers] += (-1)**inversions * term
    determinant = {k: v for k, v in determinant.items() if v}
    e = base_exponent(size)
    phis = [F.from_expr(cyclotomic_poly(d, T)) for d in range(1, size+1)]
    H = prod(phi**cyclotomic_exponent(size, d)
             for d, phi in enumerate(phis, 1))

    order = min(zero_order(v.numer)-zero_order(v.denom) for v in determinant.values())
    losses = [max(0, max(multiplicity(v.denom, phi.numer)-multiplicity(v.numer, phi.numer)
                         for v in determinant.values())) for phi in phis]
    assert order >= e
    assert all(loss <= cyclotomic_exponent(size, d) for d, loss in enumerate(losses, 1))
    at_two = [0]*(size+1)
    for (i, j, k, ell), v in determinant.items():
        normalized = (v*H)/(t**e)
        assert normalized.denom.degree() == 0
        scalar = normalized.denom.to_dict()[(0,)]
        coefficients = {m[0]: c/scalar for m, c in normalized.numer.to_dict().items()}
        assert all(c.denominator==1 for c in coefficients.values())
        at_two[ell] += 2**i*(-3)**j*sum(int(c)*2**n for n, c in coefficients.items())
    if size == 1:
        assert at_two == [-2, 1]
    return dict(K=size, common_order_at_zero=order, guaranteed=e,
                max_poles_by_order=losses, nonzero_monomials=len(determinant)), at_two


def check_five():
    F, t = field("t", QQ)
    poly_ring = t.numer.ring
    phi3 = (t*t+t+1).numer
    size = 5
    _, D, _ = prepare_poles(t, size)
    for i in range(size):
        for j in range(size):
            quotient, _ = divrem([0]*(i+j)+[1], D)
            coefficient = quotient[3] if len(quotient) > 3 else 0
            assert coefficient == (-t**15 if i==j==4 else 0)

    entries = symbolic_entries(size, F, t)
    results = []
    for parameters in [(2,-3,1,5), (3,-2,5,7)]:
        moments = [sum((a*b for a, b in zip(row, parameters)), F.zero) for row in entries]
        common = poly_ring.one
        for value in moments:
            common = common.lcm(value.denom)
        integers = [m.numer*common.exquo(m.denom) for m in moments]
        matrix = [[integers[i+j] for j in range(size)] for i in range(size)]
        numerator = DomainMatrix.from_list(matrix, poly_ring.to_domain()).det()
        delta = F(numerator)/F(common**size)
        pole = multiplicity(delta.denom, phi3)-multiplicity(delta.numer, phi3)
        order = zero_order(delta.numer)-zero_order(delta.denom)
        assert pole <= 4*2+3+1 and order >= base_exponent(size)
        results.append(dict(parameters=parameters, order_at_zero=order, pole_at_Phi3=pole))
    return dict(quotient_rank=1, only_nonzero_entry="-t^15 at (4,4)",
                determinant_specializations=results)


def cube_interval(poly, n=45):
    lo = sum(Q(2, (2**m-1)**3) for m in range(1, n+1))
    hi = lo + Q(16, 7*2**(3*n))
    a = b = 0
    for c in reversed(poly):
        values = [x*y for x in (a, b) for y in (lo, hi)]
        a, b = c+min(values), c+max(values)
    assert a*b > 0
    scale = 10**12
    return [a*scale//1, -((-b)*scale//1)]


def local_costs(max_size):
    tests = 0
    for k in range(1, max_size+1):
        for d in range(1, k+1):
            if 2*d <= k:
                m, extra = divmod(k, d)
                def cost(n):
                    return sum(max(1, n+2-2*j) for j in range(n))
                cost_at_d = extra*cost(m+1)+(d-extra)*cost(m)
            else:
                u = k-d
                cost_at_d = 4*u+3+max(0, u-1)
            assert cost_at_d <= cyclotomic_exponent(k, d)
            tests += 1
    return tests


def degree_samples(sizes):
    maximum = max(sizes)
    phi = list(range(maximum+1))
    for p in range(2, maximum+1):
        if phi[p] == p:
            for i in range(p, maximum+1, p):
                phi[i] -= phi[i]//p
    return {k: sum(cyclotomic_exponent(k,d)*phi[d] for d in range(1,k+1)) for k in sizes}


def main():
    mixtures = [(1,0,0), (0,1,0), (0,0,1), (2,-3,1)]
    official_grid = [(2,1), (2,3), (2,5), (3,5)]
    counts = [finite_cases(Q(b), k, 6, mixtures) for b, k in official_grid]
    combined = list(map(sum, zip(*counts)))
    assert combined == [96,148,40]
    negative = finite_cases(Q(-2), 5, 6, [(3,-1,2)])
    gaussian = finite_cases(QQ_I(1,1), 5, 4, [(QQ_I(1,1),-1,2)])
    symbolic = []
    for size in range(1, 4):
        result, poly = symbolic_certificate(size)
        symbolic.append(result)
    print(json.dumps(dict(finite_counts=combined, extra_negative_counts=negative,
                          extra_gaussian_counts=gaussian, symbolic=symbolic,
                          size_five_quotient_rank_and_specializations=check_five(),
                          cube_P3_at_base_2_polynomial_X=poly,
                          cube_P3_interval_in_units_of_1e_minus12=cube_interval(poly),
                          local_cost_tests_to_size_256=local_costs(256),
                          exact_H_degrees=degree_samples([100,1000,10000])), indent=2))


if __name__ == "__main__":
    main()

