# The weight-three Lambert sum at the golden ratio is irrational

**twinheterodox · 18 November 2026 · Number theory**

Every partial sum of $Z(\varphi)=\sum_{n\ge1}n^2/(\varphi^n-1)$ belongs to $\mathbb Q(\sqrt5)$, where $\varphi=(1+\sqrt5)/2$. The paper proves that its limit lies outside that field.

More generally, let $b$ be an algebraic integer with $|b|>1$. If every conjugate outside the closed unit disk is $b$ or its complex conjugate, then $Z(b)\notin\mathbb Q(b)$. This includes all Pisot and Salem numbers. At the unit-circle conjugates of a Salem number the defining series diverges. A bound on the entire closed disk lets the proof use the same integer polynomials at those conjugates.

At the Gaussian integer $1+i$, every partial sum has rational real and imaginary parts. The conclusion $Z(1+i)\notin\mathbb Q(i)$ means that at least one of these coordinates is irrational in the limit. Tachiya proved the corresponding weight-one result. The numerator here is $n^2$.

[Read the paper (PDF)](output/pdf/algebraic-base-lambert-values.pdf) · [LaTeX source](paper/algebraic-base-lambert-values.tex)

## Companion paper on the height of the base

*A height criterion for irrational weight-three Lambert sums* proves that $Z(b)$ lies outside $\mathbb Q(b)$ when $\log M(b)<(352\pi^2/2889)\,r\log|b|$. Here $M(b)$ is the Mahler measure and $r$ is one for a real base and two for a nonreal base. The criterion gives $Z(3+\sqrt3)\notin\mathbb Q(\sqrt3)$ even though $3-\sqrt3$ also exceeds one. For a nonintegral example, it applies to the larger root of $2x^2-mx+1$ for every integer $m\ge129$.

[Read the companion (PDF)](output/pdf/lambert-height.pdf) · [LaTeX source](paper/lambert-height.tex)

Run `python3 -B checks/check_height.py` from the repository root to reproduce the two height ratios, the exact comparisons in Corollary 1.2, the pole-moment identity, and the leading coefficients in Lemma 4.1 for sizes 1 through 4. [Expected results](checks/height-results.json) · [Review scope](review/lambert-height-scope.json)

## The sum of the reciprocal Mersenne cubes is irrational

Postelmans and Van Assche proved that $1$, the Erdős–Borwein constant $\sum_{n\ge1}(2^n-1)^{-1}$, and the sum of the squared reciprocals are linearly independent over $\mathbb Q$. Adjoining $\sum_{n\ge1}(2^n-1)^{-3}$ raises the dimension to four. The paper proves the same four-value independence for every integer base $b\ge2$. It also proves independence over $\mathbb Q(b)$ when $b$ is an algebraic integer of modulus greater than one and every conjugate other than $b$ and $\bar b$ lies in the closed unit disk.

[Read the paper (PDF)](output/pdf/reciprocal-mersenne-cubes.pdf) · [LaTeX source](paper/reciprocal-mersenne-cubes.tex) · [Verification programs and results](checks/mersenne/) · [Audit scope](review/reciprocal-mersenne-cubes-scope.json)

From the repository root, install the dependencies listed below and run `python3 -B checks/verify.py mersenne` to reproduce all four programs. To build the PDF, run `python3 -B tools/build_paper.py reciprocal-mersenne-cubes`.

## Reproduce the checks

Run these commands from the repository root in a Python environment.

~~~sh
python3 -m pip install -r checks/requirements.txt
python3 -B checks/verify.py
~~~

The first program reconstructs the cleared determinant over the rational function field for sizes 1 through 4. It also checks all coefficients at six integer bases through size 8 and compares two evaluations of the Lambert moments at expanding and contracting bases. The second program checks all eleven affine moments at size 6 and $t=7i/6$ by exact arithmetic in $\mathbb Q(i)$, including the polynomial quotients in Section 4. The runner compares the output of both programs with the [reference results](checks/reference-results.json).

An independent AI model audited the mathematics and exposition of the initial manuscript and supplied the short closed-disk argument used in Section 4. A collaborating model then reviewed the revised seven-page proof and supplied the check at $7i/6$. [Audit scope](review/audit-scope.json) · [Manuscript review](review/final-manuscript-review.txt)

## Build the paper

From the repository root, with pdfLaTeX, Latin Modern and Poppler's `pdfinfo` installed, run

~~~sh
python3 -B tools/build_paper.py
~~~

The bibliography is in the main TeX file. The script runs pdfLaTeX twice and writes the PDF under `output/pdf/`. [Build record](paper/build-record.json)

To build the companion, run `python3 -B tools/build_paper.py lambert-height`. [Companion build record](paper/lambert-height-build-record.json)

The manuscript and explanatory text are offered under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The software is distributed under [Apache-2.0](LICENSE).
