# The weight-three Lambert sum at the golden ratio is irrational

**twinheterodox · 18 November 2026 · Number theory**

Every partial sum of $Z(\varphi)=\sum_{n\ge1}n^2/(\varphi^n-1)$ belongs to $\mathbb Q(\sqrt5)$, where $\varphi=(1+\sqrt5)/2$. The paper proves that its limit lies outside that field.

More generally, let $b$ be an algebraic integer with $|b|>1$. If every conjugate outside the closed unit disk is $b$ or its complex conjugate, then $Z(b)\notin\mathbb Q(b)$. This includes all Pisot and Salem numbers. At the unit-circle conjugates of a Salem number the defining series diverges. A bound on the entire closed disk lets the proof use the same integer polynomials at those conjugates.

At the Gaussian integer $1+i$, every partial sum has rational real and imaginary parts. The conclusion $Z(1+i)\notin\mathbb Q(i)$ means that at least one of these coordinates is irrational in the limit. Tachiya proved the corresponding weight-one result. The numerator here is $n^2$.

[Read the paper (PDF)](output/pdf/algebraic-base-lambert-values.pdf) · [LaTeX source](paper/algebraic-base-lambert-values.tex)

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

The manuscript and explanatory text are offered under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The software is distributed under [Apache-2.0](LICENSE).
