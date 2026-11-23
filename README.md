# Sharp thresholds for anonymous recovery in Hamming graphs

Brandon Li

Paper and verification code on random landmark sets in the Hamming graph
$H(n,q)$ that identify an unknown query from the unlabeled histogram of its
distances to the landmarks (multiset resolving sets).

For the binary cube, a random landmark set of expected size $e^{An^{1/3}}$
distinguishes all $3^n$ faces with high probability when $A>C$ and fails
when $A\le C$, where

$$
C=\max_{0<\rho\le1}\left(\frac{3h(\rho)}{\sqrt{\rho/2}}\right)^{2/3}=3.4417465\ldots,
\qquad h(\rho)=-\rho\log\rho-(1-\rho)\log(1-\rho)+\rho\log2.
$$

The maximum is attained at $\rho^*=0.3655\ldots$, so the faces that decide
the threshold fix about 37% of the coordinates.

Under an adversarial error budget of $e^{Bn^{1/3}}$ with $B>0$, recovery
holds with high probability exactly when $A-2B\ge C$. The paper also finds the sharp polynomial factors at the second
scale and the faces that cause the rare failures above the threshold.

For every fixed alphabet size $q$, the paper proves sharp random thresholds
for resolving all vertices, outer vertices and edges, the critical windows
and limiting recovery probabilities for vertices and edges, and the same
results for uniformly chosen landmark sets of a fixed size. The vertex,
outer-vertex and edge multiset dimensions $D_n$ of $H(n,q)$ satisfy

$$
\frac{C_q}{4^{2/3}} \le \liminf_{n\to\infty} \frac{\log D_n}{n^{1/3}}
\le \limsup_{n\to\infty} \frac{\log D_n}{n^{1/3}} \le C_q,
\qquad C_q=\left(\frac{3\log q}{\sqrt{2\sigma^2}}\right)^{2/3},\quad \sigma^2=\frac{q-1}{q^2},
$$

which answers Questions 2 and 4 of Allikvere on the edge multiset dimension
of hypercubes. For binary edges, the same constants were found independently
in [eliottcassidy2000/math](https://github.com/eliottcassidy2000/math)
(October 2026); the paper extends them to every fixed alphabet and all three
query families and makes the random upper bound sharp.

- [Paper (PDF)](paper/main.pdf)
- [LaTeX source](paper/)
- [Verification code](checks/)

## Reproduce

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python checks/verify_constants.py
.venv/bin/python checks/verify_poisson_window.py
```

`verify_constants.py` recomputes the binary-face constants to 60 digits
with the standard library. `verify_poisson_window.py` recomputes the vertex
and edge Poisson-window constants with NumPy and SciPy. Their outputs are in
[`checks/`](checks/). The PDF builds with Tectonic
(`cd paper && tectonic main.tex`) or `pdflatex`.

## License

Code: MIT. Paper: CC BY 4.0.
