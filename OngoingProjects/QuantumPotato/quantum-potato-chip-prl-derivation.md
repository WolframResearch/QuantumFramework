# A Closed-Form Factorization Surface in the Qubit Bloch Ball

## Abstract

We identify the qubit states for which the four probabilities of a tetrahedral symmetric informationally complete measurement can be arranged as the product of two binary distributions. For a fixed arrangement of the four outcomes, classical independence is equivalent to the vanishing of a $2\times2$ determinant. In Bloch coordinates this condition reduces exactly to $y=xz/\sqrt3$. Its intersection with the Bloch ball is a compact curved surface whose boundary consists of pure states. Two other surfaces follow from the two inequivalent rearrangements of the four outcomes. We give the closed-form parametrization and show that, under the promise that the state belongs to one such surface, two Pauli measurement settings determine the full density matrix.

---

## Derivation

**Step 1—Write the qubit in Bloch form.**
A qubit state is

$$
\rho(\mathbf r)
=
\frac12
\left(
\mathbb I+\mathbf r\cdot\boldsymbol\sigma
\right),
\qquad
\mathbf r=(x,y,z),
\qquad
|\mathbf r|^2\leq1.
\tag{1}
$$

Here $\boldsymbol\sigma=(\sigma_x,\sigma_y,\sigma_z)$. The inequality in Eq. (1) is equivalent to $\rho\geq0$; equality holds for pure states.

**Step 2—Choose a four-outcome tetrahedral measurement.**
Let

$$
\begin{aligned}
\mathbf n_1&=\frac1{\sqrt3}(-1,1,-1),&
\mathbf n_2&=\frac1{\sqrt3}(1,-1,-1),\\
\mathbf n_3&=\frac1{\sqrt3}(-1,-1,1),&
\mathbf n_4&=\frac1{\sqrt3}(1,1,1).
\end{aligned}
\tag{2}
$$

These unit vectors obey

$$
\mathbf n_i\cdot\mathbf n_j=-\frac13
\quad(i\neq j),
\qquad
\sum_{i=1}^4\mathbf n_i=0.
\tag{3}
$$

The measurement effects

$$
E_i=\frac14
\left(
\mathbb I+\mathbf n_i\cdot\boldsymbol\sigma
\right)
\tag{4}
$$

are positive and satisfy $\sum_iE_i=\mathbb I$. Their outcome probabilities are

$$
p_i=\operatorname{tr}(\rho E_i)
=
\frac14
\left(
1+\mathbf r\cdot\mathbf n_i
\right).
\tag{5}
$$

Equations (2) and (5) give

$$
\begin{aligned}
p_1&=\frac14\left(1+\frac{-x+y-z}{\sqrt3}\right),&
p_2&=\frac14\left(1+\frac{x-y-z}{\sqrt3}\right),\\
p_3&=\frac14\left(1+\frac{-x-y+z}{\sqrt3}\right),&
p_4&=\frac14\left(1+\frac{x+y+z}{\sqrt3}\right).
\end{aligned}
\tag{6}
$$

**Step 3—Read the four outcomes as two binary labels.**
Arrange the probabilities as

$$
P=
\begin{pmatrix}
p_1&p_2\\
p_3&p_4
\end{pmatrix}.
\tag{7}
$$

The row and column labels define two ordinary binary variables. They are independent if and only if

$$
P=
\begin{pmatrix}
p\\
1-p
\end{pmatrix}
\begin{pmatrix}
q&1-q
\end{pmatrix},
\tag{8}
$$

or, equivalently,

$$
(p_1,p_2,p_3,p_4)
=
\bigl(
pq,\ p(1-q),\ (1-p)q,\ (1-p)(1-q)
\bigr).
\tag{9}
$$

For a normalized nonnegative $2\times2$ table, Eq. (8) is equivalent to

$$
\det P=p_1p_4-p_2p_3=0.
\tag{10}
$$

**Step 4—Evaluate the determinant.**
Define

$$
\begin{aligned}
A&=-x+y-z,&B&=x-y-z,\\
C&=-x-y+z,&D&=x+y+z.
\end{aligned}
\tag{11}
$$

Using Eq. (6),

$$
16\det P
=
\left(1+\frac{A}{\sqrt3}\right)
\left(1+\frac{D}{\sqrt3}\right)
-
\left(1+\frac{B}{\sqrt3}\right)
\left(1+\frac{C}{\sqrt3}\right).
\tag{12}
$$

The linear and quadratic terms reduce separately:

$$
A+D-B-C=4y,
\qquad
AD-BC=-4xz.
\tag{13}
$$

Therefore,

$$
16\det P
=
\frac{4y}{\sqrt3}-\frac{4xz}{3},
\qquad
\det P
=
\frac{\sqrt3\,y-xz}{12}.
\tag{14}
$$

Combining Eqs. (10) and (14) gives the factorization surface

$$
\boxed{
y=\frac{xz}{\sqrt3}.
}
\tag{15}
$$

**Step 5—Impose quantum-state positivity.**
Equation (15) describes an unbounded surface in $\mathbb R^3$. A physical qubit must also satisfy Eq. (1). Substitution of Eq. (15) yields

$$
\boxed{
x^2+z^2+\frac{x^2z^2}{3}\leq1.
}
\tag{16}
$$

The quantum potato chip is therefore

$$
\boxed{
\mathcal C_y
=
\left\{
\rho(x,y,z):
y=\frac{xz}{\sqrt3},
\quad
x^2+y^2+z^2\leq1
\right\}.
}
\tag{17}
$$

Its boundary is obtained by replacing the inequality in Eq. (16) with equality. Hence every boundary point has $|\mathbf r|=1$ and represents a pure state.

**Step 6—Obtain the other two chips.**
The four outcomes admit three inequivalent $2\times2$ arrangements. Their determinants are

$$
\begin{aligned}
p_1p_4-p_2p_3
&=\frac{\sqrt3\,y-xz}{12},\\
p_1p_3-p_2p_4
&=\frac{yz-\sqrt3\,x}{12},\\
p_1p_2-p_3p_4
&=\frac{xy-\sqrt3\,z}{12}.
\end{aligned}
\tag{18}
$$

Thus the three surfaces are

$$
\boxed{
\mathcal C_y:\ y=\frac{xz}{\sqrt3},
\qquad
\mathcal C_x:\ x=\frac{yz}{\sqrt3},
\qquad
\mathcal C_z:\ z=\frac{xy}{\sqrt3},
}
\tag{19}
$$

with the Bloch-ball condition understood in each case.

**Step 7—Find a closed parametrization.**
The inverse of Eq. (6) is

$$
\begin{aligned}
x&=\sqrt3(1-2p_1-2p_3),\\
y&=\sqrt3(1-2p_2-2p_3),\\
z&=\sqrt3(1-2p_1-2p_2).
\end{aligned}
\tag{20}
$$

Substituting the product probabilities in Eq. (9) gives

$$
\boxed{
\mathbf r(p,q)
=
\sqrt3
\left(
1-2q,\,
(2p-1)(2q-1),\,
1-2p
\right).
}
\tag{21}
$$

Let

$$
a=2p-1,
\qquad
b=2q-1.
\tag{22}
$$

Then

$$
\mathbf r=\sqrt3(-b,ab,-a),
\qquad
|\mathbf r|^2
=
3\left(a^2+b^2+a^2b^2\right).
\tag{23}
$$

The exact physical domain is therefore

$$
\boxed{
3\left(a^2+b^2+a^2b^2\right)\leq1.
}
\tag{24}
$$

The pure-state boundary satisfies

$$
\boxed{
b^2=\frac{1-3a^2}{3(1+a^2)},
\qquad
|a|\leq\frac1{\sqrt3}.
}
\tag{25}
$$

In terms of the original parameters,

$$
\boxed{
q_\pm(p)
=
\frac12
\left[
1\pm
\sqrt{
\frac{1-3(2p-1)^2}
{3\left[1+(2p-1)^2\right]}
}
\right],
}
\tag{26}
$$

where

$$
\frac12
\left(
1-\frac1{\sqrt3}
\right)
\leq p\leq
\frac12
\left(
1+\frac1{\sqrt3}
\right).
\tag{27}
$$

At $p=q=1/2$, Eq. (21) gives $\mathbf r=0$, the maximally mixed state.
At $a=\pm1/\sqrt3$, $b=0$, the boundary reaches the two $\sigma_z$ eigenstates; at $a=0$, $b=\pm1/\sqrt3$, it reaches the two $\sigma_x$ eigenstates.

**Step 8—Recover the state from two Pauli settings.**
The row and column sums of Eq. (7) are

$$
p=p_1+p_2
=
\frac12
\left(
1-\frac{z}{\sqrt3}
\right),
\qquad
q=p_1+p_3
=
\frac12
\left(
1-\frac{x}{\sqrt3}
\right).
\tag{28}
$$

Sharp Pauli measurements give

$$
\Pr(\sigma_x=\pm1)=\frac{1\pm x}{2},
\qquad
\Pr(\sigma_z=\pm1)=\frac{1\pm z}{2}.
\tag{29}
$$

For a state promised to lie in $\mathcal C_y$, these two settings determine $x$ and $z$, while Eq. (15) fixes the unmeasured coordinate. The reconstructed density matrix is

$$
\boxed{
\rho_{\mathcal C_y}
=
\frac12
\left(
\mathbb I
+x\sigma_x
+\frac{xz}{\sqrt3}\sigma_y
+z\sigma_z
\right).
}
\tag{30}
$$

This is a constrained tomography result: two Pauli settings determine a chip state, but they do not determine an arbitrary qubit.

---

## Physical interpretation

The construction identifies a two-dimensional family of qubit states whose probabilities under one fixed four-outcome measurement form a rank-one $2\times2$ table. The factorization is classical independence between two labels assigned to the measurement outcomes. It is not Hilbert-space separability: the system is a single qubit and has no bipartite tensor-product structure.

The result is also representation dependent. Rotating the tetrahedral measurement rotates the three chips, and rearranging its outcomes selects which of the three equations in Eq. (19) applies. What is invariant is the construction: a physical qubit state, a chosen four-outcome tetrahedral measurement, and a vanishing determinant.
