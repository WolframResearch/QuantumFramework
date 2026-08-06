# Quantum potato chips: closed-form answer sheet

*Answer sheet for [quantum-potato-chip-paper-question-chain.md](quantum-potato-chip-paper-question-chain.md): the closed forms of P1, the P2 reconstruction, P3's diagnostic $\phi$, and the chip-preserving noise families of P4. Everything here uses that chain's conventions (the pinned D2 tetrahedron, row coin $p=\tfrac12(1-z/\sqrt{3})$); [quantum-potato-chip-derivation-pra.md](quantum-potato-chip-derivation-pra.md) uses the sign-pattern convention $p=\tfrac12(1+z/\sqrt{3})$, so matching formulas differ by signs.*

## What is the answer in one line?

For a qubit with Bloch vector $\mathbf r=(x,y,z)$, one quantum potato chip is

$$
\boxed{
\mathcal C_y=
\left\{(x,y,z):
y=\frac{xz}{\sqrt3},
\quad
x^2+y^2+z^2\leq1
\right\}.
}
$$

The other two chips are obtained by cycling $x,y,z$:

$$
\boxed{
\mathcal C_x:\ x=\frac{yz}{\sqrt3},
\qquad
\mathcal C_z:\ z=\frac{xy}{\sqrt3}.
}
$$

The name “potato chip” refers to the shape of each curved surface patch inside the Bloch ball.

---

## 1. What is the starting point?

A general qubit state is

$$
\rho=\frac12
\left(
\mathbb I+x\sigma_x+y\sigma_y+z\sigma_z
\right),
$$

where

$$
x^2+y^2+z^2\leq1.
$$

The vector $\mathbf r=(x,y,z)$ is the Bloch vector. Equality gives a pure state; strict inequality gives a mixed state.

---

## 2. Which measurement produces the four probabilities?

Use the standard four-outcome tetrahedral qubit measurement. Its four directions are

$$
\begin{aligned}
\mathbf n_1&=\frac1{\sqrt3}(-1,1,-1),&
\mathbf n_2&=\frac1{\sqrt3}(1,-1,-1),\\
\mathbf n_3&=\frac1{\sqrt3}(-1,-1,1),&
\mathbf n_4&=\frac1{\sqrt3}(1,1,1).
\end{aligned}
$$

The four measurement operators are

$$
Q_i=\frac14\left(\mathbb I+\mathbf n_i\cdot\boldsymbol\sigma\right),
$$

and the four outcome probabilities are

$$
p_i=\operatorname{tr}(\rho Q_i)
=\frac14\left(1+\mathbf r\cdot\mathbf n_i\right).
$$

Therefore,

$$
\begin{aligned}
p_1&=\frac14\left(1+\frac{-x+y-z}{\sqrt3}\right),\\
p_2&=\frac14\left(1+\frac{x-y-z}{\sqrt3}\right),\\
p_3&=\frac14\left(1+\frac{-x-y+z}{\sqrt3}\right),\\
p_4&=\frac14\left(1+\frac{x+y+z}{\sqrt3}\right).
\end{aligned}
$$

These four numbers are nonnegative and sum to one for every physical qubit.

---

## 3. Where do two binary variables enter?

Arrange the four probabilities as a $2\times2$ table:

$$
P=
\begin{pmatrix}
p_1&p_2\\
p_3&p_4
\end{pmatrix}.
$$

This is only a relabeling of the four outcomes. The two row labels form one binary variable, and the two column labels form another.

The two binary variables are independent exactly when

$$
P=
\begin{pmatrix}
p\\
1-p
\end{pmatrix}
\begin{pmatrix}
q&1-q
\end{pmatrix}.
$$

Equivalently,

$$
(p_1,p_2,p_3,p_4)
=
\bigl(
pq,\ p(1-q),\ (1-p)q,\ (1-p)(1-q)
\bigr).
$$

For a $2\times2$ probability table, this happens exactly when its determinant vanishes:

$$
\boxed{p_1p_4-p_2p_3=0.}
$$

---

## 4. How does the determinant give the potato-chip formula?

Insert the four qubit probabilities into the determinant:

$$
\begin{aligned}
p_1p_4-p_2p_3
&=
\frac{\sqrt3\,y-xz}{12}.
\end{aligned}
$$

Independence requires this expression to vanish:

$$
\frac{\sqrt3\,y-xz}{12}=0.
$$

Hence

$$
\boxed{y=\frac{xz}{\sqrt3}.}
$$

This equation gives the full curved surface. Physical qubits must also satisfy the Bloch-ball condition. Substituting $y=xz/\sqrt3$ gives

$$
\boxed{
x^2+z^2+\frac{x^2z^2}{3}\leq1.
}
$$

Thus the reference potato chip is the part of the surface $y=xz/\sqrt3$ that lies inside the Bloch ball.

---

## 5. Why are there three potato chips?

There are three different ways to arrange four probabilities into a $2\times2$ table, apart from swapping rows or columns.

The three determinant conditions are

$$
\begin{aligned}
p_1p_4-p_2p_3&=0,\\
p_1p_3-p_2p_4&=0,\\
p_1p_2-p_3p_4&=0.
\end{aligned}
$$

Substitution gives

$$
\boxed{
y=\frac{xz}{\sqrt3},
\qquad
x=\frac{yz}{\sqrt3},
\qquad
z=\frac{xy}{\sqrt3}.
}
$$

These are the three potato chips associated with the chosen tetrahedral measurement.

---

## 6. What is the closed-form parametrization?

Let the two binary distributions be

$$
(p,1-p),
\qquad
(q,1-q).
$$

Their product gives the four probabilities

$$
\bigl(
pq,\ p(1-q),\ (1-p)q,\ (1-p)(1-q)
\bigr).
$$

Converting these probabilities back to Bloch coordinates gives

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
$$

It is useful to define

$$
a=2p-1,
\qquad
b=2q-1.
$$

Then

$$
\mathbf r=\sqrt3(-b,ab,-a),
$$

and the qubit is physical exactly when

$$
\boxed{
3\left(a^2+b^2+a^2b^2\right)\leq1.
}
$$

This restriction matters: not every pair $(p,q)\in[0,1]^2$ gives a physical qubit.

The boundary of the chip consists of pure states. It is given by

$$
b^2=\frac{1-3a^2}{3(1+a^2)},
\qquad
|a|\leq\frac1{\sqrt3}.
$$

Equivalently,

$$
\boxed{
q_\pm(p)=
\frac12
\left[
1\pm
\sqrt{
\frac{1-3(2p-1)^2}
{3\left(1+(2p-1)^2\right)}
}
\right],
}
$$

with

$$
\boxed{
\frac12\left(1-\frac1{\sqrt3}\right)
\leq p\leq
\frac12\left(1+\frac1{\sqrt3}\right).
}
$$

At $p=q=1/2$, the Bloch vector is zero, so the center of the chip is the maximally mixed state.

---

## 7. Why can two Pauli measurements reconstruct a chip state?

A general qubit needs three Bloch coordinates, so ordinary tomography measures three independent spin directions.

On the reference chip,

$$
y=\frac{xz}{\sqrt3}.
$$

Therefore, measuring only $\sigma_x$ and $\sigma_z$ determines $x$ and $z$, and the chip condition supplies $y$.

The measured probabilities are

$$
\Pr(\sigma_x=\pm1)=\frac{1\pm x}{2},
\qquad
\Pr(\sigma_z=\pm1)=\frac{1\pm z}{2}.
$$

The reconstruction is

$$
\boxed{
y=\frac{xz}{\sqrt3},
\qquad
\rho=
\frac12
\left(
\mathbb I+x\sigma_x+\frac{xz}{\sqrt3}\sigma_y+z\sigma_z
\right).
}
$$

This works only when the state is already known to lie on the chip. Two Pauli measurements do not determine an arbitrary qubit.

The two binary parameters used in the four-outcome table are

$$
p=\frac12\left(1-\frac{z}{\sqrt3}\right),
\qquad
q=\frac12\left(1-\frac{x}{\sqrt3}\right).
$$

---

## 8. How can one test whether a state lies on the chip?

From the four measured probabilities, compute

$$
D=p_1p_4-p_2p_3.
$$

Then

$$
\boxed{
D=0
\quad\Longleftrightarrow\quad
y=\frac{xz}{\sqrt3}.
}
$$

A normalized version is

$$
\boxed{
\phi=
\frac{\sqrt3\,y-xz}
{\sqrt{(3-x^2)(3-z^2)}}.
}
$$

The reference chip is exactly the set of qubit states for which $\phi=0$.

---

## 9. Which common noise processes keep a state on this chip?

Suppose the initial state satisfies $y=xz/\sqrt3$.

Three common noise maps preserve this equation:

$$
\begin{aligned}
\text{bit flip:}\quad
(x,y,z)&\longmapsto(x,\lambda y,\lambda z),\\
\text{phase flip:}\quad
(x,y,z)&\longmapsto(\lambda x,\lambda y,z),\\
\text{dephasing in the }z\text{ basis:}\quad
(x,y,z)&\longmapsto(\mu x,\mu y,z),
\end{aligned}
$$

where $-1\leq\lambda\leq1$ and $0\leq\mu\leq1$.

For example, after a bit flip,

$$
y'=\lambda y
=\lambda\frac{xz}{\sqrt3}
=\frac{x'z'}{\sqrt3}.
$$

Thus the transformed state remains on the same chip. Generic depolarizing noise and amplitude damping do not preserve the full chip.

---

## 10. What can the potato-chip construction be used for?

It provides:

1. a two-parameter family of qubit states with an exact closed-form description;
2. two-setting tomography when chip membership is known in advance;
3. a direct experimental test, $p_1p_4-p_2p_3=0$;
4. a simple model for studying which noise processes preserve a curved family of states;
5. a bridge between a four-outcome quantum measurement and two independent binary probability distributions.

These are structural and experimental simplifications. They do not by themselves imply a computational speedup.

---

## 11. What does the construction not mean?

- The qubit has not been split into two quantum subsystems.
- A chip state is not a bipartite product state.
- The condition $\phi=0$ does not say that all quantum correlations vanish.
- Two Pauli measurements reconstruct only states promised to lie on the chip.
- The construction depends on the chosen tetrahedral measurement and on how its four outcomes are arranged.
- Not every pair of binary probabilities corresponds to a physical qubit.

The precise statement is narrower:

> For one chosen four-outcome qubit measurement, a potato-chip state is a physical qubit whose four outcome probabilities can be arranged as the product of two ordinary binary distributions.
