# Trick
some numerical linear algebra tricks

## Diagonal scaling
for linear system $Ax = b$, $A \in R^{N \times N}$, $x, b \in R^N$.
denoted $\tilde{D}$ as the scaling matrix, where
$$
d_i =
\begin{cases}
\sqrt{a_{ii}}, &a_{ii} > 0 \\
1, &a_{ii} = 0 \\
\sqrt{-a_{ii}}, &a_{ii} < 0
\end{cases}
$$
after diagonal scaling, the linear system has converted to $\tilde{A} \tilde{x} = \tilde{b}$
$$
\tilde{A} = \tilde{D}^{-1} A \tilde{D}^{-1} \\
\tilde{b} = \tilde{D}^{-1} b \\
\tilde{x} = \tilde{D} x
$$
element-wise representation
$$
\tilde{a}_{ij} = \frac{a_{ij}}{d_i \ d_j} \\
\tilde{b}_i = \frac{b_i}{d_i} \\
x_i = \frac{\tilde{x}_i}{d_i}
$$