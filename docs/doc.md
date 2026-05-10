# Code Documentation

This document provides a brief description of the main modules, functions, and derived types implemented in the archived repository.

The implementation is organized into three main modules:

- `dir_der_mod`
- `test_functions_mod`
- `fdb_bell_mod`

---

## Module `dir_der_mod`

This module implements routines for computing higher-order directional derivatives, mixed partial derivatives, and Taylor-series-based integration methods using dual numbers.

---

### `dnfvector`

Computes the `n`-th order vector directional derivative along multiple directions.


#### Syntax

```fortran
vec_dir_der = dnfvector(fvecd, A, q, dimf, nd)
```

#### Arguments

- `fvecd`: vector-valued dual-number function, with interface

```fortran
interface
  function fvecdual(xd) result(frd)
    use dualzn_mod
    type(dualzn), intent(in), dimension(:)  :: xd
    type(dualzn), allocatable, dimension(:) :: frd
  end function fvecdual
end interface
```

- `A`: matrix whose rows are the direction vectors along which the directional derivative is computed. These vectors are not necessarily unit vectors.
- `q`: evaluation point.
- `dimf`: dimension of the vector-valued function `fvecd`.
- `nd`: order of the directional derivative.

#### Result

- `vec_dir_der`: value of the vector directional derivative evaluated at `q`.

#### Note

If, for example, the fifth-order vector directional derivative is to be computed along a single direction, do not use a matrix `A` with five identical rows, since this is inefficient. Instead, use a matrix with a single row or call `dnfvector_1vec`.

---

### `dnfscalar`

Computes the `n`-th order scalar directional derivative along multiple directions.

#### Syntax

```fortran
scalar_dir_der = dnfscalar(fsd, A, q, nd)
```

#### Arguments

- `fsd`: scalar-valued dual-number function, with interface

```fortran
interface
  function fsdual(xd) result(frsd)
    use dualzn_mod
    type(dualzn), intent(in), dimension(:) :: xd
    type(dualzn) :: frsd
  end function fsdual
end interface
```

- `A`: matrix whose rows are the direction vectors along which the directional derivative is computed. These vectors are not necessarily unit vectors.
- `q`: evaluation point.
- `nd`: order of the directional derivative.

#### Result

- `scalar_dir_der`: value of the scalar directional derivative evaluated at `q`.

#### Note

If, for example, the fifth-order directional derivative is to be computed along a single direction, do not use a matrix `A` with five identical rows, since this is inefficient. Instead, use a matrix with a single row or call `dnfscalar_1vec`.

---

### `dnfvector_1vec`

Computes the `n`-th order directional derivative of a vector-valued function along a single direction vector.

#### Syntax

```fortran
fr = dnfvector_1vec(fvecd, v, q, dimf, n)
```

#### Arguments

- `fvecd`: vector-valued dual-number function.
- `v`: direction vector.
- `q`: evaluation point.
- `dimf`: dimension of the vector-valued function `fvecd`.
- `n`: order of the directional derivative.

#### Result

- `fr`: value of the `n`-th order directional derivative evaluated at `q`.

---

### `dnfscalar_1vec`

Computes the `n`-th order directional derivative of a scalar-valued function along a single direction vector.

#### Syntax

```fortran
fr = dnfscalar_1vec(fsd, v, q, n)
```

#### Arguments

- `fsd`: scalar-valued dual-number function.
- `v`: direction vector.
- `q`: evaluation point.
- `n`: order of the directional derivative.

#### Result

- `fr`: value of the `n`-th order directional derivative evaluated at `q`.

---

### `basis_vectors`

Constructs a matrix whose rows are standard basis vectors.

#### Syntax

```fortran
fr = basis_vectors(indx_mat, dimEV)
```

#### Arguments

- `indx_mat`: index matrix specifying which basis vectors are generated and how many times each one is repeated.
  - First column: basis vector index.
  - Second column: repetition count.
- `dimEV`: dimension of the vector space.

#### Result

- `fr`: matrix whose rows are the generated basis vectors.

#### Example

```fortran
fr = basis_vectors([[3, 1], [1, 2]], 3)
```

This generates the matrix

```text
[e_3, e_1, e_1]
```

where `e_1` and `e_3` are standard basis vectors in `R^3`.

---

### `SPD`

Computes mixed partial derivatives.

#### Syntax

```fortran
dnfs = SPD(fsd, Aindx, q)
```

#### Arguments

- `fsd`: scalar-valued function to be differentiated.
- `Aindx`: index matrix specifying the differentiation pattern.
  - First column: variable indices with respect to which differentiation is performed.
  - Second column: corresponding derivative orders.
- `q`: point at which the derivative is evaluated.

#### Result

- `dnfs`: value of the specified mixed partial derivative evaluated at `q`.

#### Example

```fortran
Aindx = [[1, 2], [2, 1], [3, 2]]
```

This corresponds to the differentiation pattern

```text
{{x, 2}, {y, 1}, {z, 2}}
```

that is,

```text
partial^5 f / (partial x^2 partial y partial z^2)
```

---

### `TSMDD`

Implements a real-valued Taylor Series Method for solving systems of ordinary differential equations up to fourth order.

#### Syntax

```fortran
Xsol = TSMDD(RHS_ED, xinit_cond, t_vec)
```

#### Arguments

- `RHS_ED`: right-hand side of the ODE system,

```fortran
RHS_ED = Fvec([x1, x2, ..., xn, t])
```

  The function `Fvec` may explicitly depend on time.

- `xinit_cond`: initial conditions.
- `t_vec`: vector of time points where the solution is required.

#### Result

- `Xsol`: matrix containing the computed solution. It does not include the independent variable `t_vec`.

#### Dimensions

- `size(xinit_cond)`: number of equations, or dependent variables.
- `size(t_vec)`: number of time points.

---

## Module `test_functions_mod`

This module contains benchmark and auxiliary functions used in the examples and validation problems.

The module uses the following modules:

- `config_mod`
- `dualzn_mod`
- `diff_mod`

---

### `fstest(q)`

Example of a scalar dual-number function

```text
f : D^m -> D
```

used in one of the derivative-computation examples.

#### Argument

- `q`: evaluation point.

---

### `fvectest(q)`

Example of a vector-valued dual-number function

```text
f : D^m -> D^n
```

#### Argument

- `q`: evaluation point.

---

### `SinProbF(q)`

Sinusoidal problem function used as a benchmark problem.

#### Argument

- `q`: evaluation point.

---

### `rot_mat(th, ev)`

Constructs a rotation matrix over dual numbers.

#### Arguments

- `th`: rotation angle.
- `ev`: axis of rotation.

---

### `HTM(th, ev, Dr)`

Constructs a homogeneous transformation matrix over dual numbers.

#### Arguments

- `th`: rotation angle.
- `ev`: axis of rotation.
- `Dr`: displacement vector.

---

### `RCR_rD(q)`

Computes the position vector of the end effector of the RCR robot manipulator used in the examples.

#### Argument

- `q`: generalized coordinates.

---

## Module `fdb_bell_mod`

This module implements combinatorial utilities related to the Faa di Bruno formula and partial Bell polynomials.

If Wolfram Mathematica is available, this module may be omitted in favor of symbolic generation methods.

---

### `generate_fdb_bell`

Enumerates all admissible integer tuples `(m1, ..., mn)` that contribute to the `n`-th order Faa di Bruno / partial Bell polynomial expansion and returns a list of term descriptors.

#### Syntax

```fortran
terms = generate_fdb_bell(n)
```

#### Argument

- `n`: order of the expansion or derivative.

#### Result

- `terms`: allocatable array of type `term_t` containing all admissible terms for the given order `n`.

---

### `print_terms`

Prints a list of Faa di Bruno / Bell polynomial term descriptors.

#### Syntax

```fortran
call print_terms(n, terms)
```

#### Output format

```text
{[coeff, k], [[j, m], ...]}
```

where:

- `coeff`: factorial coefficient.
- `k`: total multiplicity, `sum m_j`.
- `(j, m)`: nonzero index-multiplicity pair, with `m_j > 0`.

#### Arguments

- `n`: order associated with the term list. This is used only for labeling the output.
- `terms`: array of type `term_t`, as produced by `generate_fdb_bell`.

---

## Derived Type `pair_t`

Stores a single `(j, m_j)` pair used in the multivariate Faa di Bruno / partial Bell polynomial representation.

### Fields

- `j`: index of the derivative group, corresponding to the subscript `j` in `m_j`.
- `m`: multiplicity associated with `j`. Only pairs with `m > 0` are stored in a term descriptor.

---

## Derived Type `term_t`

Describes one admissible term in the `n`-th order Faa di Bruno combinatorial expansion.

### Stored information

Each term stores:

- the factorial coefficient,

```text
coeff = n! / product_{j=1..n} [m_j! (j!)^{m_j}]
```

- the total number of directional arguments,

```text
k = sum_{j=1..n} m_j
```

- the sparse list of pairs `(j, m_j)` for which `m_j > 0`.

### Fields

- `k`: total multiplicity, or number of directional arguments.
- `coeff_int64`: exact integer coefficient if representable in `int64`. If there is an overflow risk, this value is set to `-1`.
- `log_coeff`: natural logarithm of the coefficient. This value is always valid and is computed using log-gamma functions for numerical stability.
- `pairs`: allocatable array of type `pair_t` containing only the indices with `m_j > 0`.
