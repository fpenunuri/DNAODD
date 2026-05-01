# DNAODD

## Generalized dual numbers for efficient computation of arbitrary-order directional derivatives in multiple directions

This repository provides a Fortran implementation for the efficient
computation of arbitrary-order directional derivatives in multiple
directions using generalized dual numbers.

------------------------------------------------------------------------

## Overview

This project implements a generalized dual number framework for
computing directional derivatives of arbitrary order along multiple
directions. The approach is designed to be efficient and scalable,
making it suitable for applications in kinematics, optimization, and
automatic differentiation.

The implementation supports the computation of directional derivatives
of arbitrary order along multiple independent directions, based on a
generalized dual number algebra.

------------------------------------------------------------------------

## 📁 Project Structure

-   `src/` --- Fortran source modules
-   `app/` --- Executable programs
-   `fpm.toml` --- Project configuration file

------------------------------------------------------------------------

## 📦 Requirements

-   A Fortran compiler (e.g., `gfortran`, `ifx`)
-   A recent version of `fpm`

To install `fpm`, visit:
👉 https://github.com/fortran-lang/fpm

------------------------------------------------------------------------

## 🛠️ Building with `ifx` (default: `real64`)

From the project root:

``` bash
FPM_FC=ifx fpm build
```

### Quadruple precision

To enable quadruple precision (`real128`), compile with:

``` bash
fpm run --flag "-DUSE_REAL128" <executable_name>
```

For example:

``` bash
fpm run --flag "-DUSE_REAL128" partialD_fun
```

------------------------------------------------------------------------

## ▶️ Running Examples

By default, `MAX_ORDER_DUALZN = 5`.
The example `EA1` requires a minimum order of 7:

``` bash
fpm run --flag "-DMAX_ORDER_DUALZN=7" EA1
```

Other examples can be executed directly, for instance:

``` bash
fpm run RCR_KQs
```

------------------------------------------------------------------------

## Notes

-   The maximum derivative order is controlled at compile time via the
    `MAX_ORDER_DUALZN` flag.
-   Increasing this value allows higher-order derivatives but may
    increase memory usage and computational cost.

------------------------------------------------------------------------

## Relation to the Paper

This repository supports the results presented in the article:

"Efficient Computation of Arbitrary-Order Directional Derivatives in
Multiple Directions via Generalized Dual Numbers"

The code provides implementations and examples used to validate the
proposed methodology.

------------------------------------------------------------------------

## Status

This repository is under active development and is intended to support a
research article currently in preparation/submission.

------------------------------------------------------------------------

## License

(Add a license here, e.g., MIT, BSD, GPL)
