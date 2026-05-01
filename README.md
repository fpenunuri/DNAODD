# DNAODD

## Generalized dual numbers for efficient computation of arbitrary-order
directional derivatives in multiple directions

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

- `src/` --- Fortran source modules  
- `app/` --- Executable programs  
- `fpm.toml` --- Project configuration file  

------------------------------------------------------------------------

## 📦 Requirements

- A Fortran compiler, such as `gfortran` or `ifx`
- A recent version of `fpm` is recommended

To install `fpm`, visit:

👉 https://github.com/fortran-lang/fpm

Although this project is designed to be built conveniently with `fpm`,
the use of `fpm` is not strictly required. If `fpm` is not used, the
source files must be compiled and linked manually in the correct order,
taking into account the module dependencies among the Fortran source
files. In addition, preprocessor flags (e.g., `-cpp`, `-DUSE_REAL128`,
`-DMAX_ORDER_DUALZN=<N>`) may be required depending on the desired
configuration.

------------------------------------------------------------------------

## 🛠️ Building with `fpm` (default: `real64`)

From the project root:

```bash
fpm build
```

For example, using Intel Fortran:

```bash
FPM_FC=ifx fpm build
```

### Quadruple precision

To enable quadruple precision (`real128`), compile with:

```bash
fpm run --flag "-DUSE_REAL128" <executable_name>
```

For example:

```bash
fpm run --flag "-DUSE_REAL128" partialD_fun
```

------------------------------------------------------------------------

## ▶️ Running Examples

By default, `MAX_ORDER_DUALZN = 5`.

The example `EA1` requires a minimum order of 7:

```bash
fpm run --flag "-DMAX_ORDER_DUALZN=7" EA1
```

Other examples can be executed directly, for instance:

```bash
fpm run RCR_KQs
```

------------------------------------------------------------------------

## Notes

- The maximum derivative order is controlled at compile time via the
  `MAX_ORDER_DUALZN` flag.
- Increasing this value allows higher-order derivatives but may
  increase memory usage and computational cost.

------------------------------------------------------------------------

## Status

This repository is under active development and is intended to support a
research article currently in preparation and submission.
