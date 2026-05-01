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

    FPM_FC=ifx fpm build

------------------------------------------------------------------------

## ▶️ Running Examples

By default, `MAX_ORDER_DUALZN = 5`. The example `EA1` requires **at
least 7**:

    fpm run --flag "-DMAX_ORDER_DUALZN=7" EA1

Other examples can be executed directly. For instance:

    fpm run RCR_KQs

------------------------------------------------------------------------

## Notes

-   The maximum derivative order is controlled at compile time via the
    `MAX_ORDER_DUALZN` flag.\
-   Increasing this value allows higher-order derivatives but may
    increase memory usage and computational cost.

------------------------------------------------------------------------

## Status

This repository is under active development and is intended to support a
research article currently in preparation/submission.

------------------------------------------------------------------------

