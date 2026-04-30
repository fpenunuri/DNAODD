!config_mod.F90 
module config_mod
  use iso_fortran_env, only: real64, real128
  implicit none
  private

! --- BEGIN: preprocessor setup ----------------------------------------
#ifndef MAX_ORDER_DUALZN
#define MAX_ORDER_DUALZN 5
#endif

#if (MAX_ORDER_DUALZN < 0)
#  error "MAX_ORDER_DUALZN must be >= 0"
#endif
  
! Precision can be selected by macro: -DUSE_REAL128
#ifdef USE_REAL128
  integer, parameter, public :: prec = real128
#else
  integer, parameter, public :: prec = real64
#endif
! --- END: preprocessor setup ------------------------------------------

  integer, parameter, public :: max_order_dualzn = MAX_ORDER_DUALZN
end module config_mod
