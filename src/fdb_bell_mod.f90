!====================================================================
!  fdb_bell_mod.f90
!====================================================================
!  - Computes Faa di Bruno combinatorial coefficients using log-gamma
!    for numerical stability.
!  - Enumerates all admissible integer tuples (m1,...,mn) satisfying
!    Sum j*m_j = n and builds the corresponding term descriptors.
!  - Output format compatible with the symbolic Mathematica version:
!      {[coeff,k],[[j,m],...]}
!    where:
!      coeff : factorial coefficient n! / (Prod m_j! (j!)^m_j)
!      k     : total number of directional arguments (Sum m_j)
!      (j,m) : index-multiplicity pairs with m > 0
!
!  This module provides a fully numerical alternative to the
!  symbolic construction based on partial Bell polynomials (BellY).
!  It is suitable for low-level implementations of the
!  multivariate Faa di Bruno expansion, e.g., in evaluating
!  higher-order directional derivatives of vector functions.
!====================================================================
module fdb_bell_mod
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private
  public :: term_t, pair_t, generate_fdb_bell, print_terms

  type :: pair_t
     integer :: j = 0
     integer :: m = 0
  end type pair_t

  type :: term_t
     integer        :: k = 0
     integer(int64) :: coeff_int64 = -1_int64   ! exact if fits; otherwise -1
     real(8)        :: log_coeff   = 0.0d0      ! always valid
     type(pair_t), allocatable :: pairs(:)
  end type term_t

contains
  function log_fact(n) result(lf)
    integer, intent(in) :: n
    real(8) :: lf

    lf = log_gamma(real(n+1,8))   ! log(n!) = log_gamma(n+1)
  end function log_fact

  function log_bell_coeff(n, mvec) result(lc)
    integer, intent(in) :: n
    integer, intent(in) :: mvec(:)   ! size >= n
    real(8) :: lc
    integer :: j

    lc = log_fact(n)
    do j = 1, n
       if (mvec(j) > 0) then
          lc = lc - log_fact(mvec(j)) - real(mvec(j),8)*log_fact(j)
       end if
    end do
  end function log_bell_coeff

  function fact64(n) result(f)    
    integer, intent(in) :: n
    integer(int64) :: f
    integer :: i

    f = 1_int64
    if (n <= 1) return
    do i = 2, n
       f = f * int(i, int64)
    end do
  end function fact64

  function powi64(a, e) result(p)
    integer(int64), intent(in) :: a
    integer, intent(in) :: e
    integer(int64) :: p
    integer :: i

    p = 1_int64
    if (e <= 0) return
    do i = 1, e
       p = p * a
    end do
  end function powi64

  subroutine bell_coeff_int64(n, mvec, c_out)
    integer, intent(in) :: n
    integer, intent(in) :: mvec(:)
    integer(int64), intent(out) :: c_out
    integer :: j
    integer(int64) :: denom

    c_out = fact64(n)
    denom = 1_int64
    do j = 1, n
       if (mvec(j) > 0) then
          denom = denom * fact64(mvec(j)) * powi64(fact64(j), mvec(j))
       end if
    end do
    c_out = c_out / denom
  end subroutine bell_coeff_int64

  subroutine push_term(buf, used, t)
    type(term_t), allocatable, intent(inout) :: buf(:)
    integer, intent(inout) :: used
    type(term_t), intent(in) :: t
    type(term_t), allocatable :: tmp(:)
    integer :: oldcap, newcap

    if (.not. allocated(buf)) then
       allocate(buf(8))
       used = 0
    end if

    if (used >= size(buf)) then
       oldcap = size(buf)
       newcap = max(8, 2*oldcap)
       allocate(tmp(newcap))
       if (used > 0) tmp(1:used) = buf(1:used)
       call move_alloc(tmp, buf)
    end if

    used = used + 1
    buf(used) = t
  end subroutine push_term

  function generate_fdb_bell(n) result(terms)
    integer, intent(in) :: n
    type(term_t), allocatable :: terms(:)
    type(term_t), allocatable :: buf(:)
    integer :: used
    integer :: k
    integer, allocatable :: mvec(:)

    used = 0
    allocate(mvec(n))
    mvec = 0

    do k = 1, n
       mvec = 0
       call enumerate_local(j=n, rem_n=n, rem_k=k)
    end do

    if (used == 0) then
       allocate(terms(0))
    else
       allocate(terms(used))
       terms = buf(1:used)
    end if

  contains
    subroutine emit_local(mv, k_local)
      integer, intent(in) :: mv(:)
      integer, intent(in) :: k_local

      type(term_t) :: t
      integer :: cnt, j, p
      real(8) :: lc, log_huge
      integer(int64) :: c64

      ! Build pair list
      cnt = 0
      do j = 1, n
         if (mv(j) > 0) cnt = cnt + 1
      end do

      allocate(t%pairs(cnt))
      p = 0
      do j = 1, n
         if (mv(j) > 0) then
            p = p + 1
            t%pairs(p)%j = j
            t%pairs(p)%m = mv(j)
         end if
      end do

      t%k = k_local
      lc = log_bell_coeff(n, mv)   ! log coefficient (robust)
      t%log_coeff = lc

      ! Default: not representable exactly in int64
      t%coeff_int64 = -1_int64

      ! Safe exact int64 only if n! fits (n <= 20) and the final coefficient fits
      log_huge = log(real(huge(0_int64), 8))
      if (n <= 20 .and. lc <= log_huge) then
         call bell_coeff_int64(n, mv, c64)
         t%coeff_int64 = c64
      end if

      call push_term(buf, used, t)
    end subroutine emit_local

    recursive subroutine enumerate_local(j, rem_n, rem_k)
      integer, intent(in) :: j, rem_n, rem_k
      integer :: up, val

      if (j == 0) then
         if (rem_n == 0 .and. rem_k == 0) call emit_local(mvec, k)
         return
      end if

      up = min(rem_n / j, rem_k)
      do val = up, 0, -1
         mvec(j) = val
         call enumerate_local(j-1, rem_n - j*val, rem_k - val)
      end do
      mvec(j) = 0
    end subroutine enumerate_local
  end function generate_fdb_bell

  subroutine print_terms(n, terms)
    integer, intent(in) :: n
    type(term_t), intent(in) :: terms(:)
    integer :: i, p
    real(8) :: lg10, frac, mant
    integer :: e10

    write(*,'(a,i0,a)') '--- Faa di Bruno Bell Terms for n = ', n, ' ---'
    do i = 1, size(terms)
       write(*,'(a)', advance='no') '{['

       if (terms(i)%coeff_int64 >= 0_int64) then
          write(*,'(i0)', advance='no') terms(i)%coeff_int64
       else
          lg10 = terms(i)%log_coeff / log(10.0d0)
          e10  = int(floor(lg10))
          frac = lg10 - real(e10, 8)
          mant = 10.0d0**frac

          ! Guard against rounding that could produce mant == 10.0
          if (mant >= 10.0d0) then
             mant = mant / 10.0d0
             e10  = e10 + 1
          end if

          write(*,'(f10.6,a,i0)', advance='no') mant, 'E', e10
       end if

       write(*,'(a)', advance='no') ','
       write(*,'(i0)', advance='no') terms(i)%k
       write(*,'(a)', advance='no') '],['

       do p = 1, size(terms(i)%pairs)
          if (p > 1) write(*,'(a)', advance='no') ','
          write(*,'(a)', advance='no') '['
          write(*,'(i0)', advance='no') terms(i)%pairs(p)%j
          write(*,'(a)', advance='no') ','
          write(*,'(i0)', advance='no') terms(i)%pairs(p)%m
          write(*,'(a)', advance='no') ']'
       end do

       write(*,'(a)') ']}'
    end do
  end subroutine print_terms
end module fdb_bell_mod

