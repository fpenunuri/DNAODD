!Dual numbers to arbitrary order of complex components
!This can be used to compute first, second,... nth order derivatives
!F. Pe~nu~nuri
!UADY, M'erida Yucat'an M'exico
!2025
!
module dualzn_mod
  use iso_fortran_env, only: int32, int64, real32, real64, real128
  use config_mod
  implicit none

  private

  type, public :: dualzn
     complex(prec), dimension(0:max_order_dualzn) :: f
     integer :: ord
  end type dualzn

  public :: initialize_dualzn, f_part, f_set_part
  public :: xto_dzn, xto_complex
  !public :: Dnd  !use this function to "dualize" other functions

  public :: inv, sin, cos, tan, exp, log, sqrt, asin, acos, atan, asinh
  public :: acosh, atanh, sinh, cosh, tanh, absx, atan2, conjg
  public :: matmul, sum, product
  
  !---------------------------------------------------------------------
  ! Operator overloads
  !---------------------------------------------------------------------
  public :: operator(==),  operator(/=)
  public :: operator(+), operator(-), operator(*), operator(/)
  public :: operator(**)

  ! We do not overload intrinsic assignment (=).
  ! Fortran already covers assignment between like types. For cross-type
  ! conversions we prefer an explicit constructor/function.
  !
  ! Use xto_dzn(X, n) to convert a scalar (or array) X into a dualzn of
  ! the desired order n, and then assign the result:
  !
  !   d = xto_dzn(x, n)
  !
  ! Rationale:
  ! - There is no "clean" way to pass the order n through "=".
  ! - It can introduce ambiguities with intrinsic assignment.
  ! - An explicit conversion keeps the API predictable and readable.
  !
  ! Note:
  ! - If desired, the below assign_DCRI_dzn subroutine (or a variation
  !   of it) can be used to overload the "=" operator.
  !
  !public :: assignment(=)   !intentionally not provided
  !---------------------------------------------------------------------

  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dzn
  end interface operator (==)

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dzn
  end interface operator (/=)

  !addition
  interface operator (+)
     module procedure masd      ! unary
     module procedure sumadX    ! dualzn + class
     module procedure sumaXd    ! class + dualzn      
  end interface operator (+)

  !subtraction
  interface operator (-)
     module procedure menosd    ! unary
     module procedure restadX   ! dualzn - class
     module procedure restaXd   ! class - dualzn
  end interface operator (-)

  !multiplication
  interface operator(*)
     module procedure timesdX   ! dualzn * class
     module procedure timesXd   ! class * dualzn     
  end interface operator(*)

  !division
  interface operator(/)
     module procedure divdX     ! dualzn / class
     module procedure divXd     ! class / dualzn
  end interface operator(/)

  !power
  interface operator(**)
     module procedure powerdX   ! dualzn ** class
     module procedure powerXd   ! class ** dualzn
  end interface operator(**)

  !matrix multiplication
  interface matmul
     module procedure MtimesdX
     module procedure MtimesC128d
     module procedure MtimesC64d
     module procedure MtimesC32d
     module procedure MtimesR128d
     module procedure MtimesR64d
     module procedure MtimesR32d
     module procedure MtimesI64d
     module procedure MtimesI32d
  end interface matmul

  !some matrix and vector operations
  interface sum
     module procedure sumR2dzn   !sum(dualzn_rank2,dir)
     module procedure sumR20dzn  !sum(dualzn_rank2)
     module procedure sumR1dzn   !sum(dualzn_rank1)
  end interface sum

  !product for dual 
  interface product
     module procedure prodR2dzn   !product(dual_rank2,dir)
     module procedure prodR20dzn  !product(dual_rank2)
     module procedure prodR1dzn   !product(dual_rank1)
  end interface product

  !overloaded functions
  interface sin
     module procedure sind_
  end interface sin

  interface cos
     module procedure cosd_
  end interface cos

  interface tan
     module procedure tand_
  end interface tan

  interface exp
     module procedure expd
  end interface exp

  interface log
     module procedure logd
  end interface log

  interface sqrt
     module procedure sqrtd
  end interface sqrt

  interface asin
     module procedure asind_
  end interface asin

  interface acos
     module procedure acosd_
  end interface acos

  interface atan
     module procedure atand_
     module procedure atan2d_
  end interface atan

  interface asinh
     module procedure asinhd
  end interface asinh

  interface acosh
     module procedure acoshd
  end interface acosh

  interface atanh
     module procedure atanhd
  end interface atanh

  interface sinh
     module procedure sinhd
  end interface sinh

  interface cosh
     module procedure coshd
  end interface cosh

  interface tanh
     module procedure tanhd
  end interface tanh

  interface atan2
     module procedure atan2d_
  end interface atan2

  interface conjg
     module procedure conjg_dzn
  end interface conjg
  !=====================================================================

  !interface to define a function of complex variable (and an integer)
  !returning a dual variable
  abstract interface
     pure function funzdual(z_val,nrd) result(f_result)
       use config_mod
       import :: dualzn
       complex(prec), intent(in) :: z_val
       integer, intent(in) :: nrd
       type(dualzn) :: f_result
     end function funzdual
  end interface

  ! ====== Functions ======
contains
  elemental subroutine initialize_dualzn(zdn, n)
    type(dualzn), intent(out) :: zdn
    integer, intent(in) :: n

    if (n < 0 .or. n > max_order_dualzn) then
       error stop &
            "initialize_dualzn:  0<=order(ord)<=max_order_dualzn"
    end if

    zdn%ord = n
    zdn%f(0:n) = (0.0_prec, 0.0_prec)
  end subroutine initialize_dualzn

  elemental function f_part(x,k) result(fr)
    type(dualzn), intent(in) :: x
    integer, intent(in) :: k
    complex(prec) :: fr

    if (k < 0 .or. k > x%ord .or. k > max_order_dualzn) then
       error stop &
            "f_part: 0 <= part-k <= 'x%ord' <= max_order_dualzn"
    end if

    fr = x%f(k)
  end function f_part

  elemental subroutine f_set_part(x,y,k)
    type(dualzn), intent(inout) :: x
    complex(prec), intent(in) :: y
    integer, intent(in) :: k

    if (k < 0 .or. k > x%ord .or. k > max_order_dualzn) then
       error stop &
            "f_set_part: 0 <= part-k <= 'x%ord' <= max_order_dualzn"
    end if

    x%f(k) = y
  end subroutine f_set_part

  elemental function xto_dzn(X, n) result(fr)
    use, intrinsic :: ieee_arithmetic, only: &
         ieee_value, ieee_quiet_nan
    
    class(*), intent(in) :: X
    integer,  intent(in) :: n
    type(dualzn) :: fr
    integer :: m

    call initialize_dualzn(fr, n)

    select type (X)
    type is (dualzn)
       m = min(n, X%ord)
       fr%f(0:m) = X%f(0:m)       
    type is (complex(kind=real128))
       fr%f(0) = cmplx(X, kind=prec)
    type is (complex(kind=real64))
       fr%f(0) = cmplx(X, kind=prec)
    type is (complex(kind=real32))
       fr%f(0) = cmplx(X, kind=prec)
    type is (real(kind=real128))
       fr%f(0) = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (real(kind=real64))
       fr%f(0) = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (real(kind=real32))
       fr%f(0) = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (integer(kind=int64))
       fr%f(0) = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (integer(kind=int32))
       fr%f(0) = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    class default
       fr%f(0:n) = cmplx( ieee_value(0.0_prec, ieee_quiet_nan), &
            ieee_value(0.0_prec, ieee_quiet_nan), kind=prec )
    end select
  end function xto_dzn

  elemental function xto_complex(X) result(fr)
    use, intrinsic :: ieee_arithmetic, only: &
         ieee_value, ieee_quiet_nan

    class(*), intent(in) :: X
    complex(prec) :: fr

    select type (X)
    type is (complex(kind=real128))
       fr = cmplx(X, kind=prec)
    type is (complex(kind=real64))
       fr = cmplx(X, kind=prec)
    type is (complex(kind=real32))
       fr = cmplx(X, kind=prec)
    type is (real(kind=real128))
       fr = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (real(kind=real64))
       fr = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (real(kind=real32))
       fr = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (integer(kind=int64))
       fr = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    type is (integer(kind=int32))
       fr = cmplx(real(X, kind=prec), 0.0_prec, kind=prec)
    class default
       fr = cmplx( ieee_value(0.0_prec, ieee_quiet_nan), &
            ieee_value(0.0_prec, ieee_quiet_nan), kind=prec )
    end select
  end function xto_complex
  
  !> Tests whether the numerical value of X is exactly equal to an
  !! integer. The check works for real, complex, and dual numbers.
  !!
  !! Examples:
  !!   x = 1.0                      --> numIntQ(x) = .TRUE.
  !!   x = (1.0, 0.0)               --> numIntQ(x) = .TRUE.
  !!   x = (1.0, 0.0)e0
  !!       + (0.0, 0.0)e1
  !!       + (0.0, 0.0)e2
  !!       + ... + (0.0, 0.0)en     --> numIntQ(x) = .TRUE.
  elemental function numIntQ(X) result (fr)
    class(*), intent(in) :: X
    logical :: fr

    select type(X)
    type is (dualzn)
       fr = xto_dzn(nint(real(x%f(0))),X%ord) == X
    type is (complex(kind=real128))
       fr = nint(real(X,KIND=real128)) == X
    type is (complex(kind=real64))
       fr = nint(real(X,KIND=real64)) == X
    type is (complex(kind=real32))
       fr = nint(real(X,KIND=real32)) == X
    type is (real(kind=real128))
       fr = nint(X) == X
    type is (real(kind=real64))
       fr = nint(X) == X
    type is (real(kind=real32))
       fr = nint(X) == X
    type is (integer(kind=int64))
       fr = .TRUE.
    type is (integer(kind=int32))
       fr = .TRUE.
    class default
       fr = .FALSE.
    end select
  end function numIntQ
  
  !dualzn - class  (B - X)
  elemental function restadX(B, X) result(fr)
    class(*),   intent(in) :: X  
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr

    select type (X)
    type is (dualzn)
       fr = restad(B, X)
    class default
       fr = restad(B, xto_dzn(X, B%ord))
    end select
  end function restadX
  
  !class - dual (X - B)
  elemental function restaXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr

    select type (XX)
    type is (dualzn)
       fr = restad(XX, BB)
    class default
       fr = restad(xto_dzn(XX, BB%ord),BB)
    end select
  end function restaXd
  
  !dualzn + class  (B + X)
  elemental function sumadX(B, X) result(fr)
    class(*),   intent(in) :: X  
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr

    select type (X)
    type is (dualzn)
       fr = sumad(B, X)
    class default
       fr = sumad(B, xto_dzn(X, B%ord))
    end select
  end function sumadX
  
  !class + dual (X + B)
  elemental function sumaXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr

    select type (XX)
    type is (dualzn)
       fr = sumad(XX, BB)
    class default
       fr = sumad(xto_dzn(XX, BB%ord),BB)
    end select
  end function sumaXd
  
  !dualzn*class  (B*X)
  elemental function timesdX(B, X) result(fr)
    class(*),   intent(in) :: X  
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr

    select type (X)
    type is (dualzn)
       fr = timesd(B, X)
    class default
       fr = timesd(B, xto_dzn(X, B%ord))
    end select
  end function timesdX
  
  !class*dual (X*B)
  elemental function timesXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr

    select type (XX)
    type is (dualzn)
       fr = timesd(XX, BB)
    class default
       fr = timesd(xto_dzn(XX, BB%ord),BB)
    end select
  end function timesXd
  
  !dualzn/class  (B/X)
  elemental function divdX(B, X) result(fr)
    class(*),   intent(in) :: X  
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr

    select type (X)
    type is (dualzn)
       fr = divd(B, X)
    class default
       fr = divd(B, xto_dzn(X, B%ord))
    end select
  end function divdX
  
  !class/dual (X/B)
  elemental function divXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr

    select type (XX)
    type is (dualzn)
       fr = divd(XX, BB)
    class default
       fr = divd(xto_dzn(XX, BB%ord),BB)
    end select
  end function divXd
  
  !dual**class  (B**X)
  elemental function powerdX(B, X) result(fr)
    class(*),   intent(in) :: X  
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr

    select type (X)
    type is (dualzn)
       fr = powerd(B, X)
    type is (integer(kind=int64))
       fr = power_dint64(B, X)
    type is (integer(kind=int32))
       fr = power_dint32(B, X)
    class default
       fr = powerd(B, xto_dzn(X, B%ord))
    end select
  end function powerdX
  
  !class**dual (X**B)
  elemental function powerXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr

    select type (XX)
    type is (dualzn)
       fr = powerd(XX, BB)
    class default
       fr = powerd(xto_dzn(XX, BB%ord),BB)
    end select
  end function powerXd
  
  !! A**B
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function powerd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr
    integer(real64) :: iaux    

    if(A%ord /= B%ord) error stop "(**): different orders for operands"
    
    if(numIntQ(B)) then
       iaux = nint(real(B%f(0),kind=prec))
       fr = power_dint64(A,iaux)
    else
       fr = exp(B*log(A))
    end if
  end function powerd
  
  !A**n (n integer)
  elemental function power_dint32(A,n) result(fr)
    type(dualzn), intent(in) :: A
    integer, intent(in) :: n
    type(dualzn) :: fr
    integer :: k, oA

    oA = A%ord
    call initialize_dualzn(fr,oA)

    if(A==xto_dzn(1,oA)) then
       fr=A
       return
    elseif(A==xto_dzn(0,oA) .and. n>0)then
       fr=A
       return
    end if

    !0**0 ---> 1
    if(n==0) then
       fr%f(0) = (1.0_prec,0.0_prec)
    elseif(n>=1) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,n
          fr = fr*A
       end do
    elseif(n<0) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,-n
          fr = fr*A
       end do
       fr = inv(fr)
    end if
  end function power_dint32
  
  !A**n (n integer)
  elemental function power_dint64(A,n) result(fr)
    type(dualzn), intent(in) :: A
    integer(int64), intent(in) :: n
    type(dualzn) :: fr
    integer(int64) :: k
    integer :: oA

    oA = A%ord
    call initialize_dualzn(fr,oA)

    if(A==xto_dzn(1,oA)) then
       fr=A
       return
    elseif(A==xto_dzn(0,oA) .and. n>0)then
       fr=A
       return
    end if

    !0**0 ---> 1
    if(n==0) then
       fr%f(0) = (1.0_prec,0.0_prec)
    elseif(n>=1) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,n
          fr = fr*A
       end do
    elseif(n<0) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,-n
          fr = fr*A
       end do
       fr = inv(fr)
    end if
  end function power_dint64
  
  !> Overloaded assignment (=) for type(dualzn).
  !>
  !> Supported assignments:
  !>   - From another dualzn:
  !>      Behaves like a copy constructor. Copies A%ord and
  !>      A%f(0:A%ord), and zero-fills the remaining elements.
  !>      The semantic order of A is X%ord.
  !>
  !>   - From scalar values (complex, real, or integer; any supported 
  !>      kind): Produces a dualzn of order 0. The scalar is stored 
  !>      in A%f(0). Higher coefficients are zeroed but are semantically
  !>      unused when A%ord == 0 and must not be relied upon.
  !>
  !> Notes:
  !>   - This assignment always fully initializes A%ord and A%f.
  !>   - To create a dual number with a specific order from a
  !>     scalar, use the dedicated constructor, e.g.:
  !>         A = xto_dzn(5, n)
  !>     which yields a dualzn of order n with A%f(0) = 5 and
  !>     A%f(1:n) = 0 (semantically valid up to index n).
  !>
  !>   - Do not read or depend on coefficients beyond A%ord; they are
  !>     considered semantically undefined even if zero-filled here.
  !---------------------------------------------------------------------
  ! elemental subroutine assign_DCRI_dzn(A, X)
  !   class(*), intent(in) :: X
  !   type(dualzn), intent(inout) :: A
    
  !   select type (X)
  !   type is (dualzn)
  !      A%ord = X%ord
  !      A%f(0:X%ord) = X%f(0:X%ord)
  !      A%f(X%ord+1:) = (0.0_prec,0.0_prec)
  !   type is (complex(kind=real128))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, kind=prec)
  !   type is (complex(kind=real64))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, kind=prec)
  !   type is (complex(kind=real32))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, kind=prec)
  !   type is (real(kind=real128))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, 0.0_prec, kind=prec)
  !   type is (real(kind=real64))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, 0.0_prec, kind=prec)
  !   type is (real(kind=real32))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, 0.0_prec, kind=prec)
  !   type is (integer(kind=int64))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, 0.0_prec, kind=prec)
  !   type is (integer(kind=int32))
  !      A%ord = 0
  !      A%f = (0.0_prec,0.0_prec)
  !      A%f(0) = cmplx(X, 0.0_prec, kind=prec)
  !   class default
  !      A%ord = 0
  !      A%f = (0.0_prec, 0.0_prec)
  !   end select
  ! end subroutine assign_DCRI_dzn
  !---------------------------------------------------------------------
  
  !> Logical equality operator.
  !! In this definition, two `dualzn` numbers are considered equal only
  !! if they have the same order and all corresponding components are
  !! equal. 
  elemental function eq_dzn(lhs, rhs) result(fr)
    type (dualzn), intent(in) :: lhs, rhs
    logical :: fr
    logical :: eqfk
    integer :: k

    if (lhs%ord /= rhs%ord) then
       fr = .false.
       return
    end if

    fr = .true.
    do k=0,lhs%ord
       eqfk = lhs%f(k) == rhs%f(k)
       if(.not. eqfk) then
          fr = .false.
          exit
       end if
    end do
  end function eq_dzn
  
  !Logical not equal operator
  elemental function noteq_dzn(lhs, rhs) result(f_res)
    type (dualzn), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dzn
  
  ! +dualzn (unary)
  elemental function masd(A) result(fr)
    type(dualzn), intent(in) :: A
    type(dualzn) :: fr

    fr = A
  end function masd

  !! A+B
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function sumad(A,B) result(fr)
    type(dualzn), intent(in) :: A,B
    type(dualzn) :: fr
    integer :: k

    if(A%ord /= B%ord) error stop "(+): different orders for operands"

    call initialize_dualzn(fr,A%ord)
    do k=0,A%ord
       fr%f(k) = A%f(k) + B%f(k)
    end do
  end function sumad

  ! -dualzn (unary)
  elemental function menosd(A) result(fr)
    type(dualzn), intent(in) :: A
    type(dualzn) :: fr
    integer :: k, orderA

    orderA = A%ord
    call initialize_dualzn(fr,orderA)
    do k=0,orderA
       fr%f(k) = -A%f(k)
    end do
  end function menosd

  !! A-B
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function restad(A,B) result(fr)
    type(dualzn), intent(in) :: A,B
    type(dualzn) :: fr

    if(A%ord /= B%ord) error stop "(-): different orders for operands"
    fr = -B+A
  end function restad
  
  !! A*B
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function timesd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr
    integer :: k

    if(A%ord /= B%ord) error stop "(*): different orders for operands"

    call initialize_dualzn(fr,A%ord)
    do k=0,A%ord
       fr%f(k)=timesdzn(A,B,k)
    end do
  end function timesd

  !! It is assumed A nd B of the same order. Avoid mixing orders.
  pure function timesdzn(A,B,k) result(fr)
    type(dualzn), intent(in) :: A, B
    integer, intent(in) :: k     
    complex(prec) :: fr
    integer :: i

    fr=0.0_prec
    do i=0,k
       fr=fr+binomial(k,i)*A%f(i)*B%f(k-i)
    end do
  end function timesdzn
  !---------------------------------------------------------------------

  !! A/B
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function divd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr

    if(A%ord /= B%ord) error stop "(/): different orders for operands"
    fr = A*inv(B)
  end function divd

  !! inverse multiplicative function
  elemental function inv(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(invzdn,g,k)
    end do
  end function inv

  pure function invzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, signo

    call initialize_dualzn(fr,nrd)

    signo=1
    do k=0,nrd
       fr%f(k)=signo*gamma(real(k+1,kind=prec))/(z**(k+1))
       signo=-signo
    end do
  end function invzdn
  !---------------------------------------------------------------------

  !sin function
  elemental function sind_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(sinzdn,g,k)
    end do
  end function sind_

  pure function sinzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    real(prec), parameter :: Pi = 4.0_prec*atan(1.0_prec)
    integer :: k

    call initialize_dualzn(fr,nrd)

    do k=0,nrd
       fr%f(k) = sin(z + k*Pi/2)
    end do
  end function sinzdn
  !---------------------------------------------------------------------

  ! cos function
  elemental function cosd_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    real(prec), parameter :: Pi = 4.0_prec*atan(1.0_prec)
    type(dualzn) :: gaux
    
    gaux = g
    gaux%f(0) = g%f(0) + Pi/2.0_prec
    fr = sin(gaux)
  end function cosd_
  !---------------------------------------------------------------------

  ! tan function
  elemental function tand_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = sin(g)/cos(g)
  end function tand_
  !---------------------------------------------------------------------

  ! exp function
  elemental function expd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(expzdn,g,k)
    end do
  end function expd

  pure function expzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,nrd)

    do k=0,nrd
       fr%f(k) = exp(z)
    end do
  end function expzdn
  !---------------------------------------------------------------------

  ! log function
  elemental function logd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(logzdn,g,k)
    end do
  end function logd

  pure function logzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, signo

    call initialize_dualzn(fr,nrd)
    fr%f(0) = log(z)

    signo = 1
    do k=1,nrd
       fr%f(k) = signo*gamma(real(k,kind=prec))/(z**k)
       signo = -signo
    end do
  end function logzdn
  !---------------------------------------------------------------------

  ! sqrt function
  elemental function sqrtd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = g**(0.5_prec)
  end function sqrtd
  !---------------------------------------------------------------------

  ! tanh function
  elemental function tanhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = (exp(g) - exp(-g))/(exp(g) + exp(-g))
  end function tanhd
  !---------------------------------------------------------------------
  
  ! cosh function
  elemental function coshd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = 0.5_prec*(exp(g)+exp(-g))
  end function coshd
  !---------------------------------------------------------------------

  ! sinh function
  elemental function sinhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = 0.5_prec*(exp(g)-exp(-g))
  end function sinhd
  !---------------------------------------------------------------------

  ! atanh function
  elemental function atanhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(atanhzdn,g,k)
    end do
  end function atanhd

  pure function atanhzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn

    call initialize_dualzn(fr,nrd)
    fr%f(0) = atanh(z)

    if (nrd == 0) return

    !do not 'simplify' they are complex
    den = 1.0_prec - z*z
    fr%f(1) = 1.0_prec/den
    if (nrd == 1) return

    call initialize_dualzn(auxdn,nrd)
    auxdn%f(0) = den
    auxdn%f(1) = -2.0_prec * z
    auxdn%f(2) = -2.0_prec

    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function atanhzdn
  !---------------------------------------------------------------------

  ! acosh function
  elemental function acoshd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(acoshzdn,g,k)
    end do
  end function acoshd

  pure function acoshzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn, zdn

    call initialize_dualzn(fr, nrd)
    fr%f(0) = acosh(z)
    if (nrd == 0) return

    !do not 'simplify' they are complex
    den = sqrt(z - 1.0_prec) * sqrt(1.0_prec + z)
    fr%f(1) = 1.0_prec/den
    if (nrd == 1) return

    call initialize_dualzn(zdn,nrd)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    auxdn = sqrt(zdn - 1.0_prec) * sqrt(zdn + 1.0_prec)
    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function acoshzdn
  !---------------------------------------------------------------------

  ! asinh function
  elemental function asinhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(asinhzdn,g,k)
    end do
  end function asinhd

  pure function asinhzdn(z, nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn, zdn

    call initialize_dualzn(fr,nrd)
    fr%f(0) = asinh(z)
    if (nrd == 0) return

    fr%f(1) = 1.0_prec/sqrt(1.0_prec + z**2)
    if (nrd == 1) return

    den = sqrt(1.0_prec + z*z)

    call initialize_dualzn(zdn,nrd)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    auxdn = sqrt(1.0_prec + zdn*zdn)

    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function asinhzdn
  !---------------------------------------------------------------------

  ! atan function
  elemental function atand_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(atanzdn,g,k)
    end do
  end function atand_

  pure function atanzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn

    call initialize_dualzn(fr,nrd)
    fr%f(0) = atan(z)
    if (nrd == 0) return

    fr%f(1) = 1.0_prec/(1.0_prec + z**2)
    if (nrd == 1) return

    den = 1.0_prec + z*z

    call initialize_dualzn(auxdn,nrd)
    auxdn%f(0) = den
    auxdn%f(1) = 2.0_prec * z
    auxdn%f(2) = 2.0_prec

    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function atanzdn
  !---------------------------------------------------------------------

  ! acos function
  elemental function acosd_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)   

    do k=0,g%ord
       fr%f(k) = Dnd(acoszdn,g,k)
    end do
  end function acosd_

  pure function acoszdn(z, nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, sqz
    type(dualzn) :: sqzdn, zdn

    call initialize_dualzn(fr,nrd)    
    fr%f(0) = acos(z)
    if (nrd == 0) return

    fr%f(1) = -1.0_prec/sqrt(1.0_prec - z**2)
    if (nrd == 1) return

    call initialize_dualzn(zdn,nrd)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    sqz = sqrt(1.0_prec-z*z)
    sqzdn = sqrt(1.0_prec - zdn*zdn)
    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*sqzdn%f(i)*fr%f(k-i)/sqz
       end do
       fr%f(k) = -sumterm
    end do
  end function acoszdn
  !---------------------------------------------------------------------

  ! asin function
  elemental function asind_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)

    do k=0,g%ord
       fr%f(k) = Dnd(asinzdn,g,k)
    end do
  end function asind_

  pure function asinzdn(z,nrd) result(fr)
    complex(prec), intent(in) :: z
    integer, intent(in) :: nrd
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, sqz
    type(dualzn) :: sqzdn, zdn

    call initialize_dualzn(fr,nrd)

    fr%f(0) = asin(z)
    if (nrd == 0) return

    fr%f(1) = 1.0_prec/sqrt(1.0_prec - z**2)
    if (nrd == 1) return

    call initialize_dualzn(zdn,nrd)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    sqz = sqrt(1.0_prec-z*z)
    sqzdn = sqrt(1.0_prec - zdn*zdn)

    do k=2,nrd
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*sqzdn%f(i)*fr%f(k-i)/sqz
       end do
       fr%f(k) = -sumterm
    end do
  end function asinzdn
  !---------------------------------------------------------------------

  ! conjg
  !notice tat the conjugation operation is not differentiable. In the
  !below definitions we mean (df)* not d(f*)
  elemental function conjg_dzn(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr,g%ord)
    do k=0, g%ord
       fr%f(k) = conjg(g%f(k))
    end do
  end function conjg_dzn
  !---------------------------------------------------------------------

  ! absx = sqrt(z*z) is not sqrt(z*conjg(z))
  elemental function absx(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    complex(prec) :: g0
    integer :: k

    g0 = g%f(0)

    fr = g
    do k=0,g%ord
       fr%f(k) = g%f(k) * g0/sqrt(g0*g0)
    end do
  end function absx
  !---------------------------------------------------------------------

  !! atan2d function
  !! Both operands must have the same order, otherwise execution stops
  !! with an error.
  elemental function atan2d_(y,x) result(fr)
    type(dualzn), intent(in) :: y, x
    type(dualzn) :: fr
    complex(prec) :: x0, y0

    if(y%ord /= x%ord) error stop "(atan2): different orders for operands"
    
    fr = atan(y/x)
    x0 = x%f(0)
    y0 = y%f(0)
    fr%f(0) = atan2_z(y0,x0)
  end function atan2d_

  ! Atan2 for complex arguments
  elemental function atan2_z(zy,zx) result (f_res)
    complex(prec), intent(in) :: zy, zx
    complex(prec) :: f_res
    complex(prec), parameter :: ii = cmplx(0,1,prec)
    complex(prec) :: num, den, divnd, t1, t2
    real(prec) :: x1, x2, y1, y2
    complex(prec) :: x1c, x2c, y1c, y2c

    x1 = real(zx,kind=prec)
    x2 = aimag(zx)

    y1 = real(zy,kind=prec)
    y2 = aimag(zy)

    x1c = x1
    x2c = x2
    y1c = y1
    y2c = y2

    num = x1c + ii*x2c + ii*(y1c + ii*y2c)
    den = sqrt((x1c + ii*x2c)**2 + (y1c + ii*y2c)**2)
    divnd = num/den;
    t1 = atan2(aimag(divnd),real(divnd,kind=prec))
    t2 = ii*log(sqrt((x2c + y1c)**2 + (x1c - y2c)**2)/((2.0_prec*x1c*&
         x2c +  2.0_prec*y1c*y2c)**2 + (x1c**2 - x2c**2 + y1c**2 -   &
         y2c**2)**2)**0.25_prec)

    f_res = t1 - t2
  end function atan2_z
  !---------------------------------------------------------------------

  ! mat-dualzn*class
  function MtimesdX(A,X) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A
    class(*), intent(in), dimension(:,:) :: X
    type(dualzn), dimension(size(A,1),size(X,2)) :: fr
    integer :: nA, mX, oA, oX, k
    
    nA=size(A,2)
    oA = A(1,1)%ord

    mX=size(X,1)

    if (nA /= mX) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)
    
    select type (X)
    type is (dualzn)
       oX = X(1,1)%ord
       if(oA /= oX) then
          print*,"Error: Dualzn matrices must have the same order"
          return
       end if
       fr = Mtimesd(A,X)   
    class default
        do k=0,oA
           fr%f(k) = matmul(f_part(A,k),xto_complex(X))
        end do
    end select
  end function MtimesdX
  !---------------------------------------------------------------------

  function MtimesC128d(C,A) result(fr)
    complex(real128), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)

    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)
    
    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do
  end function MtimesC128d

  function MtimesC64d(C,A) result(fr)
    complex(real64), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)

    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do
  end function MtimesC64d

  function MtimesC32d(C,A) result(fr)
    complex(real32), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)
    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do
  end function MtimesC32d

  function MtimesR128d(C,A) result(fr)
    real(real128), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)
    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do   
  end function 
  
  function MtimesR64d(C,A) result(fr)
    real(real64), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)
    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do
  end function MtimesR64d

  function MtimesR32d(C,A) result(fr)
    real(real32), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)

    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do    
  end function MtimesR32d

  function MtimesI64d(C,A) result(fr)
    integer(real64), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: mA, nC, oA, k

    nC=size(C,2)
    mA=size(A,1)
    oA = A(1,1)%ord

    if (nC /= mA) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul(xto_complex(C),f_part(A,k))
    end do
  end function MtimesI64d

  function MtimesI32d(C,A) result(fr)
    integer(real32), intent(in), dimension(:,:) :: C
    type(dualzn), intent(in), dimension(:,:) :: A
    type(dualzn), dimension(size(C,1),size(A,2)) :: fr
    integer :: oA, k

    oA = A(1,1)%ord

    if (size(C,2) /= size(A,1)) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return          
    end if

    call initialize_dualzn(fr,oA)

    do k=0,oA
       fr%f(k) = matmul( xto_complex(C), f_part(A,k) ) 
    end do
  end function MtimesI32d
  !---------------------------------------------------------------------
  
  ! direct implementation of dual matrix multiplication
  function Mtimesd(A,B) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A, B
    type(dualzn), dimension(size(A,1),size(B,2)) :: fr
    integer :: i, j, k, m, n, p, q, oA, oB

    m = size(A,1)
    n = size(A,2)
    p = size(B,1)
    q = size(B,2)
    oA = A(1,1)%ord
    oB = B(1,1)%ord
    
    if (n /= p) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return
    else if(oA /= oB) then
       print*,"Error: Dualzn matrices must have the same order"
       return  
    end if

    call initialize_dualzn(fr,oA)

    do i = 1, m
       do j = 1, q
          do k = 1, n
             fr(i,j) = fr(i,j) + A(i,k)*B(k,j)
          end do
       end do
    end do
  end function Mtimesd
  !...................

  ! sum(AR2,dir) for rank2 arrays
  ! the result is given as array of rank 1
  function sumR2dzn(AR2,dir) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    integer, intent(in) :: dir
    type(dualzn), dimension(size(AR2,2/dir)) :: fr
    complex(prec), dimension(size(AR2,2/dir)) :: sAk
    complex(prec), dimension(size(AR2,1),size(AR2,2)) :: Ak
    integer :: k, oAR2 
    
    oAR2 = AR2(1,1)%ord
    call initialize_dualzn(fr, oAR2)
    
    do k=0,oAR2
       Ak = f_part(AR2,k)
       sAk = sum(Ak,dir)
       call f_set_part(fr,sAk,k)
    end do
  end function sumR2dzn
  !---------------------------------------------------------------------

  ! sum for Rank 2 array
  function sumR20dzn(AR2) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    type(dualzn) :: fr
    complex(prec), dimension(size(AR2,1),size(AR2,2)) :: Azaux
    integer :: k, oAR2

    oAR2 = AR2(1,1)%ord
    call initialize_dualzn(fr,oAR2)

    do k=0, oAR2
       Azaux = f_part(AR2,k)
       fr%f(k) = sum(Azaux)
    end do
  end function sumR20dzn
  !---------------------------------------------------------------------

  ! sum for Rank 1 array
  function sumR1dzn(AR1) result(fr)
    type(dualzn), intent(in), dimension(:) :: AR1
    type(dualzn) :: fr
    complex(prec), dimension(size(AR1)) :: Azaux
    integer :: k, oAR1

    oAR1 = AR1(1)%ord
    call initialize_dualzn(fr,oAR1)
    do k=0,oAR1
       Azaux = f_part(AR1,k)
       fr%f(k) = sum(Azaux)
    end do
  end function sumR1dzn
  !---------------------------------------------------------------------

  ! product(AR2,dir) for Rank 2 array
  ! the result is given as array of rank 1
  function prodR2dzn(AR2,dir) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    integer, intent(in) :: dir
    type(dualzn), dimension(size(AR2,2/dir)) :: fr
    type(dualzn), dimension(size(AR2,dir)) :: vkdir
    integer :: k

    call initialize_dualzn(fr,AR2(1,1)%ord)
    
    if(dir==1) then
       do k = 1, size(AR2,2)
          vkdir = AR2(:,k)
          fr(k) = prodR1dzn(vkdir)
       end do
    else if(dir==2) then
       do k = 1, size(AR2,1)
          vkdir = AR2(k,:)
          fr(k) = prodR1dzn(vkdir)
       end do
    else 
       stop 'use 1 (2) to collapse rows (columns) in product function'
    end if
  end function prodR2dzn
  
  ! product for Rank 2 array
  function prodR20dzn(AR2) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    type(dualzn) :: fr
    integer :: k, m, orAR2

    orAR2 = AR2(1,1)%ord  

    m=size(AR2,1)
    fr = xto_dzn(1.0_prec,orAR2)
    do k=1,m
       fr = fr*prodR1dzn(AR2(k,:))
    end do
  end function prodR20dzn
  !---------------------------------------------------------------------
  
  ! product for Rank 1 array
  function  prodR1dzn(x) result(fr)
    type(dualzn), intent(in), dimension(:) :: x
    type(dualzn) :: fr
    integer :: k, orx

    orx = x(1)%ord  
    fr = xto_dzn(1.0_prec,orx)
    do k=1,size(x)
       fr = fr*x(k)
    end do
  end function prodR1dzn
  !---------------------------------------------------------------------
  
  !! chain rule, D^n (f(g))
  !! It is declared `recursive` because, when “dualizing” some
  !! functions, it can be convenient to reuse already dualized routines,
  !! even though this function does not directly call itself. 
  recursive pure function Dnd(fc,gdual,n) result(dnfc)
    procedure(funzdual) :: fc
    type(dualzn), intent(in) :: gdual
    integer, intent(in) :: n
    complex(prec) :: dnfc
    type(dualzn) :: fvd
    complex(prec) :: g0, suma
    complex(prec), allocatable, dimension(:) :: xvg
    integer :: k, j

    g0 = gdual%f(0)
    fvd = fc(g0,gdual%ord) 
    if(n==0) then
       dnfc = fvd%f(0)
    else
       suma = 0
       do k=1,n
          allocate(xvg(1:n-k+1))
          do j=1,n-k+1
             xvg(j) = gdual%f(j)
          end do
          suma = suma + fvd%f(k)*BellY(n,k,xvg)
          deallocate(xvg)
       end do
       dnfc = suma
    end if
  end function Dnd
  !---------------------------------------------------------------------

  ! Combinations
  ! m!/((m-n)! * n!) 
  pure function binomial(m, n) result(binom)
    integer, intent(in) :: m, n
    real(prec) :: binom
    integer :: j

    if (n == 0 .or. n == m) then
       binom = 1.0
    else
       binom = 1.0
       do j = 1, n
          binom = binom*(m-j+1)/j
       end do
    endif
  end function binomial

  ! Partial Bell polynomials
  pure function BellY(n, k, x) result(result_value)
    integer, intent(in) :: n, k
    complex(prec), dimension(:), intent(in) :: x
    complex(prec) :: result_value
    complex(prec), dimension(n+1,k+1) :: dp
    complex(prec) :: sum_val
    complex(prec), dimension(:), allocatable :: newx
    integer :: nn, kk, ii, LX, nx

    dp = 0.0_prec
    dp(1, 1) = 1.0_prec  

    do nn = 1, n
       dp(nn+1, 1) = 0.0_prec
    end do

    do kk = 1, k
       dp(1, kk+1) = 0.0_prec
    end do

    !special cases
    if (n == 0 .and. k == 0) then
       result_value = 1.0_prec
       return
    elseif (size(x) == 0) then
       result_value = 0.0_prec
       return
    end if

    !Main loop to compute BellY[n, k, x]
    LX = size(x)
    do nn = 1, n
       do kk = 1, k
          nx = max(nn - kk + 1, LX)
          if (nx > 0) then
             allocate(newx(nx))
             newx = 0.0_prec
             newx(1:LX) = x
             sum_val = 0.0_prec
             do ii = 0, nn - kk
                sum_val = sum_val + binomial(nn - 1, ii) * &
                     newx(ii + 1)*dp(nn - ii, kk)
             end do
             dp(nn + 1, kk + 1) = sum_val
             deallocate(newx)
          end if
       end do
    end do

    result_value = dp(n + 1, k + 1)
  end function BellY
end module dualzn_mod
