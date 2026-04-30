!local module for function used in the example
module loc_fun_mod
  use config_mod
  use dualzn_mod
  use dir_der_mod
  implicit none
  private

  public F1_rhs, r_linspace

  !Interface for a vector dual function f: D^m --> D^n
  interface
     function fvecdual(xd) result(frd)
       use dualzn_mod
       type(dualzn), intent(in), dimension(:)  :: xd
       type(dualzn), allocatable, dimension(:) :: frd
     end function fvecdual
  end interface

contains
  ! an ODE's 
  ! y1' = 2*x*exp(-x**2)*(1-x*cos y2) + sin y2
  ! y2' = -2*x*exp(-x**2)
  ! y1(0) = -2
  ! y2(0) = 2  
  ! F = f_rhs([y1,y2,x])
  !
  ! exact solution
  ! y1x(x) = x*sin(exp(-x**2) + 1) - exp(-x**2) - 1
  ! y2x(x) = exp(-x**2) + 1
  function F1_rhs(qvec) result(frvec)
    type(dualzn), intent(in), dimension(:) :: qvec
    type(dualzn), allocatable, dimension(:) :: frvec
    type(dualzn) :: y1, y2, x

    allocate(frvec(2))

    y1 = qvec(1)
    y2 = qvec(2)
    x  = qvec(3) !independent variable 

    frvec(1) = 2.0_prec*x*exp(-x**2)*(1.0_prec - x*cos(y2)) + sin(y2)
    frvec(2) = -2*x*exp(-x**2)
  end function F1_rhs
  
  ! generates n equaly spaced points in [a,b]
  function r_linspace(a,b,n) result(fr)
    real(prec), intent(in) :: a, b
    integer, intent(in) :: n
    real(prec), dimension(n) :: fr
    real(prec) :: dx
    integer :: i

    if (n < 1) then
       write(*,'(A)') "error (r_linspace) must have n > 1"
       stop 1
    elseif (n==2) then
       fr = [a,b]
       return
    end if

    dx = (b - a) / real(n - 1, prec)
    fr(1) = a
    fr(n) = b
    do i=2, n-1
       fr(i) = a + dx * real(i - 1, prec)
    end do
  end function r_linspace
 end module loc_fun_mod


!main program
program E2A
  use config_mod
  use dualzn_mod
  use dir_der_mod
  use loc_fun_mod
  implicit none

  integer, parameter :: npts = 10
  real(prec), dimension(npts) :: xvec
  real(prec), dimension(2) :: y0
  real(prec), parameter :: a = 0.0_prec, b = 3.0_prec
  real(prec), dimension(npts,size(y0)) :: ysol
  real(prec), dimension(size(y0)) :: ysol_vec
  integer :: i
  
  xvec = r_linspace(a,b,npts)
  y0 = [-2.0_prec, 2.0_prec]
  ysol = TSMDD(F1_rhs, y0, xvec)

  open(unit=15, file="ysolF1.dat", status="unknown")
  do i = 1, npts
     ysol_vec = ysol(i,:)
     write(15,"(*(F8.5,1x))") ysol_vec
  end do
  close(15)

  write(*,"(A)"), 'data saved to "ysolF1.dat"'
end program E2A
