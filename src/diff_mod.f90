!some directional derivatives and differential operators
module diff_mod
  use config_mod
  use dualzn_mod
  implicit none

  private
  public :: d1fscalar, d2fscalar, d1fvector
  public :: Gradient, Jacobian, Hessian

  !Interface for a scalar dual function f: D^m --> D (similar to
  !f: R^m --> R)
  abstract interface
     function fsdual(xd) result(frsd)
       use dualzn_mod
       type(dualzn), intent(in), dimension(:) :: xd
       type(dualzn) :: frsd
     end function fsdual
  end interface

  !Interface for a vector dual function f: D^m --> D^n
  interface
     function fvecdual(xd) result(frd)
       use dualzn_mod
       type(dualzn), intent(in), dimension(:)  :: xd
       type(dualzn), allocatable, dimension(:) :: frd
     end function fvecdual
  end interface
  
  interface d2fscalar
     module procedure d2fscalarvv
     module procedure d2fscalaruv       
  end interface d2fscalar

contains
  !Hessian operator
  function Hessian(fsd,qcmplx) result(Hmat)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: qcmplx    
    complex(prec), dimension(size(qcmplx),size(qcmplx)) :: Hmat
    complex(prec), dimension(size(qcmplx)) :: ei, ej
    complex(prec) :: hij
    integer :: i,j,m
    
    m = size(qcmplx)
    !diagonal components
    ei = 0.0_prec
    do i = 1, m
       ei(i) = 1.0_prec
       Hmat(i,i) = d2fscalar(fsd,ei,qcmplx)
       ei(i) = 0.0_prec 
    end do
    
    !off-diagonal components
    do i=1, m
       ei(i) = 1.0_prec
       ej = 0.0_prec
       do j = i+1, m
          ej(j-1) = 0.0_prec
          ej(j) = 1.0_prec
          hij = 0.5_prec * (d2fscalar(fsd, ei + ej, qcmplx) -    &
               Hmat(i,i) - Hmat(j,j))  

          Hmat(i,j) = hij
          Hmat(j,i) = hij
       end do
       ei(i) = 0.0_prec
    end do
  end function Hessian
  
  !d2fscalaruv(f,u,v,q) gives the second-order directional derivative
  !of the scalar function f: D^m --> D, where f = f(q), along vectors u
  !and v, evaluated at point q
  function d2fscalaruv(fsd,x,y,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: x, y, q
    complex(prec) :: fr

    if(all(x == y)) then
       fr = d2fscalarvv(fsd,x,q)
    else
       fr = 0.5_prec*(d2fscalarvv(fsd,x + y,q) - d2fscalarvv(fsd,x,q) -&
            d2fscalarvv(fsd,y,q))
    end if
  end function d2fscalaruv
  
  !Second-order directional derivative of a scalar function along vector
  !v, evaluated at point q
  function d2fscalarvv(fsd,v,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: v, q
    complex(prec) :: fr
    type(dualzn) :: eps1
  
    eps1 = xto_dzn(0,2)
    eps1%f(1) = 1

    fr = f_part(fsd(q + eps1*v),2)
  end function d2fscalarvv
  
  !Jacobian operator: To optimize efficiency, we include the parameter
  !'n', representing the dimension of 'fvecd'.
  function Jacobian(fvecd,qcmplx,n) result(Jmat)
    procedure(fvecdual) :: fvecd
    complex(prec), intent(in), dimension(:) :: qcmplx
    integer, intent(in) :: n
    complex(prec), dimension(n,size(qcmplx)) :: Jmat
    complex(prec), dimension(size(qcmplx))   :: ei
    integer :: i

    do i = 1,size(qcmplx)
       ei = 0
       ei(i) = 1      
       Jmat(:,i) = d1fvector(fvecd,ei,qcmplx,n)
    end do

  end function Jacobian
  
  !First-order directional derivative of a vector function along vector
  !v, evaluated at point q. To optimize efficiency, we include the
  !parameter 'n', representing the dimension of 'fvecd'.
  !
  function d1fvector(fvecd,v,q,n) result(fr)
    procedure(fvecdual) :: fvecd
    complex(prec), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(prec), dimension(n) :: fr  
    type(dualzn) :: eps1

    eps1 = xto_dzn(0,1)
    eps1%f(1) = 1.0_prec

    fr = f_part(fvecd(q + eps1*v),1)
  end function d1fvector

  !gradient operator
  function gradient(fsd,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: q
    complex(prec), dimension(size(q)) :: fr
    complex(prec), dimension(size(q)) :: ei
    integer :: i

    do i = 1,size(q)
       ei = 0
       ei(i) = 1      
       fr(i) = d1fscalar(fsd,ei,q)
    end do
  end function gradient

  !First-order directional derivative of a scalar function along vector
  !v, evaluated at point q
  function d1fscalar(fsd,v,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: v, q
    complex(prec) :: fr  
    type(dualzn) :: eps1

    !to initialize eps to "[0,1]" 
    call initialize_dualzn(eps1,1) !order one
    eps1%f(1) = 1.0_prec

    fr = f_part(fsd(q + eps1*v),1)
  end function d1fscalar
end module diff_mod
