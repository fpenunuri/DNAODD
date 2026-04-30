!module with example of functions
module test_functions_mod
  use dualzn_mod
  use config_mod
  implicit none
  private

  public :: RCR_rD, HTM, rot_mat
  public :: fstest, fvectest, SinProbF

contains
  ! Position vector of the RCR robot manipulator (dualzn version)
  ! Although fr is known to be a 3-dimensional vector, it is declared as
  ! allocatable in order to be used within the context of the
  ! dnfscalar and dnfvector functions.
  ! q = [th, phi, s, beta]
  function RCR_rD(q) result(fr)
    type(dualzn), intent(in), dimension(:) :: q
    type(dualzn), allocatable, dimension(:) :: fr
    real(prec), parameter :: BC = 3.0_prec, CD = 2.0_prec !particular case
    type(dualzn) :: th, s, phi, beta, zd
    type(dualzn), dimension(4,4) :: T1, T2, T3, T4, T5
    type(dualzn), dimension(4,4) :: fr4
    type(dualzn), dimension(3)   :: e1, e3, zerov
    type(dualzn), dimension(3)   :: Dr2, Dr3, Dr5
    integer :: ord_loc

    ord_loc = q(1)%ord
    zd = xto_dzn(0.0_prec, ord_loc)
    zerov = xto_dzn([0.0_prec, 0.0_prec, 0.0_prec], ord_loc)
    e1    = xto_dzn([1.0_prec, 0.0_prec, 0.0_prec], ord_loc)
    e3    = xto_dzn([0.0_prec, 0.0_prec, 1.0_prec], ord_loc)

    th   = q(1)
    phi  = q(2)
    s    = q(3)
    beta = q(4)

    Dr2 = s*e1
    Dr3 = BC*e3
    Dr5 = CD*e1

    T1 = HTM(th,   e3, zerov)
    T2 = HTM(phi,  e1, Dr2)
    T3 = HTM(zd,   e3, Dr3)
    T4 = HTM(beta, e3, zerov)
    T5 = HTM(zd,   e1, Dr5)

    fr4 = matmul(T1,matmul(T2,matmul(T3,matmul(T4,T5))))   

    allocate(fr(3))
    fr = fr4(1:3,4)
  end function RCR_rD

  !Homogeneous transformation matrix
  function HTM(th,eje,Dr) result(fr)
    type(dualzn), intent(in) :: th
    type(dualzn), intent(in), dimension(3) :: eje
    type(dualzn), intent(in), dimension(3) :: Dr
    type(dualzn), dimension(4,4) :: fr

    fr(1:3,1:3) = rot_mat(th,eje)
    fr(1:3,4) = Dr
    fr(4,:) = xto_dzn([0,0,0,1],th%ord) 
  end function HTM

  !rotation matrix
  function rot_mat(th,eje) result(fr)
    type(dualzn), intent(in) :: th
    type(dualzn), intent(in), dimension(3) :: eje
    type(dualzn), dimension(3) :: ejeu
    type(dualzn) :: n1, n2, n3
    type(dualzn), dimension(3,3) :: fr

    !making 'eje' a unit vector
    ejeu = eje/sqrt(sum(eje*eje));

    n1 = ejeu(1); n2 = ejeu(2); n3 = ejeu(3)

    fr(1,:) = [1 + (n2**2 + n3**2)*(cos(th) - 1), n1*n2 - n1*n2*     &
         cos(th) - n3*sin(th), n1*n3 - n1*n3*cos(th) + n2*sin(th)]

    fr(2,:) = [n1*n2 - n1*n2*cos(th) + n3*sin(th), 1 + (n1**2 +      &
         n3**2)*(cos(th) - 1), n2*n3 - n2*n3*cos(th) - n1*sin(th)]

    fr(3,:) = [n1*n3 - n1*n3*cos(th) - n2*sin(th), n2*n3 - n2*n3*    &
         cos(th) + n1*sin(th), 1 + (n1**2 + n2**2)*(cos(th) - 1)]
  end function rot_mat

  !sine problem function
  function SinProbF(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r
    type(dualzn) :: fr
    real(prec), parameter :: A = 2.5_prec, B = 5.0_prec, z = 30.0_prec

    fr = -(A*product(sin(r - z)) + product(sin(B*(r - z)))) 
  end function SinProbF
  
  !Example of scalar function f = sin(x*y*z) + cos(x*y*z)
  function fstest(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r
    type(dualzn) :: fr
    type(dualzn) :: x,y,z

    x = r(1); y = r(2); z = r(3)
    fr = sin(x*y*z) + cos(x*y*z)
  end function fstest

  !Example of vector function f = [f1,f2,f3]
  !f = fvectest(r) is a function f:D^m ---> Dn 
  function fvectest(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r
    type(dualzn), allocatable, dimension(:) :: fr
    type(dualzn) :: f1,f2,f3
    type(dualzn) :: x,y,z,w

    x = r(1); y = r(2); z = r(3); w = r(4)

    f1 = sin(x*y*z*w)
    f2 = cos(x*y*z*w)*sqrt(w/y - x/z)
    f3 = sin(log(x*y*z*w))

    allocate(fr(3))
    fr = [f1,f2,f3]
  end function fvectest
end module test_functions_mod

