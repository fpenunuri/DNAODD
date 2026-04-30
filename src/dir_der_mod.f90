! Arbitrary order directional derivatives
module dir_der_mod
  use config_mod 
  use dualzn_mod
  implicit none

  private
  public :: dnfscalar, dnfvector
  public :: dnfscalar_1vec, dnfvector_1vec
  public :: basis_vectors
  public :: SPD, TSMDD

  !Abstract interfaces for functions passed as arguments to other
  !functions.
  abstract interface
     !Interface for a scalar dual function f: D^m --> D (similar to
     !f: R^m --> R)
     function fsdual(xd) result(frsd)
       import :: dualzn
       type(dualzn), intent(in), dimension(:) :: xd
       type(dualzn) :: frsd
     end function fsdual

     !f: D^m --> D^n
     function fvecdual(xd) result(frd)
       import :: dualzn
       type(dualzn), intent(in), dimension(:) :: xd
       type(dualzn), allocatable, dimension(:) :: frd
     end function fvecdual
  end interface
contains


  ! Xsol = TSMDD(RHS_ED, xinit_cond, t_vec) computes the solution of a
  !        system of ODEs
  !
  ! RHS_ED      : Right-hand side of the ODE system,
  !               RHS_ED = Fvec([x1, x2, ..., xn, t])
  !               Note: Fvec may explicitly depend on time
  !
  ! xinit_cond  : Initial conditions
  !
  ! t_vec       : Vector of points where the solution is required
  !
  ! Xsol        : Matrix containing the solutions (does not include the
  !               independent variable t_vec)
  !
  ! size(xinit_cond) : Number of equations (dependent variables)
  !
  ! size(t_vec)      : Number of time points
  !
  ! Real-valued Taylor Series Method for solving ODEs up to 4th order
  function TSMDD(RHS_ED, xinit_cond, t_vec) result(Xsol)
    procedure(fvecdual) :: RHS_ED
    real(prec), intent(in), dimension(:) :: xinit_cond, t_vec
    real(prec), dimension(size(t_vec),size(xinit_cond)) :: Xsol

    !the dualzn numbers are implemented for complex components
    complex(prec), dimension(size(xinit_cond)+1) :: qi0p !vector
    !matrices
    complex(prec), dimension(1,size(xinit_cond)+1) :: qi1p, qi2p, qi3p
    complex(prec), dimension(2,size(xinit_cond)+1) :: Qi12p
    !
    real(prec), dimension(size(xinit_cond)) :: x1p, x2p, x3p, x4p
    real(prec), dimension(size(t_vec)-1) :: hvec
    real(prec) :: hi, ti
    integer :: ndpv, npts, ninterv, i

    ndpv = size(xinit_cond) !number of dependent variables or equations
    npts = size(t_vec)
    ninterv = npts - 1
    Xsol(1,:) = xinit_cond
    hvec = t_vec(2:npts) - t_vec(1:ninterv)

    qi1p(1,1:ndpv) = 0.0_prec 
    qi1p(1,ndpv+1) = 1.0_prec
    qi2p = 0.0_prec
    qi3p = 0.0_prec

    do i=1, ninterv
       ti = t_vec(i)
       hi = hvec(i)

       qi0p(1:ndpv) =  Xsol(i,:)
       qi0p(ndpv+1) = ti

       x1p = real(f_part(RHS_ED(xto_dzn(qi0p,0)),0),kind=prec)

       qi1p(1,1:ndpv) = x1p 

       x2p = real(dnfvector(RHS_ED,qi1p,qi0p,ndpv,1),kind=prec)
       qi2p(1,1:ndpv) = x2p 

       x3p = real(dnfvector(RHS_ED,qi1p,qi0p,ndpv,2) +  &
            dnfvector(RHS_ED,qi2p,qi0p,ndpv,1),kind=prec)

       qi3p(1,1:ndpv) = x3p 

       Qi12p(1,:) = qi1p(1,:)
       Qi12p(2,:) = qi2p(1,:)

       x4p = real(dnfvector(RHS_ED,qi1p,qi0p,ndpv,3) +      &
            3.0_prec*dnfvector(RHS_ED,Qi12p,qi0p,ndpv,2) +  &
            dnfvector(RHS_ED,qi3p,qi0p,ndpv,1),kind=prec)

       Xsol(i+1,:) =  Xsol(i,:) + hi*x1p + hi**2 * x2p/2.0_prec +  &
            hi**3 * x3p/6.0_prec + hi**4 * x4p/24.0_prec
    end do
  end function TSMDD

  ! SPD(fsd, Aindx, q) returns the partial derivative (specified by the
  ! index matrix Aindx) of the function fsd, evaluated at point q.
  !
  ! The first column of Aindx specifies the variable indices with respect
  ! to which differentiation is performed.
  ! The second column of Aindx specifies how many times the derivative
  ! is taken with respect to each corresponding variable.
  !
  ! For example:
  !  Aindx = [[1,2], [2,1], [3,2]]  corresponds to  {{x,2}, {y,1}, {z,2}},
  !  meaning: ∂⁵f / ∂x² ∂y ∂z²
  !
  ! q is the point at which the derivative is evaluated.
  function SPD(fsd,Aindx,q) result(dnfs)
    procedure(fsdual) :: fsd
    integer, intent(in), dimension(:,:) :: Aindx
    complex(prec), intent(in), dimension(:) :: q
    complex(prec) :: dnfs
    integer, allocatable, dimension(:) :: auxvec
    complex(prec), allocatable, dimension(:,:) :: StndrBasVec

    integer :: nd, mq, rowsAindx

    mq = size(q)
    rowsAindx = size(Aindx,1)
    allocate(auxvec(rowsAindx))
    auxvec = Aindx(:,2)
    nd = sum(auxvec)
    deallocate(auxvec)

    StndrBasVec = basis_vectors(Aindx,mq)
    dnfs = dnfscalar(fsd,StndrBasVec,q,nd)
  end function SPD

  !---------------------------------------------------------------------
  ! vec_dir_der = dnfvector(fvecd, A, q, dimf, nd) computes the vector
  ! directional derivative.
  !
  ! fvecd: a vector-to-vector function over dual numbers, with the
  !        following interface:
  !
  ! interface
  !    function fvecdual(xd) result(frd)
  !      use dualzn_mod  ! <-- remove if already declared in the
  !                      !     caller's scope
  !      type(dualzn), intent(in), dimension(:)  :: xd
  !      type(dualzn), allocatable, dimension(:) :: frd
  !    end function fvecdual
  ! end interface
  !
  ! A  : matrix whose rows are the direction vectors along which 
  !      the directional derivative is computed (these vectors are not
  !      necessarily unit vectors)
  ! q  : the evaluation point
  ! dimf: the dimension of fvecd
  ! nd : the order of the directional derivative
  !
  ! Important note:
  ! If, for example, the 5th-order vector directional derivative is to
  ! be computed along a single direction, do not use a matrix A with
  ! five identical rows (this would be inefficient). Instead, use a
  ! matrix with a single row or call the dnfvector_1vec function.  
  function dnfvector(fvecd,A,q,dimf,nd) result(fr)
    procedure(fvecdual) :: fvecd     
    complex(prec), intent(in), dimension(:,:) :: A
    complex(prec), intent(in), dimension(:) :: q
    integer, intent(in) :: dimf, nd
    complex(prec), dimension(dimf) :: fr

    integer, allocatable, dimension(:,:) :: index_mat
    integer :: nrowsA, ncolsA, dimq
    integer :: k,  i, signo, filas_ixmat
    complex(prec), dimension(dimf) :: suma
    complex(prec), dimension(size(q)) :: vaux

    nrowsA = size(A,1)
    ncolsA = size(A,2)
    dimq   = size(q)

    if (nrowsA == 1) then
       fr = dnfvector_1vec(fvecd,A(1,:),q,dimf,nd)
       return
    end if

    if (dimq /= ncolsA) then
       write(*,'(A,I0,A,I0,A)') "Error (dnfvector): dimq (", dimq,   &
            ") must be equal to ncolsA (", ncolsA, ")"
       stop 1
    else if (nrowsA > 1 .and. nrowsA /= nd) then
       write(*,'(A,I0,A,I0,A)') "Error (dnfvector): the order of " //&
            "directional derivative nd (", nd, ") must be equal to"//&
            " nrowsA (", nrowsA, ") when more than one direction"  //&
            " vector is provided"
       stop 1
    end if

    fr = 0.0_prec
    signo = 1
    do k = nrowsA, 1, -1
       filas_ixmat = comb(nrowsA, k)
       allocate(index_mat(filas_ixmat, k))
       call generar_combinaciones(nrowsA, k, index_mat)

       suma = 0.0_prec
       do i=1, filas_ixmat
          vaux = sumarAvecIndx(A,index_mat(i,:))
          suma = suma + dnfvector_1vec(fvecd,vaux,q,dimf,nd)
       end do
       deallocate(index_mat)
       suma = signo*suma
       signo = -signo
       fr = fr + suma
    end do
    fr = fr/gamma(nrowsA + 1d0)
  end function dnfvector

  !---------------------------------------------------------------------
  !
  ! scalar_dir_der = dnfscalar(fsd,A,q,nd) computes the directional
  ! derivative.
  !
  ! fsd: a scalar-valued function over dual numbers (D^m ---> D) , with
  ! the following interface:
  !
  ! interface
  !    function fsdual(xd) result(frsd)
  !      use dualzn_mod  ! <-- remove if already declared in the
  !                      !     caller's scope
  !      type(dualzn), intent(in), dimension(:) :: xd
  !      type(dualzn) :: frsd
  !    end function fsdual
  ! end interface
  !
  ! A  : matrix whose rows are the direction vectors along which
  !      the directional derivative is computed (these vectors are not
  !      necessarily unit vectors)
  ! q  : the evaluation point
  ! nd : the order of the directional derivative
  !
  ! Important note:
  ! If, for example, the 5th-order directional derivative is to be
  ! computed along a single direction, do not use a matrix A with five
  ! identical rows (this would be inefficient). Instead, use a matrix
  ! with a single row or call the dnfscalar_1vec function.  
  function dnfscalar(fsd,A,q,nd) result(fr)
    procedure(fsdual) :: fsd      
    complex(prec), intent(in), dimension(:,:) :: A
    complex(prec), intent(in), dimension(:) :: q
    integer, intent(in) :: nd
    complex(prec) :: fr

    integer, allocatable, dimension(:,:) :: index_mat
    integer :: nrowsA, ncolsA, dimq
    integer :: k,  i, signo, filas_ixmat
    complex(prec) :: suma
    complex(prec), dimension(size(q)) :: vaux

    nrowsA = size(A,1)
    ncolsA = size(A,2)
    dimq   = size(q)

    if (nrowsA == 1) then
       fr = dnfscalar_1vec(fsd,A(1,:),q,nd)
       return
    end if

    if (dimq /= ncolsA) then
       write(*,'(A,I0,A,I0,A)') "Error (dnfscalar): dimq (", dimq,   &
            ") must be equal to ncolsA (", ncolsA, ")"
       stop 1
    else if (nrowsA > 1 .and. nrowsA /= nd) then
       write(*,'(A,I0,A,I0,A)') "Error (dnfscalar): the order of " //&
            "directional derivative nd (", nd, ") must be equal to"//&
            " nrowsA (", nrowsA, ") when more than one direction"  //&
            " vector is provided"
       stop 1
    end if

    fr = 0.0_prec
    signo = 1
    do k = nrowsA, 1, -1
       filas_ixmat = comb(nrowsA, k)
       allocate(index_mat(filas_ixmat, k))
       call generar_combinaciones(nrowsA, k, index_mat)

       suma = 0.0_prec
       do i=1, filas_ixmat
          vaux = sumarAvecIndx(A,index_mat(i,:))
          suma = suma + dnfscalar_1vec(fsd,vaux,q,nd)
       end do
       deallocate(index_mat)
       suma = signo*suma
       signo = -signo
       fr = fr + suma
    end do
    fr = fr/gamma(nrowsA + 1d0)
  end function dnfscalar

  !---------------------------------------------------------------------
  ! sum those rows of A indicated in indxv
  function sumarAvecIndx(A,indxv) result(fr)
    complex(prec), intent(in), dimension(:,:) :: A
    integer, intent(in), dimension(:) :: indxv
    complex(prec), dimension(size(A,2)) :: fr
    integer :: i, ni

    ni = size(indxv)
    fr = 0.0_prec
    do i=1,ni
       fr = fr + A(indxv(i),:)
    end do
  end function sumarAvecIndx

  !---------------------------------------------------------------------
  ! generate subsets of size k taken from {1,2,...,n}
  subroutine generar_combinaciones(n, k, subsets_out)
    integer, intent(in) :: n, k
    integer, intent(out), allocatable, dimension(:,:) :: subsets_out
    integer, dimension(k) :: subset
    integer :: i, j, fila, total

    total = comb(n, k)
    allocate(subsets_out(total, k))

    do i = 1, k
       subset(i) = i
    end do

    fila = 1
    do
       subsets_out(fila, :) = subset

       i = k
       do while (i >= 1)
          if (subset(i) /= n - k + i) exit
          i = i - 1
       end do

       if (i == 0) exit 

       subset(i) = subset(i) + 1
       do j = i + 1, k
          subset(j) = subset(j - 1) + 1
       end do

       fila = fila + 1
    end do
  end subroutine generar_combinaciones

  !---------------------------------------------------------------------
  !n-order directional derivative of a vector function along vector
  !v, evaluated at point q. To optimize efficiency, we include the
  !parameter 'dimf', representing the dimension of 'fvecd'.
  !
  function dnfvector_1vec(fvecd,v,q,dimf,n) result(fr)
    procedure(fvecdual) :: fvecd
    complex(prec), intent(in), dimension(:) :: v, q
    integer, intent(in) :: dimf, n
    complex(prec), dimension(dimf) :: fr  
    type(dualzn) :: eps1

    eps1 = xto_dzn(0.0_prec,n)
    eps1%f(1) = 1.0_prec

    fr = f_part(fvecd(q + eps1*v),n)
  end function dnfvector_1vec

  !---------------------------------------------------------------------
  !n-order directional derivative of a scalar function along vector
  !v, evaluated at point q
  function dnfscalar_1vec(fsd,v,q,n) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(prec) :: fr  
    type(dualzn) :: eps1

    !call set_order(n)

    eps1 = xto_dzn(0.0_prec,n)
    eps1%f(1) = 1.0_prec

    fr = f_part(fsd(q + eps1*v),n)
  end function dnfscalar_1vec

  !binomial coefficient (as integer)
  function comb(m, n) result(binom)
    integer, intent(in) :: m, n
    integer :: binom
    integer :: j

    if (n < 0 .or. n > m) then
       binom = 0
       return
    endif

    if (n == 0 .or. n == m) then
       binom = 1
    else
       binom = 1
       do j = 1, n
          binom = binom*(m-j+1)/j
       end do
    endif
  end function comb

  !---------------------------------------------------------------------
  ! generates a matrix whose rows are some vector e_k from the
  ! standard basis of R^n.
  ! dimEV: the dimension of the vectors (i.e., the space dimension)
  ! indx_mat: a matrix whose first column specifies which basis vectors
  ! to generate, and whose second column specifies how many times
  ! each basis vector should be repeated
  ! Example: basis_vectors([[3,1],[1,2]], 3) generates the matrix
  ! A = [e_3, e_1, e_1] with e_1 and e_3 in R^3
  function basis_vectors(indx_mat,dimEV) result(fr)
    integer, intent(in), dimension(:,:) :: indx_mat
    integer, intent(in) :: dimEV
    real(prec), allocatable, dimension(:,:) :: fr
    integer :: i, rows_fr

    integer, allocatable, dimension(:) :: indx_vec

    rows_fr=sum(indx_mat(:,2))
    allocate(fr(rows_fr,dimEV))

    indx_vec = int_repeat(indx_mat(:,1), indx_mat(:,2))
    do i=1, rows_fr
       fr(i,:) = est_vec(indx_vec(i),dimEV)
    end do
    deallocate(indx_vec)
  end function basis_vectors

  !---------------------------------------------------------------------
  ! generates a vector e_k of dimension rn with all components zero,
  ! except for the k-th component, which is set to 1
  function est_vec(k,rn) result(fr)
    integer, intent(in) :: k, rn
    real(prec), dimension(rn) :: fr

    fr = 0.0_prec
    fr(k) = 1.0_prec
  end function est_vec

  !---------------------------------------------------------------------
  ! ivec = int_repeat(iv1, iv2), repeats the integers in vector iv1 as
  ! many times as specified by the corresponding entries in vector iv2.
  ! Example: int_repeat([3,4,1], [1,1,2]) = [3,4,1,1]
  function int_repeat(integers2repeat, nrepeat) result(fr)
    integer, intent(in), dimension(:) :: integers2repeat, nrepeat
    integer, allocatable, dimension(:) :: fr

    integer :: dimfr, cont, i, j

    dimfr = sum(nrepeat)
    allocate(fr(dimfr))

    cont = 0
    do i = 1, size(integers2repeat)
       do j=1, nrepeat(i)
          cont = cont+1
          fr(cont) = integers2repeat(i)
       end do
    end do
  end function int_repeat
end module dir_der_mod

