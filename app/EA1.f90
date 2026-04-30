!main program
program EA1
  use config_mod
  use dualzn_mod
  use dir_der_mod
  use test_functions_mod
  implicit none

  integer, parameter :: m = 3000, nd = 7
  complex(prec), allocatable, dimension(:) :: q, xvec
  complex(prec), allocatable, dimension(:,:) :: xmat
  complex(prec) :: dnfs
  integer :: i, k, mvar, tl
  real(prec) :: ir
  integer, dimension(4) :: mvec = [100, 1000, 2000, 3000]
  real(8) :: t1, t2, dt

  if (max_order_dualzn < nd) then
     write(*,"(A,I0)") "use:  fpm --flag -DMAX_ORDER_DUALZN=",nd
  end if
  
  write(*,*)
  write(*,"(A,1x,I0)") "directional derivatives of order:", nd
  do k = 1, 4
     mvar = mvec(k)
     allocate(q(mvar))
     allocate(xvec(mvar))
     allocate(xmat(1,mvar))
     
     do i=1,mvar
        ir = i
        xvec(i) = sin(ir)
        q(i) = 1/ir
     end do

     xmat(1,:) = xvec

     call cpu_time(t1)
     dnfs = dnfscalar(SinProbF,xmat,q,nd)
     call cpu_time(t2)

     write(*,"(F0.5,A,F0.5,A,I0)") real(dnfs),", time(s): ",t2-t1, &
          ", m: ",mvar
     deallocate(q,xvec,xmat)
  end do
end program EA1

