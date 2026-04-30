! main program
program main
  use config_mod
  use dualzn_mod
  use dir_der_mod
  use test_functions_mod
  implicit none

  integer, allocatable, dimension(:,:) :: indx_mat
  complex(prec), parameter :: ii = (0.0_prec,1.0_prec)
  complex(prec), dimension(3) :: q
  complex(prec) :: dnfs

  q = [0.1_prec + ii, 0.2_prec + ii, 0.3_prec + ii]

  allocate(indx_mat(3,2))
  ![[1,2],[2,1],[3,2]] <---> {{x,2},{y,1},{z,2}}
  indx_mat(:,1)=[1,2,3]
  indx_mat(:,2)=[2,1,2]

  write(*,"(A)") "--SPD--"
  dnfs = SPD(fstest,indx_mat,q)
  write(*,*) "d5f/dx2dydz2:", dnfs
  deallocate(indx_mat)
end program main

