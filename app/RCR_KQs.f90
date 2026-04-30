
!main program
program main
  use config_mod
  use dualzn_mod
  use dir_der_mod
  use test_functions_mod
  implicit none

  call run_subroutine()
  
contains
  ! Example of computing kinematic quantities. If desired, dedicated
  ! functions can be implemented; however, the code below is provided
  ! for illustrative purposes
  subroutine run_subroutine()
    real(prec), parameter :: PI_2 = 2.0_prec*atan(1.0_prec)
    integer, parameter :: loc_ord = 5
    complex(prec), dimension(4) :: q0p, q1p, q2p, q3p, q4p, q5p
    real(prec), dimension(3) :: rD1, rD2, rD3, rD4, rD5, aux
    complex(prec) :: dir1_mat(1,4), dir2_mat(2,4), dir3_mat(3,4)
    complex(prec), dimension(4,4) :: dir4_mat
    
    q0p = [PI_2,     0.0_prec, 2.0_prec, 0.0_prec]
    q1p = [1.0_prec, 5.0_prec, 1.0_prec, 1.0_prec]
    q2p = [1.0_prec, 0.0_prec, 2.0_prec, 1.0_prec]
    q3p = [1.0_prec, 2.0_prec, 3.0_prec, 4.0_prec]
    q4p = [4.0_prec, 5.0_prec, 6.0_prec, 7.0_prec]
    q5p = [1.0_prec, 3.0_prec, 5.0_prec, 7.0_prec]
    
    !velocity
    dir1_mat(1,:) = q1p
    rD1 = real(dnfvector(RCR_rD,dir1_mat,q0p,3,1),kind=prec)

    !acceleration
    aux = real(dnfvector(RCR_rD,dir1_mat,q0p,3,2),kind=prec)
    dir1_mat(1,:) = q2p    
    rD2 = aux + real(dnfvector(RCR_rD,dir1_mat,q0p,3,1),kind=prec)

    !jerk
    dir1_mat(1,:)  = q1p
    aux = real(dnfvector(RCR_rD,dir1_mat,q0p,3,3),kind=prec)
    
    dir2_mat(1,:) = q1p
    dir2_mat(2,:) = q2p

    dir1_mat(1,:)  = q3p
    
    rD3 = aux + 3*real(dnfvector(RCR_rD,dir2_mat,q0p,3,2),kind=prec) + &
         real(dnfvector(RCR_rD,dir1_mat,q0p,3,1),kind=prec)

    !jounce/snap
    dir3_mat(1,:) = q1p
    dir3_mat(2,:) = q1p
    dir3_mat(3,:) = q2p

    dir2_mat(1,:) = q1p
    dir2_mat(2,:) = q3p
    
    dir1_mat(1,:) = q1p

    aux = real(dnfvector(RCR_rD,dir1_mat,q0p,3,4),kind=prec)
    dir1_mat(1,:) = q2p
    aux =  aux + 3*real(dnfvector(RCR_rD,dir1_mat,q0p,3,2),kind=prec)
    dir1_mat(1,:) = q4p
    aux =  aux + real(dnfvector(RCR_rD,dir1_mat,q0p,3,1),kind=prec)

    rD4 = aux + 6*real(dnfvector(RCR_rD,dir3_mat,q0p,3,3),kind=prec) + &
         4*real(dnfvector(RCR_rD,dir2_mat,q0p,3,2),kind=prec)

    !D5r
    ! --- Faa di Bruno Bell Terms for n = 5 ---
    ! {[1,1],[[5,1]]}
    ! {[5,2],[[1,1],[4,1]]}
    ! {[10,2],[[2,1],[3,1]]}
    ! {[10,3],[[1,2],[3,1]]}
    ! {[15,3],[[1,1],[2,2]]}
    ! {[10,4],[[1,3],[2,1]]}
    ! {[1,5],[[1,5]]}
    !
    ! 1*d1(q5p) + 5*d2(q1p,q4p) + 10*d2(q2p,q3p) + 10*d3([q1p]^2,q3p) +
    ! 15*d3(q1p,[q2p]^2) + 10*d4([q1p]^3,q2p) + 1*d5([q1p]^5)
    !
    ! dn(x) := dn([x]^n)
    ! [x]^n := x,x,...,x; "n times"
    !
    dir4_mat(1,:) = q1p
    dir4_mat(2,:) = q1p
    dir4_mat(3,:) = q1p
    dir4_mat(4,:) = q2p

    dir3_mat(1,:) = q1p
    dir3_mat(2,:) = q1p
    dir3_mat(3,:) = q3p

    dir2_mat(1,:) = q1p
    dir2_mat(2,:) = q4p
    
    dir1_mat(1,:) = q5p
    aux = real(dnfvector(RCR_rD,dir1_mat,q0p,3,1),kind=prec)
    dir1_mat(1,:) = q1p
    
    aux = aux + real(dnfvector(RCR_rD,dir1_mat,q0p,3,5),kind=prec) +   &
         5*real(dnfvector(RCR_rD,dir2_mat,q0p,3,2),kind=prec)

    dir2_mat(1,:) = q2p
    dir2_mat(2,:) = q3p
    aux = aux + 10*real(dnfvector(RCR_rD,dir2_mat,q0p,3,2),kind=prec) +&
         10*real(dnfvector(RCR_rD,dir3_mat,q0p,3,3),kind=prec)

    dir3_mat(2,:) = q2p    
    dir3_mat(3,:) = q2p
    rD5 = aux + 15*real(dnfvector(RCR_rD,dir3_mat,q0p,3,3),kind=prec) +&
         10*real(dnfvector(RCR_rD,dir4_mat,q0p,3,4),kind=prec)     

    write(*,"(A)") NEW_LINE("A")
    write(*,"(A)") "Kinematic quantities (first- to fifth-order):"
    write(*,"(*(f9.2,1x))") rD1
    write(*,"(*(f9.2,1x))") rD2
    write(*,"(*(f9.2,1x))") rD3
    write(*,"(*(f9.2,1x))") rD4
    write(*,"(*(f9.2,1x))") rD5
  end subroutine run_subroutine  
end program main
