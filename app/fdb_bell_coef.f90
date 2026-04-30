! main program
program main
  use fdb_bell_mod
  implicit none
  integer :: n
  type(term_t), allocatable :: terms(:)

  write(*,'(a)') 'nd:'
  read(*,*) n

  terms = generate_fdb_bell(n)
  call print_terms(n, terms)
  deallocate(terms)
end program main

