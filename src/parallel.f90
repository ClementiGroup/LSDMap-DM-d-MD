!--------------------------------------------------------------------!
! E. Breitmoser, EPCC, University of Edinburgh                       !
! Module to keep track of varaibles related to MPI                   !
!--------------------------------------------------------------------!

module parallel

implicit none

integer :: size, rank, ierr, comm
integer, allocatable :: displacements(:),counts(:)

end module parallel
