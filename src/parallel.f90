!---------------------------------------------------------------!                                                           
! Function Module for LSDMap code                               !
!   to keep track of varaibles related to MPI                   !
!                                                               !
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uinversity of Edinburgh                 !                                                           
!---------------------------------------------------------------!                                                           
! Input                                                         !
!                                                               ! 
! Output                                                        !
!                                                               !
!---------------------------------------------------------------! 
! Develop log                                                   !
!                                                               !
! September 2013 Initial Release                                !
!---------------------------------------------------------------!
 

module parallel
  
  implicit none
  
  ! General MPI variables
  integer :: size    ! Total number of processors
  integer :: rank    ! Value for each processor, 0 to size-1 
  integer :: ierr    ! MPI error
  integer :: comm    ! MPI communicator

  ! Variables for MPI-Gathers
  integer, allocatable :: displacements(:)  ! Entry 'i' specifies the displacement at which to place
                                            ! incoming data from process 'i'
  integer, allocatable :: counts(:)         ! containing the number of elements that are received
                                            ! from each process
  
end module parallel
