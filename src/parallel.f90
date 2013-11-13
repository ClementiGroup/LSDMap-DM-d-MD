!---------------------------------------------------------------!                                                           
! Function Module for LSDMap code                               !
!   to keep track of varaibles related to MPI                   !
!                                                               !
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uinversity of Edinburgh                 ! 
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
  integer, parameter :: MPI_success =0 ! to check MPI error success/failure. Set default to failure

  ! Variables for MPI-Gathers
  integer, allocatable :: displacements(:)  ! Entry 'i' specifies the displacement at which to place
                                            ! incoming data from process 'i'
  integer, allocatable :: counts(:)         ! containing the number of elements that are received
                                            ! from each process

contains

  ! Subroutine to stop the code and produce error message if MPI call fails
  subroutine trace_exit(string,status)
    implicit none
    character(len=*), intent(in) :: string 
    integer, intent(in) :: status
    
    if (status.ne.MPI_success) then
       write(*,*) 'Error in MPI', string
       stop
    end if
    
  end subroutine trace_exit
  
  
end module parallel
