!---------------------------------------------------------------!                                                           
! Function Module for LSDMap code                               !
! Module to keep track of variables related to Chemistry        !
!                                                               !
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uinversity of Edinburgh                 !                                                           
!---------------------------------------------------------------!
! Develop log                                                   !
!                                                               !
! September 2013 Initial Release                                !
!---------------------------------------------------------------!
 
module data
  
  implicit none
  
  integer :: ns         ! Start point ID
  integer :: ne         ! End  point ID
  integer :: Npoints    ! 
  integer :: dim        ! 
  integer :: Nneigh     ! Number of neighbors
  integer :: Natoms     ! := dim/3
  
  ! Variables for data distribution across processors
  integer :: nloc       ! 
  integer :: nstart     ! Local loop start index
  integer :: nend       ! Local loop end index
  integer :: nlimit     !
  integer :: extra      !

  real, allocatable :: traj(:,:)  ! Input trajectory
  ! Added second dim for arrays below for book-keeping between subroutine calls
  integer, allocatable, save :: idneigh(:,:)           !
  real(kind=8), allocatable :: dist(:,:)               !

  ! Introduced to be able to avoid intermediate IO between subroutine calls
  real(kind=8), allocatable, save :: tmp_rmsd(:,:)     ! Gathers all dist into one array on all processors
  real(kind=8), allocatable, save :: FullEpsArray(:,:) ! Gathers all EpsArrays into one array on processor 0
  real(kind=8), allocatable :: EpsArray(:,:)           ! Local to each processor, keeps track of 6 values,
                                                       ! which can be written to localscale/...
  
  character(200) :: nn_traj  ! Name of the xyz-input file  
  integer :: norder     ! Job ID
  integer :: ncore      ! Number of CPUs to share the trajectory ncore
  integer :: dmds       ! Number of points for MDS
  integer :: kmin       ! Start point for MDS spectra
  integer :: dk         ! Step size for MDS spectra
  integer :: neps       ! Cutoff for the first derivative of MDS spectra
  real(kind=8) :: seps  !     "
  real(kind=8) :: deps  !     "
  
  
  character(180) :: NN_input_weight  ! Weight file name
  character(200) :: output_file      ! Output file name
  integer :: status_dmap             ! Status of LSDMap
  integer :: status_eps              !     "
  integer :: column_eps              ! Local scale file name or
  real(kind=8) :: eps0               ! constant local scale value
  real(kind=8) :: cutoff             ! Cutoff for sparse matrix
  
  logical :: write_rmsd        ! If .true. the intermediate IO files '/rmsd/...' are written
  logical :: write_neighbor    ! If .true. the intermediate IO files '/neighbor/...' are written
  logical :: write_localscale  ! If .true. the intermediate IO files '/localscale/...' are written
  
contains

  ! Subroutine to write date and time into output file 
  subroutine current_time(text)
    character(8)::date1
    character(10)::time1
    character(5)::zone
    integer::values1(8)
    character(*)::text
    call date_and_time(date1,time1,zone,values1)
    write(6,'(a,a,a)') time1,' : ',text
  end subroutine current_time
  
  
  
end module data
