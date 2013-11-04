!---------------------------------------------------------------!                                                           
! Function Module for LSDMap code                               !
! Module to keep track of variables related to Chemistry        !
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
 
module data
  
  implicit none
  
  integer :: ns         !
  integer :: ne         !
  integer :: Npoints    !
  integer :: dim        !
  integer :: Nneigh     !
  integer :: Natoms     !
  
  ! Variables for data distribution across processors
  integer :: nloc       !
  integer :: nstart     !
  integer :: nend       !
  integer :: nlimit     !
  integer :: extra      !

  real, allocatable :: traj(:,:)
  ! Added second dim for arrays below for book-keeping
  integer, allocatable, save :: idneigh(:,:)
  real(kind=8), allocatable :: dist(:,:)
  real(kind=8), allocatable, save :: tmp_rmsd(:,:)
  real(kind=8), allocatable, save :: FullEpsArray(:,:)
  real(kind=8), allocatable :: EpsArray(:,:)
  
  ! Name of the xyz-input file
  character(200) :: nn_traj
  
  integer :: norder
  integer :: ncore
  integer :: dmds
  integer :: kmin
  integer :: dk
  integer :: neps
  real(kind=8) :: seps
  real(kind=8) :: deps
  
  
  character(180) :: NN_input_weight
  character(200) :: output_file
  integer :: status_dmap
  integer :: status_eps
  integer :: column_eps
  real(kind=8) :: cutoff
  real(kind=8) :: eps0
  
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
