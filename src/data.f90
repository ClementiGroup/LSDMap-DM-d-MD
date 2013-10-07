!--------------------------------------------------------------------!                                 
! E. Breitmoser, EPCC, University of Edinburgh                       !                                 
! Module to keep track of variables related to Chemistry             !
!--------------------------------------------------------------------!                                 
module data

implicit none

integer :: ns, ne
integer :: Npoints, dim, Nneigh, Natoms

integer :: nloc,nstart,nend,nlimit
integer :: extra
real, allocatable :: tmp2(:,:),tmpTraj(:,:),traj(:,:)
! Added second dim for arrays below for book-keeping
integer, allocatable, save :: idneigh(:,:)
real(kind=8), allocatable :: dist(:,:)
real(kind=8), allocatable, save :: tmp_rmsd(:,:)
real(kind=8), allocatable, save :: FullEpsArray(:,:)
real(kind=8), allocatable :: EpsArray(:,:)

! Name of the xyz-input file
character(200) :: nn_traj

integer :: norder
integer :: ncore, dmds, kmin,dk,neps
real(kind=8) :: seps,deps


character(180) :: NN_input_weight
character(200) :: output_file
integer :: status_dmap,status_eps,column_eps
real(kind=8) :: cutoff, eps0

contains
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
