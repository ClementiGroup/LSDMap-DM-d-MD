!---------------------------------------------------------------!
! LSDMap v1.0 - Aug 01 2012 - Initial Release                   !
!	                                                        !
! Developed by						        !
!    Wenwei Zheng       wz7@rice.edu			        !
!    Mary Rohrdanz      mar3@rice.edu			        !
!    Cecilia Clementi   cecilia@rice.edu			!
!      							        !
! Please reference the papers below if you use LSDMap:	        !
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., ! 
!    J. Chem. Phys., 134, 124116, 2011			        !
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R.,     !
!    and Clementi, C., J. Phys. Chem. B, 115, 13065-13074, 2011 !
!---------------------------------------------------------------!

!-----------------------------------------------------------------------
! Serial RMSD and nearest neighboring map calculation.
! 
! Input
! 1. Name of the trajectory
! 2. start and end point ID
! 3. number of neighbors
! See "rmsd_neighbor.input_example" for an example of input.
!
! Output
! Pairwise RMSD matrix and nearest neighboring map in binary format.
! Pairwise RMSD matrix serves as the input of p_wlsdmap.
! Nearest neighboring map serves as the input of p_local_mds.
! 
! Please reference the papers below if you use LSDMap:
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., 
!    J. Chem. Phys., 134, 124116, 2011
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R., and Clementi, C.,
!    J. Phys. Chem. B, 115, 13065-13074, 2011
!
! WZ Jul 25 2012
!-----------------------------------------------------------------------


program s_rmsd_neighbor

use Qsort_Module
use iso_c_binding
use ftn_c

implicit none

integer, parameter :: real_kind=8
integer :: ns,ne
integer :: inlen
integer :: Npoints,dim,Natoms
integer :: idx,jdx
integer :: Nneigh
integer,allocatable :: idneigh(:)
real,allocatable :: traj(:,:)
real(real_kind),allocatable :: dist(:)
real(real_kind),allocatable :: xx(:,:),yy(:,:),rot(:),weight(:)

character(200) :: nn_traj, nn_output, nn_neigh
integer,parameter :: sizeofreal=4

! use quaterion characteristic polynomial to calculate RMSD
! J Comp. Chem. 31, 7, 1561; J Comp. Chem. 32, 1, 185
!real(real_kind) calcrmsdrotationalmatrix
!external calcrmsdrotationalmatrix

! read input parameters
read(5,*) nn_traj
read(5,*) ns,ne
read(5,*) Nneigh

print *,'LSDMap v1.0 - Aug 01 2012 - Initial Release'
call current_time('Program start...')

! write status
inlen=index(nn_traj,' ') - 1
open(10,file=nn_traj,status='old')
read(10,*) Npoints,dim
print *,'trajectory read from file ', nn_traj
print *,'number of points in dataset = ',Npoints
print *,'original dimension = ',dim
print *,'number of nearest neighbors = ',Nneigh
print *,'from ',ns,' to ',ne
Natoms=dim/3

! read trajectories
allocate(traj(Npoints*Natoms,3))
do idx=1,Npoints
   read(10,*) (traj((idx-1)*Natoms+jdx,1),traj((idx-1)*Natoms+jdx,2),traj((idx-1)*Natoms+jdx,3),&
        jdx=1,Natoms)
enddo
close(10)

! define output file name
write(nn_output,'(a,a,i7,a,i7)') nn_traj(1:inlen),'_rmsd_',ns+9000000,'_',ne+9000000

open(unit=51,file='rmsd/'//nn_output,form='unformatted',access='direct',status='replace',&
     recl=Npoints*sizeofreal)
! open(unit=51,file='rmsd/'//nn_output,status='replace')
if (Nneigh/=0) then
   write(nn_neigh,'(a,a,i7,a,i7)') nn_traj(1:inlen),'_neighbor_',ns+9000000,'_',ne+9000000
   open(unit=52,file='neighbor/'//nn_neigh,status='replace')
endif

! allocate matrix
allocate(xx(Natoms,3))
allocate(yy(Natoms,3))
allocate(dist(Npoints))
allocate(idneigh(Npoints))
allocate(rot(9))
allocate(weight(Natoms))
weight=1.
do idx=ns,ne
   xx=traj((idx-1)*Natoms+1:idx*Natoms,:)
   ! calculate the distance to the regarding point
   do jdx=1,Npoints
      yy=traj((jdx-1)*Natoms+1:jdx*Natoms,:)
      dist(jdx)=calcrmsdrotationalmatrix(Natoms,xx(:,1),xx(:,2),xx(:,3),yy(:,1),yy(:,2),yy(:,3),rot,weight)
   enddo
   ! output rmsd
   write(51,rec=idx-ns+1) (real(dist(jdx)),jdx=1,Npoints)
!   do jdx=1,Npoints
!      write(51,'(f10.6)',advance='no') dist(jdx)
!   end do
!   write(51,*)
   ! call quick sort
   if(Nneigh/=0) then
      idneigh=(/(jdx,jdx=1,Npoints,1)/)
      call Qsort(dist,idneigh)
      ! output nearest neighbor graph
      do jdx=1,Nneigh
         write(52,'(i7,1x)',advance='no') idneigh(jdx)
      enddo
      write(52,*)
   endif
enddo

close(51)
if(Nneigh/=0) close(52)

call current_time('Program end.')

end program s_rmsd_neighbor

subroutine current_time(text)
  character(8)::date1
  character(10)::time1
  character(5)::zone
  integer::values1(8)
  character(*)::text
  call date_and_time(date1,time1,zone,values1)
  write(6,'(a,a,a)') time1,' : ',text 
end subroutine current_time
