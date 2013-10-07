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
!---------------------------------------------------------------!                                                           
! LSDMap v1.1 - Sept 2013 - Merger Release                      !                                                           
!                                                               !                                                          
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uonversity of Edinburgh                 !                                                          !---------------------------------------------------------------!                                                           

!-----------------------------------------------------------------------
! Parallel RMSD and nearest neighboring map calculation.
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

subroutine p_rmsd_neighbor

use Qsort_Module
use iso_c_binding
use ftn_c
use parallel, only : rank, ierr, comm, displacements, counts
use data, only : Npoints,Nneigh,Natoms,nloc,nstart,nend,traj,idneigh,dist,tmp_rmsd,nn_traj,current_time


implicit none

include 'mpif.h'

integer :: inlen
integer :: idx,jdx,reci
real :: tmp
integer, allocatable :: localidneigh(:)
real(kind=8),allocatable :: tmpdist(:)
real(kind=8),allocatable :: xx(:,:),yy(:,:),rot(:),weight(:)
integer,parameter :: sizeofreal=4,sizeofinteger=4
character(200) :: nn_output, nn_neigh

! use quaterion characteristic polynomial to calculate RMSD
! J Comp. Chem. 31, 7, 1561; J Comp. Chem. 32, 1, 185
!real(real_kind) calcrmsdrotationalmatrix
!external calcrmsdrotationalmatrix

if(rank==0) call current_time('Subroutine p_rmsd_neighbor start...')

! allocate storage
allocate(xx(3,Natoms))
allocate(yy(3,Natoms))
allocate(localidneigh(Npoints))
allocate(rot(9))
allocate(weight(Natoms))
allocate(tmpdist(Npoints))

inlen=index(nn_traj,' ') - 1


if(rank==0)  print *,'start calculating pairwise rmsd and nearest neighbor graph'
weight=1.


do idx=nstart,nend
!   xx=traj((idx-1)*Natoms+1:idx*Natoms,:)
   xx=traj(:,(idx-1)*Natoms+1:idx*Natoms)
   ! calculate the distance to the regarding point
   do jdx=1,Npoints
!      yy=traj((jdx-1)*Natoms+1:jdx*Natoms,:)
      yy=traj(:,(jdx-1)*Natoms+1:jdx*Natoms)
!      dist(jdx,idx)=calcrmsdrotationalmatrix(Natoms,xx(:,1),xx(:,2),xx(:,3),yy(:,1),yy(:,2),yy(:,3),rot,weight)
      dist(jdx,idx)=calcrmsdrotationalmatrix(Natoms,xx(1,:),xx(2,:),xx(3,:),yy(1,:),yy(2,:),yy(3,:),rot,weight)
   enddo
   ! output RMSD 
   tmpdist(:) = dist(:,idx)
   ! call quick sort
   if(Nneigh/=0) then
      localidneigh=(/(jdx,jdx=1,Npoints)/)
      call Qsort(tmpdist,localidneigh)  ! dist is inout!!!
      reci = idx-nstart+1
      do jdx=1,Nneigh
         idneigh(jdx,reci) = localidneigh(jdx)
      enddo
   endif
enddo

!To be used in p_wlsdmap.f90
call MPI_ALLGATHERV(dist,Npoints*nloc,MPI_REAL8,tmp_rmsd,Npoints*counts,displacements,MPI_REAL8,comm,ierr)


if(rank==0) call current_time('Subroutine p_rmsd_neighbor end.')

deallocate(localidneigh)
deallocate(xx)
deallocate(yy)
deallocate(rot)
deallocate(weight)
deallocate(tmpdist)

end subroutine p_rmsd_neighbor

