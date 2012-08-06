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

program p_rmsd_neighbor

use Qsort_Module
use iso_c_binding
use ftn_c

implicit none

include 'mpif.h'

integer,parameter :: real_kind=8
!valuables for MPI
integer :: comm, nproc, myid, ierr

integer :: ns,ne
integer :: inlen
integer :: Npoints,dim,Natoms
integer :: idx,jdx
integer :: nloc,nstart,nend
integer :: Nneigh,extra
integer,allocatable :: idneigh(:)
real :: tmp
real,allocatable :: traj(:,:)
real(real_kind),allocatable :: dist(:)
real(real_kind),allocatable :: xx(:,:),yy(:,:),rot(:),weight(:)
character(200) :: nn_traj, nn_output, nn_neigh
integer,parameter :: sizeofreal=4,sizeofinteger=4

! use quaterion characteristic polynomial to calculate RMSD
! J Comp. Chem. 31, 7, 1561; J Comp. Chem. 32, 1, 185
!real(real_kind) calcrmsdrotationalmatrix
!external calcrmsdrotationalmatrix

call MPI_init(ierr)
comm=MPI_COMM_WORLD
call MPI_comm_size(comm,nproc,ierr)
call MPI_comm_rank(comm,myid,ierr)
if(myid==0) print *,'LSDMap v1.0 - Aug 01 2012 - Initial Release'
if(myid==0) call current_time('Program Start...')

! get parameters from input file
if (myid==0) then
   read(5,*) nn_traj
   read(5,*) ns,ne
   read(5,*) Nneigh
endif
call MPI_BCAST(nn_traj,200,MPI_CHARACTER,0,comm,ierr)
call MPI_BCAST(ns,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(ne,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(Nneigh,1,MPI_INTEGER,0,comm,ierr)

! split data
nloc=int((ne-ns+1)/nproc)
!if(myid<mod(ne-ns+1,nproc)) nloc=nloc+1
!allocate(nallloc(nproc))
!call MPI_ALLGATHER(nloc,1,MPI_INTEGER,nallloc,1,MPI_INTEGER,comm,ierr)
!nend=ns-1
!do idx=1,myid+1
!   nend=nend+nallloc(idx)
!enddo
!nstart=nend-nloc+1
if(myid<mod(ne-ns+1,nproc)) then
   nloc=nloc+1
   extra=0
else
   extra=mod(ne-ns+1,nproc)
endif
nstart=ns+myid*nloc+extra
nend=nstart+nloc-1


! write status
inlen=index(nn_traj,' ') - 1
if(myid==0) then
   open(10,file=nn_traj,status='old')
   read(10,*) Npoints,dim
   print *,'trajectory read from file ', nn_traj
   print *,'number of points in dataset = ',Npoints
   print *,'original dimension = ',dim
   print *,'number of nearest neighbors = ',Nneigh
   print *,'from ',ns,' to ',ne
   print *,'number of cores=',nproc
   print *,'points per core=',nloc,'~',nloc+1
endif
call MPI_BCAST(Npoints,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(dim,1,MPI_INTEGER,0,comm,ierr)
Natoms=dim/3

! read input trajectory
allocate(traj(Npoints*Natoms,3))
if (myid==0) then
   do idx=1,Npoints
      read(10,*) (traj((idx-1)*Natoms+jdx,1),traj((idx-1)*Natoms+jdx,2),traj((idx-1)*Natoms+jdx,3),&
           jdx=1,Natoms)
   enddo
   close(10)
endif
call MPI_BCAST(traj,Npoints*dim,MPI_REAL,0,comm,ierr)

! define output file name
write(nn_output,'(a,a,i7,a,i7)') nn_traj(1:inlen),'_rmsd_',nstart+9000000,'_',nend+9000000

open(unit=51,file='rmsd/'//nn_output,form='unformatted',access='direct',status='replace',&
     recl=Npoints*sizeofreal)
if (Nneigh/=0) then
   write(nn_neigh,'(a,a,i7,a,i7)') nn_traj(1:inlen),'_neighbor_',nstart+9000000,'_',nend+9000000
   open(unit=52,file='neighbor/'//nn_neigh,form='unformatted',access='direct',status='replace', &
   recl=Nneigh*sizeofinteger)
endif

! allocate storage
allocate(xx(Natoms,3))
allocate(yy(Natoms,3))
allocate(dist(Npoints))
allocate(idneigh(Npoints))
allocate(rot(9))
allocate(weight(Natoms))

if(myid==0)  print *,'start calculating pairwise rmsd and nearest neighbor graph'
weight=1.
do idx=nstart,nend
   xx=traj((idx-1)*Natoms+1:idx*Natoms,:)
   ! calculate the distance to the regarding point
   do jdx=1,Npoints
      yy=traj((jdx-1)*Natoms+1:jdx*Natoms,:)
      dist(jdx)=calcrmsdrotationalmatrix(Natoms,xx(:,1),xx(:,2),xx(:,3),yy(:,1),yy(:,2),yy(:,3),rot,weight)
   enddo
   ! output RMSD
   write(51,rec=idx-nstart+1) (real(dist(jdx)),jdx=1,Npoints)
   ! call quick sort
   if(Nneigh/=0) then
      idneigh=(/(jdx,jdx=1,Npoints,1)/)
      call Qsort(dist,idneigh)
      ! output nearest neighbor graph
      write(52,rec=idx-nstart+1) (idneigh(jdx),jdx=1,Nneigh)
   endif
enddo
close(51)
if(Nneigh/=0) close(52)

if(myid==0) call current_time('Program end.')

call MPI_finalize(ierr)

end program p_rmsd_neighbor

subroutine current_time(text)
  character(8)::date1
  character(10)::time1
  character(5)::zone
  integer::values1(8)
  character(*)::text
  call date_and_time(date1,time1,zone,values1)
  write(6,'(a,a,a)') time1,' : ',text 
end subroutine current_time
