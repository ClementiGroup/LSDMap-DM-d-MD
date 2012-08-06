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
! Parallel LSDMap
!
! Input
! 1. Status of LSDMap (local scale or constant, weighted or non-weighted)
! 2. Number of points
! 3. RMSD binary file name
! 4. Weight file
! 5. Output file name
! 6. Cutoff for sparse matrix
! 7. local scale file name or constant local scale value
! 8. file name of the original trajectory
! 9. file name of the new trajectory
! 10. number of points to add for each loop
! 11. weight file name for the new points
! See "wlsdmap_embed.input_example" for an example of input.
!
! Output
! Diffusion coordinates for the new points 
! 
! Please reference the papers below if you use LSDMap:
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., 
!    J. Chem. Phys., 134, 124116, 2011
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R., and Clementi, C.,
!    J. Phys. Chem. B, 115, 13065-13074, 2011
!
! WZ Jul 25 2012
!-----------------------------------------------------------------------

program wlsdmap_embed

use iso_c_binding
use ftn_c

implicit none

include 'mpif.h'
include 'debug.h'
include 'stat.h'

integer,parameter :: real_kind=8

!
! MPI valuable
!
integer :: comm, nproc, myid, ierr

integer :: NlocPoint
integer, allocatable :: NallPoint(:)
integer, allocatable :: idxlocPoint(:)
real(real_kind), allocatable :: Weight(:)
real(real_kind),allocatable :: LocalScale(:)
real(real_kind),allocatable :: EV(:,:),EVAll(:,:),EVOld(:,:)

! valuables for checking Pearson correlations
real(real_kind) :: avex,avey
integer :: tmpint
real(real_kind) :: tmpreal1
real(real_kind) :: tmpdouble
integer :: Npoints, dim, inlen, onlen, Natoms
! loop index
integer :: ndx,idx,jdx,kdx,iternum, sumi, count
! eps0: eps0 in nlocal=0 mode (input)
! cutoff: sparse matrix cutoff, the entry which is smaller than cutoff will be treated as zero
! dist: temporary entry value
real(real_kind) :: cutoff, dist, distt, eps0
! sumrow: the sum of row
! sumloc: the sum of row in a local core
! vv(:): one column of the basis in a local core
real(real_kind), allocatable :: dloc(:),sumrow(:),sumloc(:),vv(:),v_all(:),ww(:)
! tmp_rmsd: temporarily store the rmsd data from the input file
real,allocatable :: tmp_rmsd(:)
! sparse storage
integer, allocatable :: jdxloc(:),iloc(:)
! input and output file name
character(180) :: NN_input_rmsd, NN_input_eps, NN_input_weight, NN_input_loc,&
 NN_old_traj,NN_new_traj,NN_new_weight
character(200) :: output_file,output,NN_output_add
! PARACK's tol
real(real_kind), parameter ::  zero = 0.0D+0
! nspace: the ratio of the space needed for the sparse system
real(real_kind), parameter :: nspace=1.

integer :: status_dmap,status_eps,column_eps
real(real_kind),allocatable :: tmp_eps(:)

! value for add_points
integer :: Naddall,Nadd,iii,nloop,cid,Naddlast
integer :: NpointsNew,NlocPointOld
real,allocatable :: d_addcolumn(:,:)
real,allocatable :: traj_old(:,:)
real,allocatable :: traj_new(:,:)
real(real_kind),allocatable :: xx(:,:),yy(:,:),rot(:),weight_rmsd(:)

!
!     ARPACK valuable
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
integer          maxnloc, maxnev, maxncv, ldv
parameter       (maxnloc=20000, maxnev=10, maxncv=30,ldv=maxnloc )

!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
real(real_kind),allocatable :: v(:,:), workl(:),workd(:), d(:,:), resid(:),ax(:)

logical,allocatable :: select(:)
integer,allocatable :: iparam(:), ipntr(:)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
character        bmat*1, which*2
integer          ido, n, nev, ncv, lworkl, info, j, nconv, maxitr, mode, ishfts
logical          rvec
real(real_kind)             tol, sigma
!
!     %----------------------------------------------%
!     | Local Buffers needed for MPI communication |
!     %----------------------------------------------%
!
real(real_kind),allocatable ::  mv_buf(:)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
real(real_kind) pdnorm2 

external pdnorm2, daxpy
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic :: abs

!real(real_kind) calcrmsdrotationalmatrix
!external calcrmsdrotationalmatrix

call MPI_init(ierr)
comm=MPI_COMM_WORLD
call MPI_comm_size(comm, nproc, ierr)
call MPI_comm_rank(comm, myid ,ierr)

allocate(v(ldv,maxncv), workl(maxncv*(maxncv+8)),workd(3*maxnloc), d(maxncv,2), &
 resid(maxnloc),ax(maxnloc),mv_buf(maxnloc))
allocate(select(maxncv))
allocate(iparam(11), ipntr(11))

if (myid==0) write (6,*) 'LSDMap v1.0 - Aug 01 2012 - Initial Release'
if (myid==0) write (6,*) 'Program Start...'

nev=10

!!
!! Read Input Parameters
!!
if (myid==0) then
   read(5,*) status_dmap,status_eps
   read(5,*) Npoints
   read(5,*) NN_input_rmsd  
   if(status_dmap==1) then
      read(5,*) NN_input_weight
   else
      read(5,*)
   endif
   read(5,'(a)') output_file
   read(5,*) cutoff
   if(status_eps==1) then
      read(5,*) NN_input_eps,column_eps
   else
      read(5,*) eps0
   endif
   read(5,*) NN_old_traj
   read(5,*) NN_new_traj
   read(5,*) Nadd
   read(5,*) NN_new_weight
endif
call MPI_BCAST(Npoints,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(NN_input_rmsd,200,MPI_CHARACTER,0,comm,ierr)
call MPI_BCAST(output_file,200,MPI_CHARACTER,0,comm,ierr)
call MPI_BCAST(cutoff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_BCAST(status_dmap,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(status_eps,1,MPI_INTEGER,0,comm,ierr)
if(status_eps==0) call MPI_BCAST(eps0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

if (myid == 0) then
   call current_time()
   print *,'Program start...'
   if(status_dmap==1) then
      inlen = index(NN_input_weight,' ') - 1
      write(6,'(a,a)')'biased energy read from file ', NN_input_weight(1:inlen)
   endif
   inlen = index(NN_input_rmsd,' ') - 1
   write(6,'(2a)') 'RMSD read from file ', NN_input_rmsd(1:inlen)
   write(6,'(a,i12)') 'number of points in dataset = ', Npoints
   onlen = index(output_file,' ') - 1
   write(6,'(2a)') ' results output in file ', output_file(1:onlen)
   write(6,'(a,f6.3)') 'cutoff = ', cutoff
   if(status_eps==1) then
      inlen = index(NN_input_eps,' ') - 1
      write(6,'(a,a)') 'local scale read from file ',NN_input_eps(1:inlen)
   else
      write(6,'(a,f6.3)') 'use constant local scale', eps0
   endif
endif

if(myid==0)then
   open(50,file=NN_new_traj,status='old')
   read(50,*) Naddall,dim
   close(50)
   open(50,file=NN_old_traj,status='old')
   read(50,*) tmpint,dim
   if(tmpint/=Npoints) then
      print *,'Warning: The input value Npoints is not the same as the one in the trajectory file.'
      print *,'Program is terminated.'
      stop
   endif
   close(50)
   print *,'number of points to add:', Naddall
   print *,'number of points to add for one loop:',Nadd
endif
call MPI_BCAST(Naddall,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(Nadd,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(dim,1,MPI_INTEGER,0,comm,ierr)
Natoms=dim/3
nloop=ceiling(float(Naddall)/float(Nadd))
Naddlast=Naddall-(nloop-1)*Nadd
if(myid==0) then
!   if(Nadd*nloop/=Naddall) then
!      print *,'Warning: the number of points to add is not divisible by the number of points to add for one loop.'
!      print *,'Program is terminated.'
!      stop
!   endif
   print *,'number of loops:', nloop
   print *,'number of points in the last loop:',Naddlast
   print *,'dimension of the data:', dim
endif

! separate points into each cpu
NpointsNew=Npoints+Nadd
NlocPoint=int(NpointsNew/nproc)
if(myid<mod(NpointsNew,nproc)) then
   NlocPoint=NlocPoint+1
endif
allocate(NallPoint(nproc))
call MPI_ALLGATHER(NlocPoint,1,MPI_INTEGER,NallPoint,1,MPI_INTEGER,comm,ierr)
!call MPI_barrier(comm,ierr)

! build the index for NlocPoint
allocate(idxlocPoint(nproc+1))
idxlocPoint(1)=1
do idx=2,nproc+1
   idxlocPoint(idx)=idxlocPoint(idx-1)+NallPoint(idx-1)
enddo
if (idxlocPoint(myid+2)-1<=Npoints) then
NlocPointOld=NlocPoint
else if (idxlocPoint(myid+2)-1>Npoints+NlocPoint) then
   NlocPointOld=0
else
   NlocPointOld=NlocPoint-(idxlocPoint(myid+2)-1-Npoints)
endif

! generate the input file name
inlen = index(NN_input_rmsd,' ') - 1
if(idxlocPoint(myid+2)-1>Npoints) then
   write(NN_input_loc,'(a,a,i7,a,i7)') NN_input_rmsd(1:inlen),'_',&
   9000000+idxlocPoint(myid+1),'_',9000000+Npoints
else
   write(NN_input_loc,'(a,a,i7,a,i7)') NN_input_rmsd(1:inlen),'_',&
   9000000+idxlocPoint(myid+1),'_',9000000+idxlocPoint(myid+2)-1
endif

if(myid==0) then 
   write(NN_output_add,'(a,a)') output_file(1:onlen),'_embed.ev'
endif

! read local scale
allocate(LocalScale(NpointsNew))
if(myid==0) then
   if(status_eps==1)then
      call current_time()
      print *,'Read local scale...'
      allocate(tmp_eps(column_eps))
      open(20,file=NN_input_eps,status='old')
      do idx=1,Npoints
         read(20,*) (tmp_eps(jdx),jdx=1,column_eps)
         LocalScale(idx)=tmp_eps(column_eps)
      enddo
      close(20)
      deallocate(tmp_eps)
   endif
endif
if(status_eps==1) then 
   call MPI_BCAST(LocalScale(1:Npoints),Npoints,MPI_DOUBLE_PRECISION,0,comm,ierr)
else
   LocalScale=eps0
endif

! read weight
allocate(Weight(NpointsNew))
if(myid==0) then
   if(status_dmap==1)then
      call current_time()
      print *,'Read weight...'
      open(10,file=NN_input_weight,status='old')
      do idx=1,Npoints
         read(10,*) Weight(idx)
!         Weight(idx)=exp(Weight(idx)*beta)
      enddo
      close(10)
   endif
endif
if(status_dmap==1) then 
   call MPI_BCAST(Weight(1:Npoints),Npoints,MPI_DOUBLE_PRECISION,0,comm,ierr)
else
   Weight=1.
endif

! read old trajectory
if(myid==0) then
   call current_time()
   print *,'Read old trajectory...'
   open(10,file=NN_old_traj,status='old')
   read(10,*)
endif
allocate(traj_old(Npoints,dim))
if(myid==0) then
   do idx=1,Npoints
      read(10,*) (traj_old(idx,jdx),jdx=1,dim)
   enddo
   close(10)
endif
call MPI_BCAST(traj_old,Npoints*dim,MPI_REAL,0,comm,ierr)

!call MPI_barrier(comm,ierr)

!!
!! intialize values for the loop
!!
! CSR: iloc,jdxloc,dloc
!
if (myid==0) then
   call current_time()
   print *,'allocate block CSR storage...'
endif
allocate(iloc(NlocPoint+1))
sumi=ceiling(NpointsNew*nspace*NlocPoint)
allocate(dloc(sumi))
allocate(jdxloc(sumi))

if(myid==0)then
   open(30,file=NN_new_traj,status='old')
   read(30,*) 
endif
allocate(traj_new(Nadd,dim))
allocate(xx(Natoms,3))
allocate(yy(Natoms,3))
allocate(rot(9))
allocate(weight_rmsd(Natoms))
weight_rmsd=1.

if(myid==0) then
   open(200,file=NN_output_add,status='unknown')
   allocate(EVOld(NpointsNew,nev))
endif

!
! loop iii for adding points
!
do iii=1,nloop
   
   ! read new trajectory
   if(myid==0)then
      call current_time()
      print *,'start loop ',iii,' ...'
      if(iii==nloop) then
         do idx=1,Naddlast
            read(30,*) (traj_new(idx,jdx),jdx=1,dim)
         enddo
         do idx=Naddlast+1,Nadd
            traj_new(idx,:)=traj_old(idx-Naddlast,:)
         enddo
      else
         do idx=1,Nadd
            read(30,*) (traj_new(idx,jdx),jdx=1,dim)
         enddo
      endif
   endif
   call MPI_BCAST(traj_new,Nadd*dim,MPI_REAL,0,comm,ierr)

   ! read new weight
   if(myid==0) then
      if(status_dmap==1)then
         open(10,file=NN_new_weight,status='old')
         do idx=1,Nadd
            read(10,*) Weight(Npoints+idx)
         enddo
         close(10)
      endif
   endif
   if(status_dmap==1) then 
      call MPI_BCAST(Weight(Npoints+1:NpointsNew),Nadd,MPI_DOUBLE_PRECISION,0,comm,ierr)
   endif

   ! build diffusion kernel for the old points
   ! construct csr storage
   iloc=0
   count=0
   if(NlocPointOld>0) then
      open(10,file='rmsd/'//NN_input_loc,form='unformatted',status='old',&
      access='direct',recl=Npoints*4)
   endif
   allocate(tmp_rmsd(Npoints))
   do jdx=1,NlocPointOld
      ! read rmsd file
      read(10,rec=jdx) (tmp_rmsd(kdx),kdx=1,Npoints)
      iloc(jdx)=count+1
      ndx=idxlocPoint(myid+1)+jdx-1
      do kdx=1,Npoints
         dist=tmp_rmsd(kdx)**2
         distt=Weight(ndx)*Weight(kdx)*exp(-dist/Localscale(ndx)/LocalScale(kdx)/2.)
         if (distt >= cutoff) then
            count=count+1
            dloc(count)=distt
            jdxloc(count)=kdx
         endif
      enddo
      do kdx=Npoints+1,NpointsNew
         count=count+1
         jdxloc(count)=kdx
      enddo
   enddo
   do jdx=NlocPointOld+1,NlocPoint
      iloc(jdx)=count+1
      do kdx=1,NpointsNew
         count=count+1
         jdxloc(count)=kdx
      enddo
   enddo
   iloc(NlocPoint+1)=count+1
   if(NlocPointOld>0)then
      close(10)
   endif
   deallocate(tmp_rmsd)
   if(count>sumi) then
      print *,'Warning: core #',myid,' does not have enough memory to store the sparse matrix'
      print *,'core #',myid,' has used',real(count)/sumi*100, '% sparse storage'
      print *,'Please increase the value nspace in the program'
      stop
   endif
   print *,'core #',myid,' has used',real(count)/sumi*100, '% sparse storage'
   

   if (myid==0) then
      call current_time()
      print *,'rmsd new columns'
   endif
   ! calculate rmsd for the new columns
   allocate(d_addcolumn(Nadd,NlocPoint))
   do idx=1,Nadd
      do jdx=1,NlocPoint
!         xx=traj_new(idx,:)
         xx(:,1)=traj_new(idx,1:dim-2:3)
         xx(:,2)=traj_new(idx,2:dim-1:3)
         xx(:,3)=traj_new(idx,3:dim:3)
         if(jdx<=NlocPointOld) then
!            yy=traj_old(idxlocPoint(myid+1)+jdx-1,:)
            yy(:,1)=traj_old(idxlocPoint(myid+1)+jdx-1,1:dim-2:3)
            yy(:,2)=traj_old(idxlocPoint(myid+1)+jdx-1,2:dim-1:3)
            yy(:,3)=traj_old(idxlocPoint(myid+1)+jdx-1,3:dim:3)
         else
!            yy=traj_new(idxlocPoint(myid+1)-Npoints+jdx-1,:)
            yy(:,1)=traj_new(idxlocPoint(myid+1)-Npoints+jdx-1,1:dim-2:3)
            yy(:,2)=traj_new(idxlocPoint(myid+1)-Npoints+jdx-1,2:dim-1:3)
            yy(:,3)=traj_new(idxlocPoint(myid+1)-Npoints+jdx-1,3:dim:3)
         endif
!         call rmsd(dim,xx,yy,dist)
         dist=calcrmsdrotationalmatrix(Natoms,xx(:,1),xx(:,2),xx(:,3),yy(:,1),&
                 yy(:,2),yy(:,3),rot,weight_rmsd)
         d_addcolumn(idx,jdx)=dist
      enddo
   enddo
   
   if (myid==0) then
      call current_time()
      print *,'gather rmsd'
   endif
   ! gather rmsd matrix together at the cores with new points
   ! loop the number of new points here but the number of cores with new points because of the memory issue
   do idx=1,Nadd
      ! which core is the point in
      do jdx=1,nproc
         if(Npoints+idx<idxlocPoint(jdx+1)) exit
      enddo
      tmpint=jdx-1
      if(myid==tmpint)then
         allocate(tmp_rmsd(NpointsNew))
      endif
      call MPI_GATHERV(d_addcolumn(idx,:),NlocPoint,MPI_REAL,tmp_rmsd,NallPoint,&
      idxlocPoint(1:nproc)-1,MPI_REAL,tmpint,comm,ierr)
      if(myid==tmpint)then
         cid=Npoints+idx-(idxlocPoint(tmpint+1)-1)
         dloc(iloc(cid):iloc(cid+1)-1)=tmp_rmsd
         deallocate(tmp_rmsd)
      endif
   enddo
   call MPI_barrier(comm,ierr)

   if (myid==0) then
      call current_time()
      print *,'nearest neighbors'
   endif
   ! get the local scale for the new points
   ! find nearest neighbors
   if(status_eps==1)then
      do idx=NlocPointOld+1,NlocPoint
         tmpdouble=999.
         do jdx=1,Npoints
            if(dloc(iloc(idx)+jdx-1)<tmpdouble)then
               tmpint=jdx
               tmpdouble=dloc(iloc(idx)+jdx-1)
            endif
         enddo
         LocalScale(idxlocPoint(myid+1)-1+idx)=LocalScale(tmpint)
      enddo
      ! find which core starts having new points
      do jdx=1,nproc
         if(idxlocPoint(jdx+1)>Npoints) exit
      enddo
      ! share local scale and weight with the other cores
      do idx=jdx,nproc
         call MPI_BCAST(LocalScale(idxlocPoint(idx):idxlocPoint(idx+1)-1),&
         idxlocPoint(idx+1)-idxlocPoint(idx),MPI_DOUBLE_PRECISION,idx-1,comm,ierr)
      enddo
   endif

   if (myid==0) then
      call current_time()
      print *,'diffusion kernel for new points'
   endif
   ! construct the diffusion kernel of the new points
   ! new columns
   do idx=1,NlocPointOld
      kdx=idxlocPoint(myid+1)+idx-1
      do jdx=iloc(idx)+Npoints,iloc(idx+1)-1
         ndx=jdxloc(jdx)
         dloc(jdx)=Weight(kdx)*Weight(ndx)*exp(-d_addcolumn(ndx-Npoints,idx)**2&
         /LocalScale(kdx)/LocalScale(ndx)/2.)
      enddo
   enddo
   ! new rows
   do idx=NlocPointOld+1,NlocPoint
      kdx=idxlocPoint(myid+1)+idx-1
      do jdx=iloc(idx),iloc(idx+1)-1
         ndx=jdxloc(jdx)
         dloc(jdx)=Weight(kdx)*Weight(ndx)*exp(-dloc(jdx)**2/LocalScale(kdx)/LocalScale(ndx)/2.)
      enddo
   enddo
   
   deallocate(d_addcolumn)

   if (myid==0) then
      call current_time()
      print *,'start eig decomposition'
   endif
!
! start eigenvalue decomposition
!

! valueables for normalization
allocate(sumrow(NpointsNew))
allocate(sumloc(NlocPoint))
!!
!! sum over rows...
!!
sumloc=0.
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      sumloc(jdx)=sumloc(jdx)+dloc(kdx)
   enddo
enddo
!call MPI_barrier(comm, ierr)
call MPI_ALLGATHERV(sumloc,NlocPoint,MPI_DOUBLE_PRECISION,sumrow,NallPoint,&
  idxlocPoint(1:nproc)-1,MPI_DOUBLE_PRECISION,comm,ierr)
if (myid==0) then
   write(6,'(a)') 'all the cores finish sum over rows 1, sum has been gathered'
endif

!!
!! rescale matrix...
!!
if (myid==0) then
   write(6,'(a)') 'rescale matrix...'
endif
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      dloc(kdx)=dloc(kdx)*sqrt(Weight(idxlocPoint(myid+1)+jdx-1)*Weight(jdxloc(kdx))&
       /(sumrow(idxlocPoint(myid+1)+jdx-1)*sumrow(jdxloc(kdx))))
   enddo
enddo

!!
!! sum over rows...
!!
sumloc=0.
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      sumloc(jdx)=sumloc(jdx)+dloc(kdx)
   enddo
enddo
!call MPI_barrier(comm, ierr)
call MPI_ALLGATHERV(sumloc,NlocPoint,MPI_DOUBLE_PRECISION,sumrow,NallPoint,&
  idxlocPoint(1:nproc)-1,MPI_DOUBLE_PRECISION,comm,ierr)
if (myid==0) then
   write(6,'(a)') 'all the cores finish sum over rows 1, sum has been gathered'
endif
!!
!! rescale matrix...
!!
if (myid==0) then
   write(6,'(a)') 'rescale matrix...'
endif
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      dloc(kdx)=dloc(kdx)/sqrt(sumrow(idxlocPoint(myid+1)+jdx-1)*sumrow(jdxloc(kdx)))
   enddo
enddo
deallocate(sumloc)
!! sumrow is needed in basis calculation
sumrow=1/sqrt(sumrow)

call MPI_barrier(comm, ierr)

! valueables for matrix multiplication
allocate(vv(NlocPoint))
allocate(ww(NlocPoint))
allocate(v_all(NPointsNew))

!--------------------------------------------------------------
!       FROM HERE IT COMPUTES THE EIGENVALUES AND EIGENVECTORS 
!       OF THE DIFFUSION MATRIX 
!--------------------------------------------------------------

if (myid == 0) then
   write(6,'(a)') 'finished creating diffusion matrix... start diagonalization...'
endif

ndigit = -3
logfil = 6
msaupd = 1

! nev: number of eigenvalues
n=NpointsNew
!nev=10
ncv=30

!!-------------------------------------------------------------
!! can be deleted in f90 by using dynamic memory allocation
if ( NlocPoint .gt. maxnloc ) then
   print *, ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
   go to 9000
else if ( nev .gt. maxnev ) then
   print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
   go to 9000
else if ( ncv .gt. maxncv ) then
   print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
   go to 9000
end if
!!-------------------------------------------------------------

bmat = 'I'
which = 'LA'
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in PSSAUPD as       |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in PSSAUPD to start the Arnoldi        |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
lworkl = ncv*(ncv+8)
tol = zero 
info = 0
ido = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of PSSAUPD is used    |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in PSSAUPD.                         |
!     %---------------------------------------------------%
!
ishfts = 1
maxitr = 5000
mode   = 1
      
iparam(1) = ishfts 
iparam(3) = maxitr 
iparam(7) = mode

!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
iternum=0

 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine PSSAUPD and take| 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
iternum=iternum+1
call pdsaupd (comm, ido, bmat, NlocPoint, which, nev, tol, resid,ncv, v, ldv, &
  iparam, ipntr, workd, workl,lworkl, info )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
!
if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!

! get vector vv from the work vector workd(ipntr(1))
   vv(1:NlocPoint)=workd(ipntr(1):ipntr(1)+NlocPoint-1)
!   call MPI_barrier(comm, ierr)
!!
!! gather a column of vector from each core
!!
   call MPI_ALLGATHERV(vv, NlocPoint, MPI_DOUBLE_PRECISION, v_all, NallPoint, &
     idxlocPoint(1:nproc)-1, MPI_DOUBLE_PRECISION, comm,ierr)   
   call av (sumi, NlocPoint, NpointsNew,dloc, jdxloc, iloc, v_all, ww)
!   call MPI_barrier(comm,ierr)

! get the work block vector workd(ipntr(2)) from vector ww
   workd(ipntr(2):ipntr(2)+NlocPoint-1)=ww

!
!           %-----------------------------------------%
!           | L O O P   B A C K to call PSSAUPD again.|
!           %-----------------------------------------%
!
   go to 10
end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in PSSAUPD.|
!        %--------------------------%
!
   if ( myid .eq. 0 ) then
      print *, ' '
      print *, ' Error with _saupd, info = ', info
      print *, ' Check documentation in _saupd '
      print *, ' '
   endif
else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using PSSEUPD.               |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
   rvec = .true.
   if (myid==0) then
      write(6,'(a)') 'no fatal errors occurred in iteration'
      write(6,'(a)') 'start calling psseupd subroutine...'
   endif
   call pdseupd ( comm, rvec, 'All', select,d, v, ldv, sigma,bmat, NlocPoint,&
     which, nev, tol, resid, ncv, v, ldv,iparam, ipntr, workd, workl, lworkl, ierr )
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
   if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of PSSEUPD.|
!            %------------------------------------%
!
!
      if ( myid .eq. 0 ) then
         print *, ' '
         print *, ' Error with _seupd, info = ', ierr
         print *, ' Check the documentation of _seupd. '
         print *, ' '
      endif
   else
      nconv =  iparam(5)
      if (myid==0) then
         write(6,'(a)')'compute the residual norm...'
      endif
      do j=1, nconv
!      
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!


         !!
         !! gather a column of vector from each core
         !!
         call MPI_ALLGATHERV(v, NlocPoint, MPI_DOUBLE_PRECISION, v_all, NallPoint, &
           idxlocPoint(1:nproc)-1, MPI_DOUBLE_PRECISION, comm,ierr)   

         call av (sumi, NlocPoint, NpointsNew,dloc, jdxloc, iloc, v_all, ax)

         call daxpy(NlocPoint, -d(j,1), v(1,j), 1, ax, 1)
         d(j,2) = pdnorm2( comm, NlocPoint, ax, 1 )
      enddo



!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
      call pdmout(comm, 6, nconv, 2, d, maxncv, -6,'Ritz values and direct residuals')
   endif
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
   if (myid .eq. 0)then
      if ( info .eq. 1) then
         print *, ' '
         print *, ' Maximum number of iterations reached.'
         print *, ' '
      else if ( info .eq. 3) then
         print *, ' ' 
         print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
         print *, ' '
      end if
      print *, ' '
      print *, '_SDRV1 '
      print *, '====== '
      print *, ' '
      print *, ' Size of the matrix is ', n
      print *, ' The number of processors is ', nproc
      print *, ' The number of Ritz values requested is ', nev
      print *, ' The number of Arnoldi vectors generated',' (NCV) is ', ncv
      print *, ' What portion of the spectrum: ', which
      print *, ' The number of converged Ritz values is ', nconv 
      print *, ' The number of Implicit Arnoldi update',' iterations taken is ', iparam(3)
      print *, ' The number of OP*x is ', iparam(9)
      print *, ' The convergence criterion is ', tol
      print *, ' '
   endif
endif
!
!     %---------------------------%
!     | Done with program pssdrv1.|
!     %---------------------------%
!
9000 continue

!--------------------------------------------------------------
!	FROM HERE IT USES THE EIGENVALUES AND EIGENVECTORS 
!	OF DIFFUSION MATRIX TO OBTAIN THE DIFFUSION WAVELETS
!--------------------------------------------------------------
call MPI_barrier(comm, ierr)

nconv =  iparam(5)

if (myid==0) then
   write(6,'(a)') 'finish diagonalization... start postprocessing....'
endif

! store the eigenvectors
allocate(EV(NlocPoint,nev))

do jdx=1,nconv
   do idx=1,NlocPoint
      EV(idx,jdx)=sumrow(idxlocPoint(myid+1)-1+idx)*v(idx,jdx)
   enddo
enddo

if (myid==0) then
   allocate(EVAll(NpointsNew,nev))
endif

do jdx=1,nconv
   call MPI_GATHERV(EV(:,jdx),NlocPoint,MPI_DOUBLE_PRECISION,EVAll(:,jdx),NallPoint,&
     idxlocPoint(1:nproc)-1,MPI_DOUBLE_PRECISION,0,comm,ierr)
enddo

! normalization
if(myid==0) then
   do jdx=1,nconv
      tmpreal1=0.
      do kdx=1,Npoints
         tmpreal1=tmpreal1+EVAll(kdx,jdx)**2
      enddo
      tmpreal1=sqrt(tmpreal1)
      do kdx=1,NpointsNew
         EVAll(kdx,jdx)=EVAll(kdx,jdx)/tmpreal1
      enddo
      if(iii==1) then
         if(EVAll(1,jdx)<0) then
            EVAll(:,jdx)=-EVAll(:,jdx)
         endif
      else 
         if(jdx==nconv) then
            if(EVAll(1,nconv)<0) then
               EVAll(:,nconv)=-EVAll(:,nconv)
            endif
         else
            ! check Pearson correlation
            avex=0.0
            avey=0.0
            do kdx=1,Npoints
               avex=avex+EVOld(kdx,jdx)
               avey=avey+EVAll(kdx,jdx)
            enddo
            avex=avex/Npoints
            avey=avey/Npoints
            tmpreal1=0.0
            do kdx=1,Npoints
               tmpreal1=tmpreal1+(EVOld(kdx,jdx)-avex)*(EVAll(kdx,jdx)-avey)
            enddo
            if(tmpreal1<0) then
               EVAll(:,jdx)=-EVAll(:,jdx)
            endif
         endif
      endif
      EVOld(:,jdx)=EVAll(:,jdx)
   enddo
endif

! write output
if(myid==0) then
   if (iii==nloop) then
      do idx=1,Naddlast
         write(200,78) (EVAll(Npoints+idx,nev-jdx),jdx=0,nev-1)
      enddo
   else
      do idx=1,Nadd
         write(200,78) (EVAll(Npoints+idx,nev-jdx),jdx=0,nev-1)
      enddo
   endif
endif

! write the eigenvectors to the file
if(myid==0) then
   write(output,'(a,a,i5,a)') output_file(1:onlen),'_',iii+90000,'.ev'
   open(50,file='embed/'//output,status='unknown')
   do idx=1,NpointsNew
      write(50,78) (EVAll(idx,nev-jdx),jdx=0,nev-1)
   enddo
   close(50)
78 format(10(d15.6,1x))
endif

! write the eigenvalues to the file
if(myid==0) then
   write(output,'(a,a,i5,a)') output_file(1:onlen),'_',iii+90000,'.eg'
   open(50,file='embed/'//output,status='unknown')
   do idx=1,nconv
      write(50,'(10(f9.6,1x))',advance='no') d(idx,1)
   enddo
   write(50,*)
   close(50)
endif

! deallocate all the valuables in the loop
if(myid==0) then
   deallocate(EVAll)
endif
deallocate(sumrow)
deallocate(EV)
deallocate(vv)
deallocate(ww)
deallocate(v_all)

!
! end loop iii for adding points
!
enddo

if(myid==0)then
   close(30)
   close(200)
   deallocate(EVOld)
endif
deallocate(LocalScale)
deallocate(Weight)
deallocate(iloc)
deallocate(dloc)
deallocate(jdxloc)
deallocate(NallPoint)
deallocate(idxlocPoint)



call MPI_finalize(ierr)

end program wlsdmap_embed

!!-------------------------------------------------------------------
!!     matrix vector multiplication subroutines
!!-------------------------------------------------------------------


subroutine av(sumi, NlocCluster, Ncluster, dloc, jdxloc, iloc, v, w)
!!
!! calculate w <--- Mv    M is a block matrix with CSR sparse storage and block size in vector cl
!!
  integer :: sumi, NlocCluster, Ncluster, idx, kdx
  integer :: jdxloc(sumi), iloc(NlocCluster+1)
  real*8 :: v(Ncluster), w(NlocCluster), dloc(sumi)
  real*8 :: sum
  do idx=1,NlocCluster
     sum=0.
     do kdx=iloc(idx),iloc(idx+1)-1
        sum=sum+dloc(kdx)*v(jdxloc(kdx))
     enddo
     w(idx)=sum
  enddo
  return
end subroutine av

subroutine current_time()
  character(8)::date1
  character(10)::time1
  character(5)::zone
  integer::values1(8)
  call date_and_time(date1,time1,zone,values1)
  write(6,'(12(a))',advance='no') date1(1:4),'/',date1(5:6),'/',date1(7:8),' ',&
   time1(1:2),':',time1(3:4),' ',time1(5:10),' : ' 
end subroutine current_time
