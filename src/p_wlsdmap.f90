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
! See "wlsdmap.input_example" for an example of input.
!
! Output
! Eigenvalues and eigenvectors (diffusion coordinates) of Fokker-Planck operator
! 
! Please reference the papers below if you use LSDMap:
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., 
!    J. Chem. Phys., 134, 124116, 2011
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R., and Clementi, C.,
!    J. Phys. Chem. B, 115, 13065-13074, 2011
!
! WZ Jul 25 2012
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       Develop log
! 
!       WZ 2/2/2011
!       Weighted LSDMap
!
!       WZ 11/04/2010
!       Cluster LSDMap
!
!       WZ 2/5/2010
!       Use double precision.
!
!       WZ 5/22/2009
!       Sparse matrix version.
!       The matrix are stored in Compressed Sparse Row(CSR) structure.
!
!       WZ 5/12/2009
!       Parallel version
!       PARPACK package is used to find the eigenvalues of the matrix.
!
!	CC 12/10/2007 
!	Take nearest neighbors data (pre-computed) and 
!	obtain diffusion coordinates 
!	Adapted from matlab code of M.Maggioni @Duke
!	The math foundation is from papers of Coifman et al.
!	(see PNAS 102 (2005), pp. 7426-7431)
!	ARPACK soubroutine DSSIMP is used to find eigenvalues 
!	and eigenvectors of a NxN real symmetric matrix
!	--> see example on how to use DSSIMP at the webpage:
!	http://www.caam.rice.edu/software/ARPACK/UG/node22.html#SECTION00630000000000000000
!
!----------------------------------------------------------------------


subroutine Weighted_LSDMap

use parallel, only :  size, rank, ierr, comm
use data, only :  Npoints,nloc,nstart,idneigh,tmp_rmsd,FullEpsArray,NN_input_weight,output_file,    & 
                  status_dmap,status_eps,column_eps,cutoff,eps0,current_time


implicit none

include 'mpif.h'
! valuables for ARPACK
include 'debug.h'
include 'stat.h'

integer :: NlocPoint,i
integer, allocatable :: NallPoint(:)
integer, allocatable :: idxlocPoint(:)
integer, allocatable :: displacements2(:),counts2(:)
real(kind=8), allocatable :: Weight(:)
real(kind=8), allocatable :: LocalScale(:),lScale(:)
real(kind=8), allocatable :: EV(:,:),EVAll(:,:)
real(kind=8) :: tmpreal1
integer :: inlen, onlen
integer :: ndx,idx,jdx,kdx, iternum, sumi, count,reci
! eps0: eps0 in nlocal=0 mode (input)
! cutoff: sparse matrix cutoff, the entry which is smaller than cutoff will be set as zero
real(kind=8) :: distLoc, distt
! sumrow: the sum of row
! sumloc: the sum of row in one core
! vv(:): one column of the basis in one core
real(kind=8), allocatable :: dloc(:),sumrow(:),sumloc(:),vv(:),v_all(:),ww(:)
! tmp_rmsd: temporarily store the rmsd data from the input file
real,allocatable :: rmsd(:)
! sparse storage
integer, allocatable :: jdxloc(:),iloc(:)
! input and output file name
character(180) :: NN_input_rmsd, NN_input_eps, NN_input_loc
character(200) :: output
! PARACK's tol
real(kind=8), parameter ::  zero = 0.0D+0
! nspace: the ratio of the space needed for the sparse matrix
real(kind=8), parameter :: nspace=1.

real(kind=8), allocatable :: tmp_eps(:)

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
real(kind=8),allocatable :: v(:,:), workl(:),workd(:), d(:,:), resid(:),ax(:)

logical,allocatable :: select(:)
integer,allocatable :: iparam(:), ipntr(:)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
character        bmat*1, which*2
integer          ido, n, nev, ncv, lworkl, info, j,nconv, maxitr, mode, ishfts
logical          rvec
real(kind=8)  tol, sigma
!
!     %----------------------------------------------%
!     | Local Buffers needed for MPI communication |
!     %----------------------------------------------%
!
real(kind=8),allocatable ::  mv_buf(:)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
real(kind=8) pdnorm2 
external pdnorm2, daxpy
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic :: abs


allocate(v(ldv,maxncv), workl(maxncv*(maxncv+8)),workd(3*maxnloc), d(maxncv,2), resid(maxnloc),ax(maxnloc),mv_buf(maxnloc))
allocate(select(maxncv))
allocate(iparam(11), ipntr(11))

!if (rank==0) write(6,*) 'Subroutine p_wlsdmap Start...'
if (rank==0) call current_time('Subroutine p_wlsdmap Start...')

NN_input_rmsd = 'aladip.xyz_rmsd'
NN_input_eps  = 'aladip.xyz_eps'

nev=10

if (rank == 0) then
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

! read weight
allocate(Weight(Npoints))
if(rank==0) then
   if(status_dmap==1)then
      open(10,file=NN_input_weight,status='old')
      do idx=1,Npoints
         read(10,*) Weight(idx)
!         Weight(idx)=exp(Weight(idx)*beta)
      enddo
      close(10)
   endif
endif
if(status_dmap==1) then 
   call MPI_BCAST(Weight,Npoints,MPI_DOUBLE_PRECISION,0,comm,ierr)
else
   Weight=1.
endif

! split data
NlocPoint=int(Npoints/size)
if(rank<mod(Npoints,size)) then
   NlocPoint=NlocPoint+1
endif
allocate(NallPoint(size))
call MPI_ALLGATHER(NlocPoint,1,MPI_INTEGER,NallPoint,1,MPI_INTEGER,comm,ierr)

! build the index for NlocPoint
allocate(idxlocPoint(size+1))
idxlocPoint(1)=1
do idx=2,size+1
   idxlocPoint(idx)=idxlocPoint(idx-1)+NallPoint(idx-1)
enddo

! generate the input file name
inlen = index(NN_input_rmsd,' ') - 1
write(NN_input_loc,'(a,a,i7,a,i7)') NN_input_rmsd(1:inlen),'_',9000000+idxlocPoint(rank+1),'_',9000000+idxlocPoint(rank+2)-1

! read local scale
allocate(LocalScale(Npoints))

if(rank==0) then
   if(status_eps==1)then
      allocate(tmp_eps(column_eps))
      do idx=1,Npoints
         LocalScale(idx)=FullEpsArray(column_eps,idx)
      enddo
      deallocate(tmp_eps)
   endif
endif

deallocate(FullEpsArray)

if(status_eps==1) then 
   call MPI_BCAST(LocalScale,Npoints,MPI_DOUBLE_PRECISION,0,comm,ierr)
else
   LocalScale=eps0
endif
call MPI_barrier(comm,ierr)

!!
!! allocate block CSR storage
!!
! CSR: iloc,jdxloc,dloc
!
if (rank==0) then
   write(6,'(a)') 'allocate block CSR storage and create similarity matrix...'   
endif
allocate(iloc(NlocPoint+1))
sumi=ceiling(Npoints*nspace*NlocPoint)
allocate(dloc(sumi))
allocate(jdxloc(sumi))
iloc=0
count=0
do jdx=1,nloc
   reci = jdx + nstart -1
   iloc(jdx)=count+1
   ndx=idxlocPoint(rank+1)+jdx-1
   do kdx=1,Npoints
      distLoc=tmp_rmsd(kdx,reci)**2  
      distt=Weight(ndx)*Weight(kdx)*exp(-distLoc/Localscale(ndx)/LocalScale(kdx)/2.)
      if (distt >= cutoff) then
         count=count+1
         dloc(count)=distt
         jdxloc(count)=kdx
      endif
   enddo
enddo

iloc(NlocPoint+1)=count+1
if(count>sumi) then
   print *,'Warning: core #',rank,' does not have enough memory to store the sparse matrix'
   print *,'core #',rank,' has used',count,'/',sumi, ' sparse storage'
   print *,'Please increase the value nspace in the program'
   stop
endif
print *,'core #',rank,' has used',count,'/',sumi, '% sparse storage'

! valuables for normalization
allocate(sumrow(Npoints))
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

call MPI_ALLGATHERV(sumloc,NlocPoint,MPI_DOUBLE_PRECISION,sumrow,NallPoint,idxlocPoint(1:size)-1,MPI_DOUBLE_PRECISION,comm,ierr)
if (rank==0) then
   write(6,'(a)') 'all the cores finish sum over rows 1, sum has been gathered'
endif

!!
!! rescale matrix...
!!
if (rank==0) then
   write(6,'(a)') 'rescale matrix...'
endif
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      dloc(kdx)=dloc(kdx)*sqrt(Weight(idxlocPoint(rank+1)+jdx-1)*Weight(jdxloc(kdx))/(sumrow(idxlocPoint(rank+1)+jdx-1)*sumrow(jdxloc(kdx))))
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
call MPI_ALLGATHERV(sumloc,NlocPoint,MPI_DOUBLE_PRECISION,sumrow,NallPoint,idxlocPoint(1:size)-1,MPI_DOUBLE_PRECISION,comm,ierr)
if (rank==0) then
   write(6,'(a)') 'all the cores finish sum over rows 1, sum has been gathered'
endif
!!
!! rescale matrix...
!!
if (rank==0) then
   write(6,'(a)') 'rescale matrix...'
endif
do jdx=1,NlocPoint
   do kdx=iloc(jdx),iloc(jdx+1)-1
      dloc(kdx)=dloc(kdx)/sqrt(sumrow(idxlocPoint(rank+1)+jdx-1)*sumrow(jdxloc(kdx)))
   enddo
enddo
deallocate(sumloc)
!! sumrow is needed in basis calculation
sumrow=1/sqrt(sumrow)

!! 
!! deallocate all the valueables that are not needed in eigenvalue decomposition
!!
deallocate(LocalScale)

call MPI_barrier(comm, ierr)

! valuables for matrix multiplication
allocate(vv(NlocPoint))
allocate(ww(NlocPoint))
allocate(v_all(NPoints))

!--------------------------------------------------------------
!       FROM HERE IT COMPUTES THE EIGENVALUES AND EIGENVECTORS 
!       OF THE DIFFUSION MATRIX 
!--------------------------------------------------------------

if (rank == 0) then
   write(6,'(a)') 'finished creating diffusion matrix... start diagonalization...'
endif

ndigit = -3
logfil = 6
msaupd = 1

! nev: number of eigenvalues
n=Npoints
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
call pdsaupd (comm, ido, bmat, NlocPoint, which, nev, tol, resid,ncv, v, ldv, iparam, ipntr, workd, workl,lworkl, info )

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
   call MPI_ALLGATHERV(vv, NlocPoint, MPI_DOUBLE_PRECISION, v_all, NallPoint, idxlocPoint(1:size)-1, MPI_DOUBLE_PRECISION, comm,ierr)   

   call av (sumi, NlocPoint, Npoints,dloc, jdxloc, iloc, v_all, ww)
   
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
   if ( rank .eq. 0 ) then
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
   if (rank==0) then
      write(6,'(a)') 'no fatal errors occurred in iteration'
      write(6,'(a)') 'start calling psseupd subroutine...'
   endif
   call pdseupd ( comm, rvec, 'All', select,d, v, ldv, sigma,bmat, NlocPoint, which, nev, tol, resid, ncv, v, ldv,iparam, ipntr, workd, workl, lworkl, ierr )
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
      if ( rank .eq. 0 ) then
         print *, ' '
         print *, ' Error with _seupd, info = ', ierr
         print *, ' Check the documentation of _seupd. '
         print *, ' '
      endif
   else
      nconv =  iparam(5)
      if (rank==0) then
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
         call MPI_ALLGATHERV(v, NlocPoint, MPI_DOUBLE_PRECISION, v_all, NallPoint, idxlocPoint(1:size)-1, MPI_DOUBLE_PRECISION, comm,ierr)   

         call av (sumi, NlocPoint, Npoints,dloc, jdxloc, iloc, v_all, ax)

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
   if (rank .eq. 0)then
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
      print *, ' The number of processors is ', size
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

if (rank==0) then
   write(6,'(a)') 'finish diagonalization... start postprocessing....'
endif

! gather the eigenvectors
allocate(EV(NlocPoint,nev))
do jdx=1,nconv
   do idx=1,NlocPoint
      EV(idx,jdx)=sumrow(idxlocPoint(rank+1)-1+idx)*v(idx,jdx)
   enddo
enddo
if (rank==0) then
   allocate(EVAll(Npoints,nev))
endif
do jdx=1,nconv
   call MPI_GATHERV(EV(:,jdx),NlocPoint,MPI_DOUBLE_PRECISION,EVAll(:,jdx),NallPoint,idxlocPoint(1:size)-1,MPI_DOUBLE_PRECISION,0,comm,ierr)
enddo

! normalization
if(rank==0) then
   do jdx=1,nconv
      tmpreal1=0.
      do kdx=1,Npoints
         tmpreal1=tmpreal1+EVAll(kdx,jdx)**2
      enddo
      tmpreal1=sqrt(tmpreal1)
      if (EVAll(1,jdx)<0) then
         tmpreal1=-tmpreal1
      endif
      do kdx=1,Npoints
         EVAll(kdx,jdx)=EVAll(kdx,jdx)/tmpreal1
      enddo
   enddo
endif

! write the eigenvectors to the file
if(rank==0) then
   write(output,'(a,a)') output_file(1:onlen),'.ev'
   open(50,file=output,status='unknown')
   do idx=1,Npoints
      write(50,78) (EVAll(idx,nev-jdx),jdx=0,nev-1)
   enddo
   close(50)
78 format(10(d15.6,1x))
endif

! write the eigenvalues to the file
if(rank==0) then
   write(output,'(a,a)') output_file(1:onlen),'.eg'
   open(50,file=output,status='unknown')
   do idx=1,nconv
      write(50,'(10(f9.6,1x))',advance='no') d(idx,1)
   enddo
   write(50,*)
   close(50)
endif


! deallocate all the valuables in the loop
if(rank==0) then
   deallocate(EVAll)
endif
deallocate(sumrow)
deallocate(EV)
deallocate(Weight)
deallocate(NallPoint)
deallocate(idxlocPoint)
deallocate(iloc)
deallocate(dloc)
deallocate(jdxloc)
deallocate(vv)
deallocate(ww)
deallocate(v_all)

if(rank==0) call current_time('Subroutine p_wlsdmap ends...')

end subroutine Weighted_LSDMap

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

