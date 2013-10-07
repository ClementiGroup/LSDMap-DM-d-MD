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
! LSDMap v1.1 - Sept 2013 - Merger Release                      !                                                         
!                                                               !                                                          
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uonversity of Edinburgh                 !
!---------------------------------------------------------------!

! Parallel Local PCA analysis
!
! Input
! 1. Name of the trajectory
! 2. Name of the nearest neighboring map from p_rmsd_neighbor
! 3. Start and end point ID
! 4. Job ID
! 5. Number of neighbors, number of point for MDS,
!    start point for the MDS spectra, the step size for the MDS spectra
! 6. cutoff for the first derivative of the MDS spectra
! 7. number of CPUs to share the trajectory
!
! See "local_mds.input_sample" for an example of input.
!
! Output
! Local scale as the input of p_wlsdmap
! 
! Please reference the papers below if you use LSDMap:
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., 
!    J. Chem. Phys., 134, 124116, 2011
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R., and Clementi, C.,
!    J. Phys. Chem. B, 115, 13065-13074, 2011
!
! WZ Jul 25 2012
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! Develop log
!
! WZ Aug 30 2011
! Take the original trajectory and nearest neighbor graph as the input.
! Share the trajectory in every <ncore> CPUs as a small group.
!
! WZ Jul 28 2011
! 1. Implement random projection multidimensional scaling.
! The old classical MDS function is still at the end of the program.
! 2. Implement the idea to use k1 of the first k neighbors for the 
! multidimensional scaling as an approximation of the spectra of 
! all the first k neighbors.  The accuracy depends on the ratio k1/k.
! If k1 equals to k, the results are the same as the old version.
!
! WZ Jan 22 2010
! Add polynomial fit to the spectra.
! Picking local scale with the help of spectra gap and first derivative
! of spectra.
! 
!-----------------------------------------------------------------------

subroutine p_local_mds

use Qsort_Module
use iso_c_binding
use ftn_c
use data, only :  ns,ne,Npoints,dim,Nneigh,Natoms,nloc,nend,nlimit,traj,idneigh,EpsArray,   &
                 FullEpsArray,ncore,dmds,kmin,dk,neps,seps,deps, current_time

use parallel, only :  size, rank, ierr, comm, counts

implicit none

include 'mpif.h'

integer :: ntotal,color1,key1,subcomm1,color2,key2,subcomm2
integer :: kmax,inlen,it_rpmds,d_rpmds
integer :: mdx,idx,jdx,kdx,stepsize
real(kind=8) :: sumall
integer :: tmp_int1
!integer,allocatable :: nlocall(:),idxlocPoint(:)
real(kind=8),allocatable :: distLocal(:,:)
real(kind=8),allocatable :: anapca(:,:)
integer,allocatable :: cneighbor(:),idneighbor(:),tmpneighbor(:),neighbor(:,:)
real(kind=8),allocatable :: eigval(:)
real(kind=8),allocatable :: coor(:,:)
real(kind=8),allocatable :: tdist(:,:)
character(30) :: tmp_text
integer,parameter :: dspec=30,dnev=30
integer :: smax
real(kind=8),dimension(2) :: res_real
integer,dimension(3) :: res_int
real(kind=8),allocatable :: vx(:),xx(:,:),yy(:,:),rot(:),weight(:)
real(kind=8):: ceps
integer :: nloctraj,nstraj,netraj,nloc_neighbor
integer,allocatable :: nstrajall(:),point_neighbor(:),nlocall_neighbor(:)
real,allocatable :: loctraj(:,:)
real,allocatable :: locpicktraj(:,:),picktraj(:,:),ctraj(:,:)
integer,parameter :: sizeofinteger=4


integer :: mynstart, mynend,mynloc,i, myeps
integer, allocatable :: dspls(:)

if (rank==0) call current_time('subroutine p_local_mds start')

kmax=dmds
smax=(kmax-kmin)/dk+1

! split the processors into small groups
! the first type of spliting
! <0 1 2> <3 4 5> <6 7 8>...
color1=int(rank/ncore/1.)
key1=mod(rank,ncore)
call MPI_COMM_SPLIT(comm,color1,key1,subcomm1,ierr)
! the second type of spliting
! <0 3 6> <1 4 7> <2 5 8>...
color2=mod(rank,ncore)
key2=int(rank/ncore/1.)
call MPI_COMM_SPLIT(comm,color2,key2,subcomm2,ierr)

! split data
ntotal=ne-ns+1
if(rank>=mod(ntotal,size)) then
   mynloc=floor(1.*ntotal/size)
   mynend=ne-(size-rank-1)*mynloc
   mynstart=mynend-mynloc+1
else
   mynloc=floor(1.*ntotal/size)+1
   mynstart=ns+rank*mynloc
   mynend=mynstart+mynloc-1
endif

myeps = 5*neps + 1

if(rank==0)then
   allocate(dspls(0:size-1))  ! deallocate later
   dspls(0) = 0
   do i = 1, size-1
      dspls(i) = dspls(i-1) + counts(i-1)*myeps
   enddo
endif


allocate(neighbor(Nneigh,nloc))
do idx = 1,nloc
   do jdx = 1,Nneigh
      neighbor(jdx,idx) = idneigh(jdx,idx)
   enddo
enddo

! write status
if(rank==size-1) then
   print *,'From subroutine p_local_mds:'
   print *,'number of points in dataset = ',Npoints
   print *,'original dimension = ',dim
   print *,'number of nearest neighbors used = ', Nneigh
   print *,'from ',ns,' to ',ne
   print *,'dmds=',dmds,'kmin=',kmin,' kmax=',kmax,' dk=',dk
   print *,'number of cores = ',size
   print *,'points per core = ',nloc,'~',nloc+1
   print *,'number of cores per group = ',ncore
   write(*,'(a,f5.2,a,f5.2,a,f5.2)')'cutoff = ',seps,':',deps,':',seps+(neps-1)*deps
endif

! read trajectory in every <ncore> CPUs
if(rank==0) call current_time('start reading trajectory')
nloctraj=floor(1.*Npoints/ncore)
if(key1>=mod(Npoints,ncore)) then
   netraj=Npoints-(ncore-key1-1)*nloctraj
   nstraj=netraj-nloctraj+1
else
   nloctraj=nloctraj+1
   nstraj=1+key1*nloctraj
   netraj=nstraj+nloctraj-1
endif
allocate(nstrajall(size))
call MPI_ALLGATHER(nstraj,1,MPI_INTEGER,nstrajall,1,MPI_INTEGER,comm,ierr)

allocate(loctraj(dim,nloctraj))

!if(rank<ncore)then  !EYB:: Checl what happens if ncore > 1
!   do idx=1,nstraj-1
!   enddo
   do idx=nstraj,netraj
      do jdx = 1, dim
         loctraj(jdx,idx-nstraj+1) = traj(mod((jdx-1),3)+1,((idx-1)*dim/3)+((jdx-1)/3)+1)
      enddo
   enddo
!endif

!call MPI_barrier(comm,ierr)
!call MPI_BCAST(loctraj,dim*nloctraj,MPI_REAL,0,subcomm2,ierr)

! allocate storage
allocate(distLocal(dmds,dmds))
allocate(cneighbor(dmds))
allocate(idneighbor(dmds))
allocate(point_neighbor(ncore))
allocate(nlocall_neighbor(ncore))
allocate(picktraj(dim,dmds))
allocate(ctraj(dim,dmds))
allocate(xx(Natoms,3))
allocate(yy(Natoms,3))
allocate(rot(9))
allocate(weight(Natoms))

weight=1.

it_rpmds=1
d_rpmds=2
stepsize=floor(Nneigh/dmds/1.)

allocate(EpsArray(myeps,mynstart:mynend))

do idx=mynstart,mynstart+nlimit-1
   if (rank==0) then
      write(tmp_text,'(a,i7)') 'start processing point ',idx
      call current_time(tmp_text)
   endif
   if(idx<=nend) then
      ! pick <dmds> neighbors from <Nneigh> neighbors
      do jdx=1,dmds
         cneighbor(jdx)=neighbor(stepsize*(jdx-1)+1,idx-mynstart+1)
      enddo
      ! get trajectory
      ! sort the neighbor list
      idneighbor=(/(jdx,jdx=1,dmds,1)/)
      call iQsort(cneighbor,idneighbor)

   else
      ! if the current CPU ends, to avoid errors in the next data sharing block,
      ! set all the neighbors the first point
      cneighbor=1
   endif
   do jdx=0,ncore-1
      allocate(tmpneighbor(dmds))
      tmpneighbor=cneighbor
      call MPI_BCAST(tmpneighbor,dmds,MPI_INTEGER,jdx,subcomm1,ierr)
      ! get the pointer to the neighbor list
      tmp_int1=ncore ! the current rank of core in the subgroup to find the beginning
      do kdx=dmds,1,-1
         do while(nstrajall(tmp_int1)>tmpneighbor(kdx))
            point_neighbor(tmp_int1)=kdx+1
            tmp_int1=tmp_int1-1
         enddo
      enddo
      point_neighbor(1)=1
      do kdx=1,ncore-1
         nlocall_neighbor(kdx)=point_neighbor(kdx+1)-point_neighbor(kdx)
      enddo
       nlocall_neighbor(ncore)=dmds-point_neighbor(ncore)+1
      nloc_neighbor=nlocall_neighbor(key1+1)
      allocate(locpicktraj(dim,nloc_neighbor))
      do kdx=1,nloc_neighbor
         locpicktraj(:,kdx)=loctraj(:,tmpneighbor(point_neighbor(key1+1)+kdx-1)-nstraj+1)
      enddo
      call MPI_GATHERV(locpicktraj,dim*nloc_neighbor,MPI_REAL,picktraj,dim*nlocall_neighbor, &
           dim*(point_neighbor-1),MPI_REAL,jdx,subcomm1,ierr)
      deallocate(locpicktraj)
      deallocate(tmpneighbor)
   enddo
   if(idx<=nend)then
      ! get trajectory in the original neighbor list order
      do jdx=1,dmds
         ctraj(:,idneighbor(jdx))=picktraj(:,jdx)
      enddo
      ! calculate rmsd
      do jdx=1,dmds
         do kdx=jdx+1,dmds
            xx(:,1)=ctraj(1:dim-2:3,jdx)
            xx(:,2)=ctraj(2:dim-1:3,jdx)
            xx(:,3)=ctraj(3:dim:3,jdx)
            yy(:,1)=ctraj(1:dim-2:3,kdx)
            yy(:,2)=ctraj(2:dim-1:3,kdx)
            yy(:,3)=ctraj(3:dim:3,kdx)
            sumall=calcrmsdrotationalmatrix(Natoms,xx(:,1),xx(:,2),xx(:,3),yy(:,1),&
                 yy(:,2),yy(:,3),rot,weight)
            sumall=sumall**2
            distLocal(jdx,kdx)=sumall
            distLocal(kdx,jdx)=sumall
         enddo
      enddo
      do jdx=1,dmds
         distLocal(jdx,jdx)=0.
      enddo
      ! MDS
      allocate(anapca(smax,dspec))
      allocate(vx(smax))
      anapca=0.
      do jdx=1,smax
         vx(jdx)=sqrt(distLocal(1,(jdx-1)*dk+kmin))
      enddo
      do jdx=kmin,kmax,dk
         allocate(eigval(dnev))
         allocate(coor(jdx,dnev))
         allocate(tdist(jdx,jdx))
         tdist=distLocal(1:jdx,1:jdx)
         call rpmds(tdist,jdx,dnev,eigval,it_rpmds,d_rpmds)
         deallocate(tdist)
         mdx=(jdx-kmin)/dk+1
         if (jdx<dnev .and. jdx<dspec) then
            do kdx=1,jdx
               anapca(mdx,kdx)=eigval(kdx)/sqrt(real(jdx))
            enddo
         else
            do kdx=1,dspec
               anapca(mdx,kdx)=eigval(kdx)/sqrt(real(jdx))
            enddo
         endif
         deallocate(eigval)
         deallocate(coor)
      enddo
      EpsArray(1,idx) = idx
      do jdx=1,neps
         ceps=seps+(jdx-1)*deps
         call pickeps(anapca,vx,smax,dspec,res_real,res_int,ceps)
         EpsArray((jdx-1)*5+2,idx) = ((res_int(1)-1)*dk+kmin-1)*stepsize+1
         EpsArray((jdx-1)*5+3,idx) = res_real(1)
         EpsArray((jdx-1)*5+4,idx) = res_real(2)
         EpsArray((jdx-1)*5+5,idx) = res_int(2)
         EpsArray((jdx-1)*5+6,idx) = res_int(3)
      enddo
      deallocate(vx)
      deallocate(anapca)
   endif ! matching line 319 if(id<=nend)
enddo

allocate(FullEpsArray(myeps,Npoints))
call MPI_GATHERV(EpsArray,myeps*nloc,MPI_REAL8,FullEpsArray,myeps*counts,dspls,MPI_REAL8,0,comm,ierr)
deallocate(EpsArray)

if (rank == size-1 .and. nend ==Npoints) then
   write(*,'(a,f10.6,a,f10.6,a,f10.6)') '% the first derivative cutoff = ', seps,':',deps,':',&
        seps+(neps-1)*deps
   write(*,'(a)') '% id/k_1/eps_1/eiggap_1/intdim1_1/intdim2_1/k_2/eps_2/eiggap_2/intdim1_2/intdim2_2/...'
endif

!close(51)


if (rank==0) call current_time('Subroutine p_local_mds end')
if(rank==0)then 
   deallocate(dspls)
endif

end subroutine p_local_mds

! random projection multidimensional scaling
! if it<=1, it is faster than the classical mds function using ARPACK
! 
subroutine rpmds(distLocal,dmds,dnev,eigval,it,d)
  implicit none
  integer :: idx,jdx,kdx,dmds,dnev,dl,it,d,lworkqr,lworksvd,info
  real(kind=8),dimension(dmds,dmds) :: distLocal
  real(kind=8),dimension(dnev) :: eigval
  real(kind=8),dimension(dmds,1) :: sumid !,ones
  real(kind=8) :: sumall,r1
  real(kind=8),dimension(dmds,dnev+d) :: mat_rand,mat_h,mat_h1
  real(kind=8),dimension(dmds,(dnev+d)*(it+1)) :: mat_f,mat_f1
  real(kind=8),dimension((dnev+d)*(it+1)) :: tau,workqr
  real(kind=8),allocatable :: worksvd(:),U(:,:),VT(:,:),S(:)

  ! center the data
  sumid(1:dmds,1)=sum(distLocal,dim=1)/dmds
  sumall=sum(sumid)/dmds
  do jdx=1,dmds
     do kdx=1,dmds
        distLocal(jdx,kdx)=(distLocal(jdx,kdx)-sumid(jdx,1)-sumid(kdx,1)+sumall)*(-0.5)
     enddo
  enddo

!!
!! Here is the matlab code for centering the point.  
!! The code above is two times faster than 
!! directly translating the matlab code to fortran.
!!
!  ones=1.
!  dist=dist-matmul(matmul(dist,ones),transpose(ones))/dmds
!  dist=-0.5*(dist1-matmul(ones,matmul(transpose(ones),dist1))/dmds)

  ! SVD directly if (it+1)*dl>=dmds/1.25
  dl=dnev+d
  if((it+1)*dl>=dmds/1.25)then
     lworksvd=5*dmds
     allocate(worksvd(lworksvd))
     allocate(S(dmds))
     allocate(U(dmds,dmds))
     allocate(VT(dmds,dmds))
     call dgesvd("N","N",dmds,dmds,distLocal,dmds,S,U,dmds,VT,dmds,worksvd,lworksvd,info)  
     eigval=sqrt(S(1:dnev))
     return
  endif
  ! generate random matrix
  do idx=1,dmds
     do jdx=1,dl
        call random_number(r1)
        mat_rand(idx,jdx)=r1*2-1
     enddo
  enddo
  ! calculate matrix H=A*mat_rand
  !mat_h=matmul(dist,mat_rand)
  call dgemm('N','N',dmds,dnev+d,dmds,1.d0,distLocal,dmds,mat_rand,dmds,0.d0,mat_h,dmds)
  ! calculate matrix F=[H A'AH (A'A)^2*H ...]
  mat_f(:,1:dl)=mat_h
  do idx=1,it
!     mat_h1=transpose(matmul(transpose(mat_h),dist))
     call dgemm('T','N',dmds,dnev+d,dmds,1.d0,distLocal,dmds,mat_h,dmds,0.d0,mat_h1,dmds)
!     mat_h=matmul(dist,mat_h1)
     call dgemm('N','N',dmds,dnev+d,dmds,1.d0,distLocal,dmds,mat_h1,dmds,0.d0,mat_h,dmds)
     mat_f(:,dl*idx+1:dl*(idx+1))=mat_h
  enddo
  ! QR factorization of F
  lworkqr=dl*3
  call dgeqrf(dmds,dl*(it+1),mat_f,dmds,tau,workqr,lworkqr,info)
  call dorgqr(dmds,dl*(it+1),dl*(it+1),mat_f,dmds,tau,workqr,lworkqr,info)
  ! SVD for A'*Q
!  mat_f1=transpose(matmul(transpose(mat_f),dist))
  call dgemm('N','N',dmds,(dnev+d)*(it+1),dmds,1.d0,distLocal,dmds,mat_f,dmds,0.d0,mat_f1,dmds)
  lworksvd=max(dmds+3*dl*(it+1),5*dl*(it+1))
  allocate(worksvd(lworksvd))
  allocate(S(dl*(it+1)))
  allocate(U(dmds,dl*(it+1)))
  allocate(VT(dl*(it+1),dl*(it+1)))
  call dgesvd("N","N",dmds,dl*(it+1),mat_f1,dmds,S,U,dmds,VT,dl*(it+1),worksvd,lworksvd,info)  
  eigval=sqrt(S(1:dnev))
  return
end subroutine rpmds

! subroutine smooth(kmax,x,y,n)
!   integer :: kmax,n,jdx,ndx
!   real,dimension(kmax) :: x, y
! ! analyze the local pca data and get epsilon
!   real,allocatable :: epsn(:),pcan(:)
!   allocate(epsn(kmax+1))
!   allocate(pcan(kmax+2))
! ! interative methods to smooth the curve
!   epsn(1)=x(2)-x(1)
!   do jdx=1,kmax-1
!      epsn(jdx+1)=x(jdx+1)-x(jdx)
!   enddo
!   epsn(kmax+1)=x(kmax)-x(kmax-1)
!   do jdx=1,kmax+1
!      if (epsn(jdx) == 0) then
!         ! set 1/epsn large enough when epsn=0
!         epsn(jdx)=10000000.
!      else
!         epsn(jdx)=1/abs(epsn(jdx))
!      endif
!   enddo
!   do ndx=1,n
!      pcan(1)=y(1)
!      do jdx=1,kmax
!         pcan(jdx+1)=y(jdx)
!      enddo
!      pcan(kmax+2)=y(kmax)
!      do jdx=1,kmax
!         y(jdx)=(pcan(jdx)*epsn(jdx)+pcan(jdx+2)*epsn(jdx+1)+pcan(jdx+1))/(epsn(jdx)+epsn(jdx+1)+1.)
!      enddo
!   enddo
!   deallocate(epsn)
!   deallocate(pcan)
!   return
! end subroutine smooth

subroutine pickeps(anapca,vx,smax,dspec,res_real,res_int,cuteps)
  implicit none
! begin analyze the singular value spectra
  integer :: smax,dspec,found
  real(kind=8),dimension(smax,dspec) :: anapca
  real(kind=8),dimension(smax,dspec) :: ndata
  real(kind=8),dimension(smax,dspec) :: dUdX
  real(kind=8),dimension(smax,dspec-1):: dU
  real(kind=8),dimension(dspec,6) :: fit
  real(kind=8),dimension(smax) :: vx,vy
! 1 k 2 noise 3 new_noise
  real(kind=8),dimension(2)::res_real
! 1 epsilon 2 gap
  integer,dimension(3)::res_int
  integer,dimension(3) :: mp
  integer :: mc,ct,noi,noise,ndx,jdx,kdx
  real(kind=8),parameter :: spe=0.5
  real(kind=8) :: cuteps
  integer, dimension(3,dspec-6)::nstat
  real(kind=8) :: z1,z2,z3,z4,z5
  integer :: d
! do 5th order polynomial fit to all the spectra
  do ndx=1,dspec
     do jdx=1,smax
        vy(jdx)=anapca(jdx,ndx)
     enddo
     d=5
     call polyfit(smax,vx,vy,fit(ndx,:),d)
  enddo
! generate new spectra with ploynomial fit results
  do ndx=1,dspec
     do jdx=1,smax
        ndata(jdx,ndx)=fit(ndx,6)*vx(jdx)**5+fit(ndx,5)*vx(jdx)**4+fit(ndx,4)*vx(jdx)**3+ &
             fit(ndx,3)*vx(jdx)**2+fit(ndx,2)*vx(jdx)+fit(ndx,1)
     enddo
  enddo
! calculate the first derivative of the spectra
  do ndx=1,dspec
     do jdx=1,smax
        dUdX(jdx,ndx)=5*fit(ndx,6)*vx(jdx)**4+4*fit(ndx,5)*vx(jdx)**3+3*fit(ndx,4)*vx(jdx)**2+ &
             2*fit(ndx,3)*vx(jdx)+fit(ndx,2)
     enddo
  enddo
! calculate the spectra gap
  do ndx=1,dspec-1
     do jdx=1,smax
        dU(jdx,ndx)=abs(ndata(jdx,ndx)-ndata(jdx,ndx+1))
     enddo
  enddo
! pick three scales: 3/7, 1/2 and 4/7
  do ndx=1,smax
     if(vx(ndx) > (3./7.*(vx(smax)-vx(1))+vx(1)) ) then
        mp(1)=ndx
        exit
     endif
  enddo
  do ndx=mp(1),smax
     if(vx(ndx) > (1./2.*(vx(smax)-vx(1))+vx(1)) ) then
        mp(2)=ndx
        exit
     endif
  enddo
  do ndx=mp(2),smax
     if(vx(ndx) > (4./7.*(vx(smax)-vx(1))+vx(1)) ) then
        mp(3)=ndx
        exit
     endif
  enddo
! get status matrix
  nstat=0
  do ndx=1,3
     mc=mp(ndx);
     do jdx=1,dspec-6
        z1=spe*dU(mc,jdx)-dU(mc,jdx+1);
        z2=spe*dU(mc,jdx)-dU(mc,jdx+2);
        z3=spe*dU(mc,jdx)-dU(mc,jdx+3);
        z4=spe*dU(mc,jdx)-dU(mc,jdx+4);
        z5=spe*dU(mc,jdx)-dU(mc,jdx+5);
        if( z1>0 .and. z2>0 .and. z3>0 .and. z4>0 .and. z5>0 ) then
           nstat(ndx,jdx)=1;
        endif
     enddo
  enddo
! find where is the beginning of noise spectra
  ct=-1;
  noise=0;
  do ndx=1,dspec-6
     if ( nstat(1,ndx) /= 0 .or. nstat(2,ndx) /= 0 .or. nstat(3,ndx) /= 0) then
        ct=0
     else if (ct /= -1) then
        ct=ct+1
        if(ct == 3) then
           noise=ndx-2
           res_int(2)=noise
           exit
        endif
     endif
  enddo
  if (noise==0) then
!     print *,'Warning: cannot find noise with step 2!!!'
     noise=dspec-6
     res_int(2)=noise
  endif
! check the first derivative of the noise spectra to get the local scale
  found=0
  do ndx=noise,dspec
     noi=ndx
     do jdx=1,smax
        mc=1
        do kdx=noi,dspec
           if(dUdX(jdx,kdx)>cuteps) then
              mc=0
           endif          
        enddo
        if(mc==1) then
           exit
        endif
     enddo
     if(jdx<smax) then
        res_real(1)=vx(jdx)
        res_int(1)=jdx
        res_int(3)=ndx
        res_real(2)=dU(jdx,1)
        found=1
        exit
     endif
  enddo
  if(found==0)then
!     print *,'Warning: cannot find the local scale!!!'
     res_real(1)=vx(smax)
     res_int(1)=smax
     res_int(3)=dspec
     res_real(2)=dU(smax,1)
  endif
  return
end subroutine pickeps

! code from http://rosettacode.org/wiki/Polynomial_Fitting
subroutine polyfit(smax,vx, vy, fit, d)
  implicit none
  integer :: d,smax
  real(kind=8), dimension(smax) :: vx, vy
  real(kind=8), dimension(d+1) :: fit, tmp_fit
  real(kind=8), dimension(:,:), allocatable :: X, XTX,XTXsvd
!  real(kind=8), dimension(:,:), allocatable :: XT
  real(kind=8), dimension(:,:), allocatable :: XC
  
  integer :: i, j, k
  integer     :: n, lda, lwork
  integer :: info
  integer, dimension(:), allocatable :: ipiv
  real(kind=8), dimension(:), allocatable :: work
  real(kind=8),dimension(:,:),allocatable :: UU,VV
  real(kind=8),dimension(:),allocatable :: S
  fit=0.
  do k=1,d
     n = d-k+2
!     allocate(XT(n, smax))
     allocate(XC(smax,n))
     allocate(X(smax, n))
     allocate(XTX(n, n))
     allocate(XTXsvd(n,n))
     ! prepare the matrix
     do i = 0, n-1
        do j = 1, smax
           X(j, i+1) = vx(j)**i
        end do
     end do
!     XT  = transpose(X)
!     XTX = matmul(XT, X)
     XC=X
     call dgemm('T','N',n,n,smax,1.d0,XC,smax,X,smax,0.d0,XTX,n)
     XTXsvd=XTX
     lwork=5*n
     allocate(work(lwork))
     allocate(S(n))
     allocate(UU(n,n))
     allocate(VV(n,n))
     call DGESVD('N','N',n,n,XTXsvd,n,S,UU,n,VV,n,work,lwork,info) 
     deallocate(UU)
     deallocate(VV)
     deallocate(work)
     deallocate(XTXsvd)
!    check the magnitude of the condition number
     if( log10(S(1))-log10(S(n)) <=12. ) then
        deallocate(S)
        goto 2000
     endif
!     deallocate(XT)
     deallocate(XC)
     deallocate(X)
     deallocate(XTX)
     deallocate(S)
  enddo
2000 continue
  lda = n
  lwork = n
  d= n-1
  allocate(ipiv(n))
  allocate(work(lwork))
  ! calls to LAPACK subs DGETRF and DGETRI
  call DGETRF(n, n, XTX, lda, ipiv, info)
  if ( info /= 0 ) then
     print *, "problem"
     return
  end if
  call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
  if ( info /= 0 ) then
     print *, "problem"
     return
  end if
  
!  fit(1:n) = matmul( matmul(XTX, XT), vy)
!  fit(1:n) = matmul( matmul(XTX, transpose(X)), vy)
  call dgemm('T','N',n,1,smax,1.d0,X,smax,vy,smax,0.d0,tmp_fit(1:n),n)
  call dgemm('N','N',n,1,n,1.d0,XTX,n,tmp_fit(1:n),n,0.d0,fit(1:n),n)
  
  deallocate(ipiv)
  deallocate(work)
  deallocate(X)
!  deallocate(XT)
  deallocate(XC)
  deallocate(XTX)
  
  return
end subroutine polyfit


!subroutine classicalmds(dist,dmds,dnev,eigval,coor,status)
!  integer :: dmds,dnev,status
!  real*8,dimension(dmds,dmds) :: dist
!  real*8,dimension(dnev) :: eigval
!  real*8,dimension(dmds,dnev) :: coor
!  integer :: idx,jdx,kdx
!  real*8,dimension(dmds,1) :: sumid
!  real*8 :: sumall
!  !-----------------------------------------------------------------------
!  ! values for ARPACK
!
!  integer, parameter :: maxn=5000, maxnev=60, maxncv=90,ldv=maxn 
!
!  !     %--------------%
!  !     | Local Arrays |
!  !     %--------------%
!
!  real*8 :: v(ldv,maxncv), workl(maxncv*(maxncv+8)), &
!       &                 workd(3*maxn), d(maxncv,2), resid(maxn), &
!       &                 ax(maxn), vv(ldv,maxncv)
!  logical :: select(maxncv)
!  integer :: iparam(11), ipntr(11)
!  
!  !     %---------------%
!  !     | Local Scalars |
!  !     %---------------%
!  
!  character :: bmat*1, which*2
!  integer :: ido, n, nev, ncv, lworkl, info, ierr, &
!       &                 j, nx, ishfts, maxitr, mode1, nconv
!  logical :: rvec
!  real*8 :: tol, sigma
!  
!  !     %------------%
!  !     | Parameters |
!  !     %------------%
!  
!  real*8, parameter ::  zero = 0.0D+0
!  
!  !     %-----------------------------%
!  !     | BLAS & LAPACK routines used |
!  !     %-----------------------------%
!
!  real*8 :: dnrm2
!  external :: dnrm2, daxpy
!  
!  !     %--------------------%
!  !     | Intrinsic function |
!  !     %--------------------%
!  
!  intrinsic :: abs
!!!
!!! center the data
!!!
!  sumid(1:dmds,1)=sum(dist,dim=1)/dmds
!  sumall=sum(sumid)/dmds
!  do jdx=1,dmds
!     do kdx=1,dmds
!        dist(jdx,kdx)=(dist(jdx,kdx)-sumid(jdx,1)-sumid(kdx,1)+sumall)*(-0.5)
!     enddo
!  enddo
!  
!  !--------------------------------------------------------------
!  !       FROM HERE IT COMPUTES THE EIGENVALUES AND EIGENVECTORS 
!  !       OF THE MDS MATRIX 
!  !--------------------------------------------------------------
!  
!  !     %-------------------------------------------------%
!  !     | The following sets dimensions for this problem. |
!  !   %-------------------------------------------------%
!  
!  n=dmds
!  
!  !     %-----------------------------------------------%
!  !     |                                               | 
!  !     | Specifications for ARPACK usage are set       | 
!  !     | below:                                        |
!  !     |                                               |
!  !     |    1) NEV = 10  asks for 10 eigenvalues to be |  
!  !     |       computed.                               | 
!  !     |                                               |
!  !     |    2) NCV = 30 sets the length of the Arnoldi |
!  !     |       factorization                           |
!  !     |                                               |
!  !     |    3) This is a standard problem              |
!  !     |         (indicated by bmat  = 'I')            |
!  !     |                                               |
!  !     |    4) Ask for the NEV eigenvalues of          |
!  !     |       largest magnitude                       |
!  !     |         (indicated by which = 'LM')           |
!  !     |       See documentation in DSAUPD for the     |
!  !     |       other options SM, LA, SA, LI, SI.       | 
!  !     |                                               |
!  !     | Note: NEV and NCV must satisfy the following  |
!  !     | conditions:                                   |
!  !     |              NEV <= MAXNEV                    |
!  !     |          NEV + 1 <= NCV <= MAXNCV             |
!  !     %-----------------------------------------------%
!   
!  nev   = dnev
!  ncv   = min(nev+60,dmds)
!
!  bmat  = 'I'
!  which = 'LA'   !--- larger algebraic
  
!  if ( n .gt. maxn ) then
!     print *, ' ERROR with _SSIMP: N is greater than MAXN '
!     go to 9000
!  else if ( nev .gt. maxnev ) then
!     print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
!     go to 9000
!  else if ( ncv .gt. maxncv ) then
!     print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
!     go to 9000
!  end if
!   
!   !     %-----------------------------------------------------%
!   !     |                                                     |
!   !     | Specification of stopping rules and initial         |
!   !     | conditions before calling DSAUPD                    |
!   !     |                                                     |
!   !     | TOL  determines the stopping criterion.             |
!   !     |                                                     |
!   !     |      Expect                                         |
!   !     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!   !     |               computed   true                       |
!   !     |                                                     |
!   !     |      If TOL .le. 0,  then TOL <- macheps            |
!   !     |           (machine precision) is used.              |
!   !     |                                                     |
!   !     | IDO  is the REVERSE COMMUNICATION parameter         |
!   !     |      used to specify actions to be taken on return  |
!   !     |      from DSAUPD. (See usage below.)                |
!   !     |                                                     |
!   !     |      It MUST initially be set to 0 before the first |
!   !     |      call to DSAUPD.                                | 
!   !     |                                                     |
!   !     | INFO on entry specifies starting vector information |
!   !     |      and on return indicates error codes            |
!   !     |                                                     |
!   !     |      Initially, setting INFO=0 indicates that a     | 
!   !     |      random starting vector is requested to         |
!   !     |      start the ARNOLDI iteration.  Setting INFO to  |
!   !     |      a nonzero value on the initial call is used    |
!   !     |      if you want to specify your own starting       |
!   !     |      vector (This vector must be placed in RESID.)  | 
!   !     |                                                     |
!   !     | The work array WORKL is used in DSAUPD as           | 
!   !     | workspace.  Its dimension LWORKL is set as          |
!   !     | illustrated below.                                  |
!   !     |                                                     |
!   !     %-----------------------------------------------------%
! 
!   lworkl = ncv*(ncv+8)
!   tol = zero
!   info = 0
!   ido = 0
!   iternum = 0
!   
!   !     %---------------------------------------------------%
!   !     | Specification of Algorithm Mode:                  |
!   !     |                                                   |
!   !     | This program uses the exact shift strategy        |
!   !     | (indicated by setting PARAM(1) = 1).              |
!   !     | IPARAM(3) specifies the maximum number of Arnoldi |
!   !     | iterations allowed.  Mode 1 of DSAUPD is used     |
!   !     | (IPARAM(7) = 1). All these options can be changed |
!   !     | by the user. For details see the documentation in |
!   !     | DSAUPD.                                           |
!   !     %---------------------------------------------------%
!   
!   ishfts = 1
!   maxitr = 5000 
!   mode1 = 1
!   
!   iparam(1) = ishfts
!   iparam(3) = maxitr
!   iparam(7) = mode1
!   
!   !      %------------------------------------------------%
!   !      | M A I N   L O O P (Reverse communication loop) |
!   !      %------------------------------------------------%
!   
! 10 continue
!   
!   !        %---------------------------------------------%
!   !        | Repeatedly call the routine DSAUPD and take | 
!   !        | actions indicated by parameter IDO until    |
!   !        | either convergence is indicated or maxitr   |
!   !        | has been exceeded.                          |
!   !        %---------------------------------------------%
! 
!   !   write(6,'(a)') 'call ssaupd subroutine.....'
! 
!   call dsaupd ( ido, bmat, n, which, nev, tol, resid,            &
!        &                 ncv, v, ldv, iparam, ipntr, workd, workl,    &
!        &                 lworkl, info )
!   iternum = iternum + 1
! 
! !   write(6,'(a,i12)') 'ssaupd subroutine has been used.....', iternum
!    
!   if (ido .eq. -1 .or. ido .eq. 1) then
!      
!      !           %--------------------------------------%
!      !           | Perform matrix vector multiplication |
!      !           |              y <--- OPI*x             |
!      !           | The user should supply his/her own   |
!      !           | matrix vector multiplication routine |
!      !           | here that takes workd(ipntr(1)) as   |
!      !           | the input, and return the result to  |
!      !           | workd(ipntr(2)).                     |
!      !           %--------------------------------------%
!      
!      workd(ipntr(2):ipntr(2)+dmds-1)=matmul(dist,workd(ipntr(1):ipntr(1)+dmds-1))
! !     call av(dmds, dist, workd(ipntr(1)), workd(ipntr(2)) )
!      
!      !  call av (nx, workd(ipntr(1)), workd(ipntr(2)))
!      !  workd(ipntr(2)) = MATRIX*workd(ipntr(1)
!      
!      !           %-----------------------------------------%
!      !           | L O O P   B A C K to call DSAUPD again. |
!      !           %-----------------------------------------%
!       
!      go to 10
!      
!   end if
!   
!   !     %----------------------------------------%
!   !     | Either we have convergence or there is |
!   !     | an error.                              |
!   !     %----------------------------------------%
!   
!   if ( info .lt. 0 ) then
!       
!      !        %--------------------------%
!      !        | Error message. Check the |
!      !        | documentation in DSAUPD. |
!      !        %--------------------------%
!      
!      print *, ' '
!      print *, ' Error with _saupd, info = ', info
!      print *, ' Check documentation in _saupd '
!      print *, ' '
!      
!   else 
!       
!      !        %-------------------------------------------%
!      !        | No fatal errors occurred.                 |
!      !        | Post-Process using DSEUPD.                |
!      !        |                                           |
!      !        | Computed eigenvalues may be extracted.    |  
!      !        |                                           |
!      !        | Eigenvectors may be also computed now if  |
!      !        | desired.  (indicated by rvec = .true.)    | 
!      !        |                                           |
!      !        | The routine DSEUPD now called to do this  |
!      !        | post processing (Other modes may require  |
!      !        | more complicated post processing than     |
!      !        | mode1.)                                   |
!      !        |                                           |
!      !        %-------------------------------------------%
!      
!      rvec = .true.
!      
!      call dseupd ( rvec, 'All', select, d, v, ldv, sigma,      &
!           &         bmat, n, which, nev, tol, resid, ncv, v, ldv,   &
!           &         iparam, ipntr, workd, workl, lworkl, ierr )
!      
!      !         %----------------------------------------------%
!      !         | Eigenvalues are returned in the first column |
!      !         | of the two dimensional array D and the       |
!      !         | corresponding eigenvectors are returned in   |
!      !         | the first NCONV (=IPARAM(5)) columns of the  |
!      !         | two dimensional array V if requested.        |
!      !         | Otherwise, an orthogonal basis for the       |
!      !         | invariant subspace corresponding to the      |
!      !         | eigenvalues in D is returned in V.           |
!      !         %----------------------------------------------%
!      
!      if ( ierr .ne. 0) then
!         
!         !            %------------------------------------%
!         !            | Error condition:                   |
!         !            | Check the documentation of DSEUPD. |
!         !            %------------------------------------%
!          
!         print *, ' '
!         print *, ' Error with _seupd, info = ', ierr
!         print *, ' Check the documentation of _seupd. '
!         print *, ' '
!         
!  !    else
!         
!  !       nconv =  iparam(5)
!  !        do j=1, nconv
!             
!             !               %---------------------------%
!             !               | Compute the residual norm |
!             !               |                           |
!             !               |   ||  A*x - lambda*x ||   |
!             !               |                           |
!             !               | for the NCONV accurately  |
!             !               | computed eigenvalues and  |
!             !               | eigenvectors.  (iparam(5) |
!             !               | indicates how many are    |
!             !               | accurate to the requested |
!             !               | tolerance)                |
!             !               %---------------------------%
!             
! !            call av(dmds,dist, v(1,j), ax)
! !            call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
! !            d(j,2) = dnrm2(n, ax, 1)
! !            d(j,2) = d(j,2) / abs(d(j,1))
!             
! !         end do
!          
!          !            %-----------------------------%
!          !            | Display computed residuals. |
!          !            %-----------------------------%
!          
! !         call dmout(6, nconv, 2, d, maxncv, -6, &
! !              &            'Ritz values and relative residuals')
!      end if
!       
!       !         %-------------------------------------------%
!       !         | Print additional convergence information. |
!       !         %-------------------------------------------%
!       
!      if ( info .eq. 1) then
!         print *, ' '
!         print *, ' Maximum number of iterations reached.'
!         print *, ' '
!      else if ( info .eq. 3) then
!         print *, ' ' 
!         print *, ' No shifts could be applied during implicit',  &
!              &                ' Arnoldi update, try increasing NCV.'
!         print *, ' '
!      end if
!      
! !      print *, ' '
! !      print *, ' _SSIMP '
! !      print *, ' ====== '
! !      print *, ' '
! !      print *, ' Size of the matrix is ', n
! !      print *, ' The number of Ritz values requested is ', nev
! !      print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
! !      print *, ' What portion of the spectrum: ', which
! !      print *, ' The number of converged Ritz values is ', nconv 
! !      print *, ' The number of Implicit Arnoldi update', &
! !           &             ' iterations taken is ', iparam(3)
! !      print *, ' The number of OP*x is ', iparam(9)
! !      print *, ' The convergence criterion is ', tol
! !      print *, ' '
!       
!   end if
! 
! 9000 continue
!    
!   nconv=iparam(5)
!    
! ! output eigenvalues of MDS matrix
! !   print *,'finish diagonalization... start postprocessing...'
!    do jdx=1,nconv
!       eigval(jdx)=sqrt(d(nconv-jdx+1,1))
!    enddo
!    
! ! output MDS coordinates
! !   print *,'output MDS coordinates...'
!   if (status==1) then
!      do jdx=1,dmds
!         do kdx=1,nconv
!            coor(jdx,kdx)=v(jdx,kdx)*sqrt(d(kdx,1))
!         enddo
!      enddo
!   endif
!   return
! end subroutine classicalmds

!subroutine current_time(text)
!  character(8)::date1
!  character(10)::time1
!  character(5)::zone
!  integer::values1(8)
!  character(*)::text
!  call date_and_time(date1,time1,zone,values1)
!  write(6,'(a,a,a)') time1,' : ',text 
!end subroutine current_time
