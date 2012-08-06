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
! Split RMSD binary file
!
! Input
! 1. RMSD binary file name
! 2. Number of points in the original trajectory
! 3. Number of points to add for each loop
! 4. Number of CPUs for the calculation
! 5. Output file name prefix
! See "split_rmsd.input_example" for an example of input.
!
! Output
! RMSD binary file for p_wlsdmap_embed
! 
! Please reference the papers below if you use LSDMap:
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., 
!    J. Chem. Phys., 134, 124116, 2011
! 2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R., and Clementi, C.,
!    J. Phys. Chem. B, 115, 13065-13074, 2011
!
! WZ Jul 25 2012
!-----------------------------------------------------------------------

program split_rmsd

implicit none

integer :: idx,jdx,kdx,npoints,inlen,outlen,nproc,tmp_int,nadd,npointsnew

integer,allocatable :: idxlocpoint(:)

real,allocatable :: tmp_rmsd(:)

character(200) :: nn_rmsd, nn_output_prefix, nn_output

logical :: lexist

print *,'LSDMap v1.0 - Aug 01 2012 - Initial Release'

read(5,*) nn_rmsd
read(5,*) npoints
read(5,*) nadd
read(5,*) nproc
read(5,*) nn_output_prefix

inlen=index(nn_rmsd,' ')-1
outlen=index(nn_output_prefix,' ')-1

call system('rm rmsd/'//nn_output_prefix(1:outlen)//'_*')
inquire(file='rmsd/'//nn_rmsd(1:inlen)//'_all',exist=lexist)
if(.not. lexist) then
   print *,'merge rmsd file...'
   call system('cat rmsd/'//nn_rmsd(1:inlen)//'_9* >rmsd/'//nn_rmsd(1:inlen)//'_all')
endif

npointsnew=npoints+nadd
allocate(idxlocpoint(nproc+1))
tmp_int=int(npointsnew/nproc)
idxlocpoint(1)=1
do idx=2,nproc+1
   if(idx<=mod(npointsnew,nproc)+1)then
      idxlocpoint(idx)=idxlocpoint(idx-1)+tmp_int+1
   else
      idxlocpoint(idx)=idxlocpoint(idx-1)+tmp_int
   endif
   if (idxlocpoint(idx)>npoints) then
      tmp_int=idx
      idxlocpoint(idx)=npoints+1
      exit
   endif
enddo

allocate(tmp_rmsd(npoints))
open(10,file='rmsd/'//nn_rmsd(1:inlen)//'_all',form='unformatted',status='old',access='direct',recl=npoints*4)

do idx=1,tmp_int-1
   write(nn_output,'(a,a,i7,a,i7)') nn_output_prefix(1:outlen),'_',9000000+idxlocpoint(idx),'_',9000000+idxlocpoint(idx+1)-1
   open(50,file='rmsd/'//nn_output,form='unformatted',status='unknown',access='direct',recl=npoints*4)
   do jdx=idxlocpoint(idx),idxlocpoint(idx+1)-1
      read(10,rec=jdx) (tmp_rmsd(kdx),kdx=1,npoints)
      write(50,rec=jdx-idxlocpoint(idx)+1) (tmp_rmsd(kdx),kdx=1,npoints)
   enddo
   close(50)
enddo
close(10)


end program split_rmsd
