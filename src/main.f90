!---------------------------------------------------------------!                                                           
! LSDMap v1.1 - Sept 2013 - Merger Release                      !                                                           
!                                                               !                                                          
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uonversity of Edinburgh                 !                                                           
!---------------------------------------------------------------!                                                           
program MainRoutine

use parallel, only : size, rank, ierr, comm, displacements, counts
use data, only : ns,ne,Npoints,dim,Nneigh,Natoms,nloc,nstart,nend,nlimit,extra,traj,    &
         idneigh,dist,tmp_rmsd,nn_traj,norder,ncore,dmds,kmin,dk,neps,seps,   &
          deps,NN_input_weight,output_file,status_dmap,status_eps,column_eps,cutoff,eps0, current_time, &
          write_rmsd, write_neighbor, write_localscale


implicit none
include 'mpif.h'

integer :: idx, jdx,i
integer :: inlen
integer, allocatable :: displacements2(:)
!real, allocatable :: traj(:,:)
integer, dimension(1:2) :: order2 = (/2,1/)

comm = MPI_COMM_WORLD

call MPI_init(ierr)
call MPI_comm_size(comm,size,ierr)
call MPI_comm_rank(comm,rank,ierr)

if(rank==0) print *,'LSDMap v1.1 - Sept 01 2013 - Initial Merger Release'
if(rank==0) call current_time('Program start...')

write_rmsd = .false.
write_neighbor = .false.
write_localscale = .false.

! get parameters from xyz-input file                                                                                           
if (rank==0) then
   read(5,*) write_rmsd, write_neighbor, write_localscale
   read(5,*) nn_traj
   read(5,*) ns,ne
   read(5,*) Nneigh
endif
call MPI_BCAST(write_rmsd,1,MPI_LOGICAL,0,comm,ierr)
call MPI_BCAST(write_neighbor,1,MPI_LOGICAL,0,comm,ierr)
call MPI_BCAST(write_localscale,1,MPI_LOGICAL,0,comm,ierr)
call MPI_BCAST(nn_traj,200,MPI_CHARACTER,0,comm,ierr)
call MPI_BCAST(ns,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(ne,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(Nneigh,1,MPI_INTEGER,0,comm,ierr)


! split data                                                                                                               
nloc=int((ne-ns+1)/size)

if(rank<mod(ne-ns+1,size)) then
   nloc=nloc+1
   extra=0
else
   extra=mod(ne-ns+1,size)
endif
nstart=ns+rank*nloc+extra
nend=nstart+nloc-1

! For second subroutine
nlimit = ceiling(1.*(ne-ns+1)/size)

! write status                                                                                                             
inlen=index(nn_traj,' ') - 1
if(rank==0) then
   open(10,file=nn_traj,status='old')
   read(10,*) Npoints,dim
   print *,'trajectory read from file ', nn_traj
   print *,'number of points in dataset = ',Npoints
   print *,'original dimension = ',dim
   print *,'number of nearest neighbors = ',Nneigh
   print *,'from ',ns,' to ',ne
   print *,'number of cores=',size
   print *,'points per core=',nloc,'~',nloc+1
   if(write_rmsd) then 
      print *,'rmsd-files will be written.'
   else
      print *,'rmsd-files will NOT be written.'
   endif
   if(write_neighbor) then
      print *,'neighbor-files will be written.'
   else
      print *,'neighbor-files will NOT be written.'
   endif
   if(write_localscale) then
      print *,'localscale-files will be written.'
   else
      print *,'localscale-files will NOT be written.'
   endif
endif
call MPI_BCAST(Npoints,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(dim,1,MPI_INTEGER,0,comm,ierr)
Natoms = dim/3

! read input trajectory                                                                                                    
! To be used in p_rmsd_neighbor and p_local_mds, better Fortran ordering/mem access
allocate(traj(3,Npoints*Natoms))

if (rank==0) then
   do idx=1,Npoints
      read(10,*) (traj(1,(idx-1)*Natoms+jdx),traj(2,(idx-1)*Natoms+jdx),traj(3,(idx-1)*Natoms+jdx), jdx=1,Natoms)
   enddo
   close(10)
endif
call MPI_BCAST(traj,Npoints*dim,MPI_REAL,0,comm,ierr)


! Set-up for call to first subroutine
allocate(dist(Npoints,nstart:nend))
allocate(tmp_rmsd(Npoints,Npoints))
allocate(idneigh(Npoints,nloc))
allocate(counts(0:size-1))
allocate(displacements(0:size-1))

call MPI_ALLGATHER(nloc,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)

 displacements(0) = 0
 do i = 1, size-1
  displacements(i) = displacements(i-1) + counts(i-1)*Npoints
 enddo

! call first subroutine p_rmsd_neighbor.f90
call p_rmsd_neighbor

! Read-in for call to second subroutine
! Continue reading input parameters                                                                                                    
if (rank==0) then
   read(5,*) norder
   read(5,*) dmds,kmin,dk
   read(5,*) seps,deps,neps
   read(5,*) ncore
endif
call MPI_BCAST(norder,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(dmds,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(kmin,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(dk,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(seps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_BCAST(deps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_BCAST(neps,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(ncore,1,MPI_INTEGER,0,comm,ierr)


! Set-up for call to second subroutine
! call second subroutine p_local_mds.f90
call p_local_mds

! Read-in data for third subroutine
if (rank==0) then
   read(5,*) status_dmap,status_eps
   if(status_dmap==1) then
      read(5,*) NN_input_weight
   else
      read(5,*)
   endif
   read(5,'(a)') output_file
   read(5,*) cutoff
   if(status_eps==1) then
      read(5,*) column_eps
   else
      read(5,*) eps0
   endif
endif
call MPI_BCAST(output_file,200,MPI_CHARACTER,0,comm,ierr)
call MPI_BCAST(cutoff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
call MPI_BCAST(status_dmap,1,MPI_INTEGER,0,comm,ierr)
call MPI_BCAST(status_eps,1,MPI_INTEGER,0,comm,ierr)
if(status_eps==0) call MPI_BCAST(eps0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

! call third subroutine p_wlsdmap.f90
call Weighted_LSDMap

deallocate(dist)
deallocate(tmp_rmsd)
deallocate(idneigh)
deallocate(counts)
deallocate(displacements)
deallocate(traj)

if(rank==0) call current_time('Program end.')

call MPI_FINALIZE(ierr)

end program MainRoutine
