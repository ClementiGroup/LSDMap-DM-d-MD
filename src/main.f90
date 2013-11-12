!---------------------------------------------------------------!                                                           
! Function Main routine for LSDMap code                         !
!                                                               !
! Calls p_rmsd_neighbor, p_local_mds, p_wlsdmap                 !
! to calculate a weighted LSDMap in parallel                    ! 
!                                                               !
! Developed by                                                  !                                                          
!   E.Breitmoser, EPCC, Uinversity of Edinburgh                 !                                                           
!---------------------------------------------------------------!                                                           
! Input                                                         !
!                                                               ! 
! 1. See 'ExampleInputFile', created by running                 ! 
!    'prepare_p_lsdmap.sh' for an example input                 !
!    beforehand                                                 !
!    a) If intermediate files should e written:                 !
!       write_rmsd, write_neighbor, write_localscale            !
!    b) Name of the trajectory    nn_traj                       !
!    c) Start and end point ID    ns, ne                        !
!    d) Number of neighbors       Nneigh                        !
!    e) Job ID                    norder                        !
!    f) Number of point for MDS,          dmds                  !
!       start point for the MDS spectra,  kmin                  !
!       the step size for the MDS spectra dk                    !  
!    g) Cutoff for the first derivative of the MDS spectra      !
!                                         seps, deps, neps      !
!    h) Number of CPUs to share the trajectory  ncore           !           
!    i) Status of LSDMap (local scale or constant,              !
!       weighted or non-weighted)       status_dmap, status_eps !
!    j) Weight file name                NN_input_weight         !
!    k) Output file name                output_file             !
!    l) Cutoff for sparse matrix        cutoff                  !
!    m) Local scale file name or constant local scale value     !
!                                       column_eps, eps0        !
! 2. The xyz-file listed in ExampleInputFile                    !
!                                                               !
!                                                               !
! Output                                                        !
!                                                               !
! 1. Writes running info into file                              !
!    (name as set in ExampleBatchScript)                        !
!---------------------------------------------------------------! 
! Develop log                                                   !
!                                                               !
! September 2013 Initial Release                                !
!---------------------------------------------------------------!
! Please reference the papers below if you use LSDMap:          ! 
! 1. Rohrdanz, M.A., Zheng, W., Maggioni, M., and Clementi, C., !
!    J. Chem. Phys., 134, 124116, 2011                          ! 
!    2. Zheng, W., Rohrdanz, M.A., Caflisch, A., Dinner, A.R.,  !
!       and Clementi, C.,                                       !                                                
!    J. Phys. Chem. B, 115, 13065-13074, 2011                   !
!                                                               !                       
! WZ Jul 25 2012                                                !
!---------------------------------------------------------------!

program MainRoutine
  
  use parallel, only : size, rank, ierr, comm, displacements, counts,       &
       trace_exit
  use data, only : ns,ne,Npoints,dim,Nneigh,Natoms,nloc,nstart,nend,        &
       nlimit,extra,traj,idneigh,dist,tmp_rmsd,nn_traj,norder,ncore,dmds,   &
       kmin,dk,neps,seps, deps,NN_input_weight,output_file,status_dmap,     &
       status_eps,column_eps,cutoff,eps0, current_time,write_rmsd,          &
       write_neighbor, write_localscale
  

  implicit none
  include 'mpif.h'
  
  integer :: idx     ! loop counter
  integer :: jdx     ! index counter
  integer :: i       ! loop counter
  integer :: inlen   !
  integer, allocatable :: displacements2(:)      ! displacements relative to receiving buffer for MPI
  
  comm = MPI_COMM_WORLD
  
  call MPI_init(ierr)
  call trace_exit("MPI_Init failed",ierr)
  call MPI_comm_size(comm,size,ierr)
  call trace_exit("MPI_Size failed",ierr)
  call MPI_comm_rank(comm,rank,ierr)
  call trace_exit("MPI_Rank failed",ierr)  

  if(rank==0) write(*,*)'LSDMap v1.1 - Sept 01 2013 - Initial Merger Release'
  if(rank==0) call current_time('Program start...')
  
  ! Defaults for three logical variables which determine, if intermediate IO is needed or not
  write_rmsd = .false.
  write_neighbor = .false.
  write_localscale = .false.
  
  ! Get parameters from ExampleInputFile on processor 0     
  if (rank==0) then
     read(5,*) write_rmsd, write_neighbor, write_localscale
     read(5,*) nn_traj
     read(5,*) ns,ne
     read(5,*) Nneigh
  endif

  ! Broadcast to all processors
  ! Three logical variables, if intermediate IO is needed or not
  call MPI_BCAST(write_rmsd,1,MPI_LOGICAL,0,comm,ierr)
  call trace_exit("MPI_Bcast write_rmsd failed",ierr)
  call MPI_BCAST(write_neighbor,1,MPI_LOGICAL,0,comm,ierr)
  call trace_exit("MPI_Bcast write_neighbor failed",ierr)
  call MPI_BCAST(write_localscale,1,MPI_LOGICAL,0,comm,ierr)
  call trace_exit("MPI_Bcast write_localscale failed",ierr)
  ! Name of the input file
  call MPI_BCAST(nn_traj,200,MPI_CHARACTER,0,comm,ierr)
  call trace_exit("MPI_Bcast nn_traj failed",ierr)
  ! 
  call MPI_BCAST(ns,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast ns failed",ierr)
  call MPI_BCAST(ne,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast ne failed",ierr)
  call MPI_BCAST(Nneigh,1,MPI_INTEGER,0,comm,ierr)
    call trace_exit("MPI_Bcast Nneigh failed",ierr)
  
  ! Split data across processors
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
  
  ! Write status
  !                                                                                                             
  inlen=index(nn_traj,' ') - 1
  if(rank==0) then
     ! Open the xyz-file
     open(10,file=nn_traj,status='old')
     read(10,*) Npoints,dim
     write(*,*)'trajectory read from file ', nn_traj
     write(*,*)'number of points in dataset = ',Npoints
     write(*,*)'original dimension = ',dim
     write(*,*)'number of nearest neighbors = ',Nneigh
     write(*,*)'from ',ns,' to ',ne
     write(*,*)'number of cores=',size
     write(*,*)'points per core=',nloc,'~',nloc+1
     if(write_rmsd) then 
        write(*,*)'rmsd-files will be written.'
     else
        write(*,*)'rmsd-files will NOT be written.'
     endif
     if(write_neighbor) then
        write(*,*)'neighbor-files will be written.'
     else
        write(*,*)'neighbor-files will NOT be written.'
     endif
     if(write_localscale) then
        write(*,*)'localscale-files will be written.'
     else
        write(*,*)'localscale-files will NOT be written.'
     endif
  endif ! Closing if(rank==0), open xyz-file

  ! Broadcast to all processors
  !
  call MPI_BCAST(Npoints,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast Npoints failed",ierr)
  call MPI_BCAST(dim,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast dim failed",ierr)
  Natoms = dim/3
  
  ! Read input trajectory                                                                                                     ! To be used in p_rmsd_neighbor and p_local_mds, better Fortran ordering/mem access
  allocate(traj(3,Npoints*Natoms))
  
  if (rank==0) then
     do idx=1,Npoints
        read(10,*) (traj(1,(idx-1)*Natoms+jdx),traj(2,(idx-1)*Natoms+jdx),traj(3,(idx-1)*Natoms+jdx), jdx=1,Natoms)
     enddo
     close(10)
  endif
  ! Broadcast to all processors 
  !
  call MPI_BCAST(traj,Npoints*dim,MPI_REAL,0,comm,ierr)
    call trace_exit("MPI_Bcast traj failed",ierr)
  
  ! Set-up for call to first subroutine
  allocate(dist(Npoints,nstart:nend))
  allocate(tmp_rmsd(Npoints,Npoints))
  allocate(idneigh(Npoints,nloc))
  allocate(counts(0:size-1))
  allocate(displacements(0:size-1))

  ! All processors must know how many elements of nloc are on each processor
  ! This info is gathered in the array counts  
  call MPI_ALLGATHER(nloc,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
  call trace_exit("MPI_Allgather of nloc failed",ierr)  

  displacements(0) = 0
  do i = 1, size-1
     displacements(i) = displacements(i-1) + counts(i-1)*Npoints
  enddo
  
  ! Call first subroutine p_rmsd_neighbor.f90
  call p_rmsd_neighbor
  
  ! Read-in for call to second subroutine
  ! Continue reading input parameters from ExampleInputFile
  if (rank==0) then
     read(5,*) norder
     read(5,*) dmds,kmin,dk
     read(5,*) seps,deps,neps
     read(5,*) ncore
  endif
  ! Broadcast to all processors                             
  !   
  call MPI_BCAST(norder,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast norder failed",ierr)
  call MPI_BCAST(dmds,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast dmds failed",ierr)
  call MPI_BCAST(kmin,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast kmin failed",ierr)
  call MPI_BCAST(dk,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast dk failed",ierr)
  call MPI_BCAST(seps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call trace_exit("MPI_Bcast seps failed",ierr)
  call MPI_BCAST(deps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call trace_exit("MPI_Bcast deps failed",ierr)
  call MPI_BCAST(neps,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast neps failed",ierr)
  call MPI_BCAST(ncore,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast ncore failed",ierr)
  
  ! Set-up for call to second subroutine
  ! Call second subroutine p_local_mds.f90
  call p_local_mds
  
  ! Read-in data for third subroutine
  ! Still from ExampleInputFile
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
  ! Broadcast to all processors                                 
  !   
  call MPI_BCAST(output_file,200,MPI_CHARACTER,0,comm,ierr)
  call trace_exit("MPI_Bcast output_file failed",ierr)
  call MPI_BCAST(cutoff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call trace_exit("MPI_Bcast cutoff failed",ierr)
  call MPI_BCAST(status_dmap,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast status_dmap failed",ierr)
  call MPI_BCAST(status_eps,1,MPI_INTEGER,0,comm,ierr)
  call trace_exit("MPI_Bcast status_eps failed",ierr)

  if(status_eps==0) then 
     call MPI_BCAST(eps0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
     call trace_exit("MPI_Bcast eps0 failed",ierr)
  endif

  ! Call third subroutine p_wlsdmap.f90
  call Weighted_LSDMap
  
  deallocate(dist)
  deallocate(tmp_rmsd)
  deallocate(idneigh)
  deallocate(counts)
  deallocate(displacements)
  deallocate(traj)
  
  if(rank==0) call current_time('Program end.')
  
  call MPI_FINALIZE(ierr)
  call trace_exit("MPI_Finalize failed",ierr)  

end program MainRoutine
