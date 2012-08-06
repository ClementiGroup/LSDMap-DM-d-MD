module ftn_c
  interface
  real(c_double) function calcrmsdrotationalmatrix (Natoms,x1,x2,x3,&
       y1,y2,y3,rot,weight) BIND(C,name='calcrmsdrotationalmatrix')
    use iso_c_binding
    implicit none
    integer,parameter :: real_kind=8
    integer(c_int) :: Natoms
    real(c_double) :: x1(Natoms,3)
    real(c_double) :: x2(Natoms,3)
    real(c_double) :: x3(Natoms,3)
    real(c_double) :: y1(Natoms,3)
    real(c_double) :: y2(Natoms,3)
    real(c_double) :: y3(Natoms,3)
    real(c_double) :: rot(9)
    real(c_double) :: weight(Natoms)
  end function calcrmsdrotationalmatrix
  end interface
end module ftn_c
