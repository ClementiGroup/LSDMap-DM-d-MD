! From rosettacode.org/wiki/Sorting_algorithms/Quicksort
MODULE Qsort_Module
 
IMPLICIT NONE
 
CONTAINS
 
RECURSIVE SUBROUTINE Qsort(a, id)
  INTEGER,parameter :: REAL_KIND = 8
  REAL(KIND=REAL_KIND), INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(IN OUT) :: id(:)
  INTEGER :: left, right
  IF(size(a) > 1) THEN
     CALL Partition(a, id, left, right)
     CALL Qsort(a(:right), id(:right))
     CALL Qsort(a(left:), id(left:))
  END IF
 
END SUBROUTINE Qsort
 
SUBROUTINE Partition(a, id, left, right)
  INTEGER,parameter :: REAL_KIND = 8
  REAL(KIND=REAL_KIND), INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(IN OUT) :: id(:)
  INTEGER, INTENT(OUT) :: left, right
  REAL(KIND=REAL_KIND) :: tmp_real, pivot
  INTEGER :: tmp_int
  INTEGER :: idx

  pivot = a(int(size(a)/2.))
  left = 1                         
  right = size(a)
  
  DO WHILE (left <= right)
     DO WHILE (a(right) > pivot)
        right = right-1
     END DO
     DO WHILE (a(left) < pivot)
        left = left + 1
     END DO
     IF (left <= right) THEN 
        tmp_real = a(left)
        a(left) = a(right)
        a(right) = tmp_real
        tmp_int = id(left)
        id(left) = id(right)
        id(right) = tmp_int
        left = left + 1
        right = right - 1
     END IF
  END DO
  RETURN
END SUBROUTINE Partition

RECURSIVE SUBROUTINE iQsort(a, id)
  INTEGER, INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(IN OUT) :: id(:)
  INTEGER :: split
  IF(size(a) > 1) THEN
     CALL iPartition(a, id, split)
     CALL iQsort(a(:split-1), id(:split-1))
     CALL iQsort(a(split:), id(split:))
  END IF
 
END SUBROUTINE iQsort
 
SUBROUTINE iPartition(a, id, marker)
  INTEGER, INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(IN OUT) :: id(:)
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: tmp_real, pivot
  INTEGER :: left, right, tmp_int
  INTEGER :: idx

  pivot = (a(1) + a(size(a))) / 2.  ! Average of first and last elements to prevent quadratic 
  left = 0                         ! behavior with sorted or reverse sorted data
  right = size(a) + 1
  
  DO WHILE (left < right)
     right = right - 1
     DO WHILE (a(right) > pivot)
        right = right-1
     END DO
     left = left + 1
     DO WHILE (a(left) < pivot)
        left = left + 1
     END DO
     IF (left < right) THEN 
        tmp_real = a(left)
        a(left) = a(right)
        a(right) = tmp_real
        tmp_int = id(left)
        id(left) = id(right)
        id(right) = tmp_int
     END IF
  END DO


  IF (left == right) THEN
     marker = left + 1
  ELSE
     marker = left
  END IF

END SUBROUTINE iPartition
 
END MODULE Qsort_Module
