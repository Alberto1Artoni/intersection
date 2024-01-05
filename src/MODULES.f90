module globalVariables

    ! collects the global variables
    ! mpi variables

    implicit none

    type mapIntegers               ! I am going to use it as a map for the intersections of elems
        integer*4 , dimension(:), allocatable :: map
    end type

    integer*4 :: mpi_id, mpi_np, mpi_comm, mpi_ierr
    
end module


module QSORT_C_MODULE

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  integer, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module QSORT_C_MODULE
