program test
    implicit none
    integer(8), allocatable :: a(:)
    integer i

    allocate(a(20))
    do i = 1, 10
        a(i) = i*2
    end do
    print *, a
end program
