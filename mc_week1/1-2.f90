subroutine init_random_seed()

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
end subroutine

real(8) function func(x)
    implicit none
    real(8), intent(in) :: x

    func = 1.d0/2.d0/sqrt(x)
end function 

program Mean_Intrgration
    implicit none
    integer :: n, i 
    real(8) :: f, f2, mean, var, r
    real(8),external :: func

999 continue
    call init_random_seed()

    print*, "Input the number of random numbers: "
    read *, n

    f = 0
    f2 = 0
    mean = 0
    var = 0

    do i = 1, n
        call random_number(r)
        f = f + func(r)
        f2 = f2 + (func(r))**2
    end do
    mean = f/n
    var = f2/n - mean**2

    print *, mean
go to 999
end program
