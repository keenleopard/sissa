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

program iteration
    implicit none
    integer :: i, j
    real(8) :: r

    call init_random_seed()

    do i = 1, 10
        print '(i4, $)', i
        do j = i+1, 10 
            print '(i4,$)', j
        end do 
        print *, " "
    end do

    do i = 1, 10
        !call random_number(r)
        !print '(i3, $)', int(r*100)+1
        if (i .eq. 5) cycle
        print '(i3, $)', i
    end do
    print *, " "
        
end program 
        
