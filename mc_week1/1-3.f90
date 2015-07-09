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

program non_uniform
!   Non-uniform distribution
    implicit none

    integer :: n, i, iter
    real(8), allocatable :: r(:)
    real(8) :: psi(10000), psi2(10000)
    real(8) :: psi_avg, psi2_avg,  var

    
999 continue
    psi = 0.d0
    psi2 = 0.d0
    psi_avg = 0.d0
    psi2_avg = 0.d0
    var = 0.d0

    call init_random_seed()
    print *, "Input the number of random numbers: "
    read *, n
    allocate (r(n))
    
    do iter = 1, 10000 
        do i = 1, n
            call random_number(r(i))
            psi(iter) = psi(iter) + r(i)**2
        end do

    psi(iter) = psi(iter)/n
    psi2(iter) = psi(iter)**2
    end do

    deallocate(r)

    do iter = 1, 10000
        psi_avg = psi_avg + psi(iter)
        psi2_avg = psi2_avg + psi2(iter)
    end do
    psi_avg = psi_avg/10000
    psi2_avg = psi2_avg/10000
    var = psi2_avg - psi_avg**2

    print *, psi_avg, var
    
    print *, var*45*n/4

    go to 999
    close(11)
end program 

        

    
