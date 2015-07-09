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

subroutine walker(N, pos, D)
    implicit none
    integer, intent(in) :: N
    integer, intent(out) :: pos(0:N), D

    integer :: i
    real(8) :: r

    pos = 0
    do i = 1, N
        call random_number(r)
        print *, r
        if (r > 0.5) then
            pos(i) = pos (i-1) + 1 ! moving right one step
        else
            pos(i) = pos (i-1) - 1 ! moving left one step
        end if
    end do
    D = pos(N)
    print *, pos
end subroutine walker

subroutine average(n, M, Dn, D2n)
    implicit none
    integer, intent(in) ::  n, M
    real(8), intent(out) :: Dn, D2n ! the average of D(n) and D^2 (n) over M walkers

    integer :: i
    real(8) :: s, s2

    integer :: pos(0:n)
    integer :: D_single

    s = 0
    s2 = 0
    do i = 1, M
        call walker (n, pos , D_single)
        s = s + D_single
        s2 = s2 + D_single**2
    end do
    Dn = s/M
    D2n = s2/M
end subroutine average

subroutine histogram (phi, M, N)
    implicit none
    integer, intent(in) :: M, N
    integer, intent(inout) :: phi(-N:N)

    integer, allocatable :: pos(:)
    integer :: D
    integer :: i

    allocate(pos(0:N))
    phi = 0

    do i = 1, M ! loop over M walkers
        call walker (N, pos, D)
        phi(D) = phi (D) + 1
    end do
    deallocate(pos)

    open (unit=23, file="walkers.dat", status="unknown")
    do i = -N, N, 2
        write (23, *), i, phi(i)
    end do
    close(23)

    print *, phi
end subroutine histogram


program main
    implicit none
    real(8) :: Dn, D2n
    integer :: nwalker, nstep
    integer, allocatable :: phi (:)
    call init_random_seed()

    print *, "Input the number of Walkers:"
    read *, nwalker

    print *, "Input the number of Steps:"
    read *, nstep

    allocate(phi(-nstep:nstep))

    !call average(nstep, nwalker, Dn, D2n)
    !print *, Dn
    !print *, D2n

    call histogram (phi, nwalker, nstep)
    deallocate(phi)
end program main
