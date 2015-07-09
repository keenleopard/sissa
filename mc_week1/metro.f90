module para_module
    implicit none
    save
    real(8), parameter :: eps12 =1.d-12
    real(8),parameter  :: rho=0.5, T=2.d0
end module

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

subroutine lennard_jones(pts,n, e, rc, L)
! compute the total energy for a given configuration
    use para_module
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: pts(3,n), rc, L
    real(8), intent(out) :: e 
    real(8) :: r2, rc2, dx(3)
    integer :: i, j, dimen
    
    e = 0.d0
    rc2 = rc**2
    do i = 1, n
        do j = i+1, n
            r2 = 0.d0
            do dimen = 1, 3
                dx(dimen) = abs(pts(dimen,i) - pts(dimen,j))
                dx(dimen) = min(dx(dimen), L-dx(dimen))
                r2 = r2 + dx(dimen)**2
            end do
            if (r2<eps12) r2=eps12
            !if (r2<= rc2) e = e + ((1/r2)**6 - (1/r2)**3)
            e = e + ((1/r2)**6 - (1/r2)**3)
        end do
    end do
    e = 4*e
end subroutine

subroutine energy_diff(r_all, r_move, mv, n, du, rc, L)
    use para_module
! compute the energy difference before and after moving one particle
    implicit none
    real(8), intent(in) :: r_move(3), r_all(3, n), rc, L
    integer, intent(in) :: n, mv
    real(8), intent(out) :: du
    real(8) :: u1, u2, r2, dx(3), rc2
    integer :: i, dimen

    u1 = 0.d0
    u2 = 0.d0
    du = 0.d0
    rc2 = rc**2
    do i = 1,n
        r2 = 0.d0
        if (i .eq. mv) cycle
        do dimen = 1, 3
            dx(dimen) = abs(r_all(dimen, i) - r_all(dimen, mv))
            dx(dimen) = min(dx(dimen), L-dx(dimen))
            !print *, i, "dx",dimen,dx(dimen)
            r2 = r2 + dx(dimen)**2
        end do
        if( r2 < eps12) r2=eps12
        u1 = u1 + ((1/r2)**6 - (1/r2)**3)
        !print *, "r2",r2
    end do

    do i = 1,n
        r2 = 0.d0
        if (i .eq. mv) cycle
        do dimen = 1, 3
            dx(dimen) = abs(r_all(dimen, i) - r_move(dimen))
            dx(dimen) = min(dx(dimen), L-dx(dimen))
            r2 = r2 + dx(dimen)**2
        end do
        if( r2 < eps12) r2=eps12
        u2 = u2 + ((1/r2)**6 - (1/r2)**3)
    end do
    du = 4.d0 * (u2 - u1)
end subroutine

program metropolis
    use para_module
    implicit none
    integer,parameter :: n=100 !number of particles
    real(8) :: r_old(3,n), r_new(3,n)
    real(8) :: rand(3), acc, rand_acc, rand_mv
    real(8) :: del, uold, unew, du
    integer :: mv, iter, dimen, i, accepted
    real(8) :: L, rc
    !logical :: newposi

    call init_random_seed()

    acc = 0.d0
    L=(n/rho)**(1.d0/3.d0) ! Size of box L=(n/rho)^(1/3)
    del = 0.1*L
    rc = L/2.d0    !cutoff

    do i = 1, n
        do dimen = 1, 3
            !r_old(dimen,i) = L/n * i   ! evenly distributed initial condition
            call random_number(rand(dimen))
            r_old(dimen,i) = rand(dimen) * L
        end do
    end do
    
    !newposi = .true.
    accepted =0

    do iter = 1, 1000000
        call random_number(rand_mv)
        rand_mv = rand_mv * n
        mv = int(rand_mv) + 1

        do dimen = 1, 3
            call random_number(rand(dimen))
            r_new(dimen, mv) = r_old(dimen,mv) + del * (rand(dimen) -0.5)
            if (r_new(dimen,mv) .lt. 0.d0) then
                r_new(dimen,mv) = r_new(dimen,mv) + L
            elseif (r_new(dimen,mv) .gt. L) then
                r_new(dimen,mv) = r_new(dimen,mv) - L
            end if
        end do

        !unew = 0.d0
        !call lennard_jones(r_new,n, unew, rc, L)
        !print *, mv, r_new(:,mv)
        call energy_diff(r_old, r_new(:,mv), mv, n, du, rc, L)
        acc = min(1.d0, exp(-1.d0*du/T))
        call random_number(rand_acc)
        if (rand_acc < acc) then
            r_old(:,mv) = r_new(:,mv)
            call lennard_jones(r_old, n, uold, rc, L)
            accepted = accepted +1
            print *, accepted,iter," u=", uold, " acc=", acc, " du= ", du
        else 
            !print *, "rejected"
        end if
    end do
end program
