module para_module
    implicit none
    save
    real(8), parameter :: eps12 =1.d-12
    real(8),parameter  :: rho=0.5, T=2.d0, P=3.0
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
! Compute the energy difference before and after moving one particle.
! The volume is not changed, so we only need to caculate du with 
! respect to selected particle for faster computation. 
    use para_module
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
    do i = 1,n                          !the cycle to compute u1
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

    do i = 1,n                          !the cycle to compute u2
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

subroutine energy_diff_v2(r1, r2, n, du, rc, L1, L2)
! Compute the energy difference before and after CHANGING L.
! The energy difference for npt ensemble now takes a more complex formula
! du = unew - uold
!    = U(s,V') - U(s,V) + P(V'-V) - N*ln(V'/V)/beta
! where s is recaled unit ri = L si so si always lies in range [0,1].
! 
! Since the distances between particles are rescaled, we need to compute
! the difference between TOTAL energies. 
    use para_module
    implicit none
    real(8), intent(in) :: r1(3,n), r2(3, n), rc, L1, L2
    integer, intent(in) :: n 
    real(8), intent(out) :: du
    real(8) :: u1, u2, rc2, v1, v2

    u1 = 0.d0
    u2 = 0.d0
    du = 0.d0
    rc2 = rc**2
    v1 = L1**3
    v2 = L2**3
    call lennard_jones(r1, n, u1, rc, L1)
    call lennard_jones(r1, n, u2, rc, L2)

    du = (u2 - u1) + P*(v2-v1) - n*T*log(v2/v1)
end subroutine

subroutine trial_move(n,r_old, r_new, L_old, L_new, mv )
! this subroutine is to propose a trial partical move OR a volume move
! on average 1 out of 100 total trial moves is volume
! and on average the prob(left) = prob(right) for volume move
    implicit none
    integer, intent(in) ::n
    real(8), intent(in) :: r_old(3,n), L_old
    real(8), intent(out) :: r_new(3,n), L_new
    integer, intent(out) :: mv
    real(8), parameter :: ratio=1.d0/100.d0
    real(8) :: rand_choice, rand_lf, rand_L, rand_mv, rand(3), dl, delta
    integer :: dimen

    dl = 0.005*L_old
    delta=0.1*L_old
    r_new = r_old 
    L_new = L_old
    mv = 0
    call random_number(rand_choice)
    if (rand_choice < ratio) then       !move volume here
        mv = -1
        call random_number(rand_lf) 
        call random_number(rand_L)
        if (rand_lf < 0.5) then         !decrease volume
            L_new = L_old - dl*rand_L
        else                            !increase volume
            L_new = L_old + dl*rand_L
        end if 
        r_new = r_old/L_old * L_new     !rescale positions here
    else                                !move one particle here
        call random_number(rand_mv)
        mv = int(rand_mv * n) + 1
        do dimen = 1,3
            call random_number(rand(dimen))
            r_new(dimen,mv) = r_new(dimen,mv) + delta*(rand(dimen) - 0.5)
            if (r_new(dimen,mv) .lt. 0) r_new(dimen,mv) = r_new(dimen,mv) +L_old
            if (r_new(dimen,mv) .gt. L_old) r_new(dimen,mv) = r_new(dimen,mv) -L_old
        end do
    end if
end subroutine trial_move

subroutine metropolis
    use para_module
    implicit none
    integer,parameter :: n=100 !number of particles
    real(8) :: r_old(3,n), r_new(3,n)
    real(8) :: initial(3), acc, rand_acc
    real(8) :: u, du
    integer :: mv, iter, dimen, i, accepted
    real(8) :: L, Lnew, rc, acceptance, density
    !logical :: newposi

    call init_random_seed()

    acc = 0.d0
    L=(n/rho)**(1.d0/3.d0)              ! initial size of box L=(n/rho)^(1/3)
    rc = L/2.d0                         ! cutoff of interaction
    accepted = 0
    acceptance = 0.d0
    density = 0.d0
    !newposi = .true.

    u = 0.d0
    du = 0.d0

    do i = 1, n
        do dimen = 1, 3
            !r_old(dimen,i) = L/n * i   ! evenly distributed initial condition
            call random_number(initial(dimen))
            r_old(dimen,i) = initial(dimen) * L    ! randomly distributed IC
        end do
    end do
    r_new = r_old
    

    do iter = 1, 100000
        call trial_move(n, r_old, r_new, L, Lnew, mv)
        if (mv .eq. -1) then            ! volume is trial changed 
            call energy_diff_v2(r_old, r_new, n, du, rc, L, Lnew)
        elseif (0 < mv .and. mv < (n+1)) then       !particle 'mv' is moved
            call energy_diff(r_old, r_new(:,mv), mv, n, du, rc, L)
        else 
            print *, "Something is wrong in trial move"
        end if

        acc = min(1.d0, exp(-1.d0*du/T))
        call random_number(rand_acc)
        if (rand_acc < acc) then                    !accept the move
            r_old = r_new
            L = Lnew
            call lennard_jones(r_old, n, u, rc, L)
            accepted = accepted +1
            acceptance = dfloat(accepted) / iter
            density = N/ (L**3)
            print *, accepted,iter, acceptance, "u=", u, "du= ", du, "rho=", density
        else 
            !print *, "rejected"
        end if
    end do
end subroutine

program main
    implicit none
    !real(8) :: energy(2000)
    !integer :: i
    call metropolis

    !do i = 1, 2000
    !    print * , energy(i)
    !end do
end program
    

