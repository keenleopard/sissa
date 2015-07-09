program correlation
    implicit none
    !real(8),dimension(424758) :: u(424758)
    real(8), allocatable :: u(:), corr(:)
    integer :: i, dt
    real(8) :: avg, var0
    integer,parameter :: n=10000, dtmax=4000

    allocate(u(n), corr(0:n-1))

    open(unit=15, file='tmp.dat',status='old',action='read')
    read(15,'(f11.6)'), u
    close(15)


    corr = 0.d0
    avg = 0.d0
    do i = 1, n
        avg = avg + u(i)
    end do 
    avg = avg/n

    do dt = 0, dtmax-1
        do i = 1, n/2
            corr(dt) = corr(dt) + (u(i)-avg) * (u(i+dt) - avg)
        end do
    end do

    do dt = 0, dtmax-1
        print *, dt, corr(dt)/(n/2)
    end do
    deallocate(u,corr)


end program 
