program correlation
    implicit none
    !real(8),dimension(424758) :: u(424758)
    real(8) :: u(400000), corr(400)
    integer :: dt, nbinn
    real(8) :: avg0, var0

    open(unit=15, file='tmp.dat',status='old',action='read')
    read(15,'(f11.6)'), u
    close(15)


    corr = 0.d0
    avg0 = 0.d0
    do nbinn = 0, 999
        avg0 = avg0 + u(1+nbinn*400)
        do dt = 1, 400
            corr(dt) = corr(dt) + u(dt+nbinn*400) * u(1+nbinn*400)
        end do 
    end do
    avg0 = avg0 / 1000.d0
    corr = corr / 1000.d0

    corr = corr - avg0**2

    do dt = 1,400
        print *, dt, corr(dt)
    end do


end program 
