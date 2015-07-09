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

real(8) function f(x)
    implicit none
    real(8), intent(in) :: x

    f = 2.d0*sqrt(1.d0-x**2)/3.1415926
end function

    
program rejection_method
    implicit none
    real(8) :: r1, r2, x, y
    integer :: i, histo(20)
    real(8), external :: f
    real(8), parameter :: pi = 3.1415926d0
    

    call init_random_seed()
    histo = 0
    open(unit=1, file="points.dat", status="unknown", position="append")

    do i=1,1000000
        call random_number(r1)
        call random_number(r2)
        x = 2.d0 * r1 -1.d0
        y = 2.d0 * r2/pi
        if (y<= f(x)) then
            write(1, '(f11.4, f11.4)'), x, y
            if ( -1.d0<= x .and. x <= (-1.d0+0.1)) histo(1) = histo(1) + 1
            if ( (-1.d0+0.1)< x .and. x <= (-1.d0+0.1*2)) histo(2)=histo(2) +1
            if ( (-1.d0+0.1*2)< x .and. x <= (-1.d0+0.1*3)) histo(3)=histo(3)+1
            if ( (-1.d0+0.1*3)< x .and. x <= (-1.d0+0.1*4)) histo(4)=histo(4)+1
            if ( (-1.d0+0.1*4)< x  .and. x<= (-1.d0+0.1*5)) histo(5)=histo(5)+1
            if ( (-1.d0+0.1*5)< x  .and. x<= (-1.d0+0.1*6)) histo(6)=histo(6)+1
            if ( (-1.d0+0.1*6)< x  .and. x<= (-1.d0+0.1*7)) histo(7)=histo(7)+1
            if ( (-1.d0+0.1*7)< x  .and. x<= (-1.d0+0.1*8)) histo(8)=histo(8)+1
            if ( (-1.d0+0.1*8)< x  .and. x<= (-1.d0+0.1*9)) histo(9)=histo(9)+1
            if ( (-1.d0+0.1*9)< x  .and. x<= (-1.d0+0.1*10)) histo(10) =histo(10) +1
            if ( (-1.d0+0.1*10)< x .and. x<= (-1.d0+0.1*11)) histo(11) = histo(11)+1
            if ( (-1.d0+0.1*11)< x .and. x<= (-1.d0+0.1*12)) histo(12) = histo(12)+1
            if ( (-1.d0+0.1*12)< x .and. x<= (-1.d0+0.1*13)) histo(13) = histo(13)+1
            if ( (-1.d0+0.1*13)< x .and. x<= (-1.d0+0.1*14)) histo(14) = histo(14)+1
            if ( (-1.d0+0.1*14)< x .and. x<= (-1.d0+0.1*15)) histo(15) = histo(15)+1
            if ( (-1.d0+0.1*15)< x .and. x<= (-1.d0+0.1*16)) histo(16) = histo(16)+1
            if ( (-1.d0+0.1*16)< x .and. x<= (-1.d0+0.1*17)) histo(17) = histo(17)+1
            if ( (-1.d0+0.1*17)< x .and. x<= (-1.d0+0.1*18)) histo(18) = histo(18)+1
            if ( (-1.d0+0.1*18)< x .and. x<= (-1.d0+0.1*19)) histo(19) = histo(19)+1
            if ( (-1.d0+0.1*19)< x .and. x<= (-1.d0+0.1*20)) histo(20) = histo(20)+1
        else 
            continue
        end if
    end do

    do i = 0,19
        print '(f11.2 f11.2, i9)', -1.d0+0.1*i, -1.d0 + 0.1*(i+1), histo(i+1)
    end do
    close(1)
end program 


