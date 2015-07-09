program newton_method
    implicit none
    real(8):: xmin, xmax, delta, e
    real(8):: x0, x1
    real(8):: g
    integer:: i

    i = 0
    x0 = 0.5
    do
        x1 = -(x0-cos(x0)) / (1+sin(x0)) + x0
        delta = x1-x0
        e = abs(delta)/x0
        if (e>=1e-6) then
            x0 = x1
            !print *, x0
            i = i + 1
        else 
            exit
        end if
    end do
    print *,"solution is", x0 
    print *, "its cosine is" , cos(x0)
    print * , "iteration step is ", i
end program newton_method


