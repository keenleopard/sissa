real(8) function func(x)
    implicit none
    real(8), intent(in) :: x
    func = x - cos(x)
end function func

program bisection_method
    implicit none
    real(8):: xmin, xmax, delta, e
    real(8):: xold, xnew
    real(8):: g
    integer:: i
    real(8), external :: func

    xmin = 0
    xmax = 3.14/2

    xold = 1
    i=0

    do 
        xnew = (xmin + xmax)/2 
        delta = xnew - xold
        e = abs(delta)/xold
        write(1,*) i, e 
        if (e>=1e-6) then 
            if (func(xmin) * func(xnew) < 0) then 
                xmax = xnew
            else 
                xmin = xnew
            end if 
            xold = xnew
            i = i+1
        else 
            exit 
        end if 
    end do

    print *,"solution is", xnew
    print *, "its cosine is" , cos(xnew)
    print * , "iteration step is ", i
end program bisection_method

    
