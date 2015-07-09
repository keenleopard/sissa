program iteration_method
    implicit none
    real(8):: xmin, xmax, delta, e
    real(8):: x0, x1
    integer:: i

    xmin = 0
    xmax = 1

    x0 = 0.9
    i=0

    
    do 
        x1 = cos(x0)
        delta = x1 - x0
        e = abs(delta)/x0 
        if ( e >= 1e-6) then
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
end program iteration_method
        
            


