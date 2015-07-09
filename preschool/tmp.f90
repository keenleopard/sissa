program test
    real(8) :: du, a
    du = 10715927240.885038-3467497187.6015377
    a = exp(-1.d0*du/2.d0)
    print *, du, a
end program
