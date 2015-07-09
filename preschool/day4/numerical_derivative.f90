program numerical_derivative
!----------------------------------------------------------------------------
   implicit none

   real(8), parameter :: pi = 3.14159265358979323846d0

   real(8) ::  x, dx, fx, fp, fs, df, df_sym, d2f
   real(8) :: f, fprime, fsecond

!
!  Input and grid initialization
!
   write (*,'(a,$)') " imput value for x >> "
   read*, x
   fx = f(x) ; fp = fprime(x) ; fs = fsecond(x)
   write (*,*) " x      =", x
   write (*,*) " f(x)   =", fx
   write (*,*) " f'(x)  =", fp
   write (*,*) " f''(x) =", fs
1  continue
! 
!  Numerical derivatives starts here
!
   write (*,'(a,$)') " imput value for dx (stop <=0) >> "
   read*, dx
   if (dx <= 0.d0) stop
!
! asymmetric difference
!
   df = (f(x+dx) - fx)/dx
   open (unit=11, file="2asy.dat", status="unknown", position="append")
   write (11,*) "asymmetric", dx, df, abs(df-fp)
   close (11)
!
! symmetric difference
!
   df_sym = (f(x+dx) -f(x-dx))/2/dx

! ...
   open (unit=12, file="2sym.dat", status="unknown", position="append")
   write (12,*) "symmetric", dx, df_sym, abs(df_sym-fp)
   close(12)
!
! second derivative
!

! ...
   d2f = (f(x+dx) + f(x-dx) - 2*f(x))/dx**2

   open (unit=13, file="2sec.dat", status="unknown", position="append")
   write (13,*) " second", dx, d2f, abs(d2f-fs)
   close(13)
!
! go back for another increment
!
  go to 1

end program numerical_derivative

real(8) function f(x)
   real(8), intent(in) :: x
   f = sin(x)
end function

real(8) function fprime(x)
   real(8), intent(in) :: x
   fprime = cos(x)
end function

real(8) function fsecond(x)
   real(8), intent(in) :: x
   fsecond = -sin(x)
end function

