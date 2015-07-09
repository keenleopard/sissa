program main
   call numerical_integration
end program

subroutine numerical_integration
!----------------------------------------------------------------------
   implicit none

   real(8), parameter :: pi = 3.14159265358979323846d0

   integer :: n,i
   real(8) :: rectangular,trapezoidal,simpson
   real(8) :: x, dx, fsave, f
   external f

1  continue    ! this is where the code will came back for the next input

!
!  read number of sub-intervals
!
   write (*,'(a,$)') "number of sub-intervals (stop if <=0) >> "
   read*, n
   if (n<=0) stop
   dx = 0.5*pi/n
   write (*,'(a,i10,a,1pe15.8)') " n = ",n, ", dx =", dx

!
!  perform integration by the rectangular rule here
!
   rectangular = 0.d0
   x=0.d0
   do i = 0, n-1
      x = i * dx
      rectangular = rectangular + f(x)
   end do
   rectangular = rectangular * dx

   open (unit=1, file="rect.dat", status="unknown", position = "append")
   write (1,*) "rectangular:", n, rectangular, rectangular - 1.d0
   close (1)

!
!  perform integration by the trapezoidal rule here
!
   trapezoidal = 0.d0
   x=0.d0
!   do i = 1, n-1
!      x = i * dx
!      trapezoidal = trapezoidal + f(x)
!   end do 
!   trapezoidal = (trapezoidal + f(1.d0)/2) * dx
   do i = 0, n -1
      x = i * dx
      trapezoidal = trapezoidal + (f(x)+f(x+dx))/2.d0
   end do
   trapezoidal = trapezoidal * dx

   open (unit=2, file="trap.dat", status="unknown", position="append")
   write (2,*) "trapezoidal:", n,trapezoidal, abs(trapezoidal - 1.d0)
   close (2)

!
!  perform integration by the simpson rule here 
!  BEWARE : n must be even here  [i.e. n+1 must be odd]
!
   if ( mod(n,2) == 1) then
      write (*,*) "Simpson integration requires an even number of sub-intervals"
      go to 1
   end if

!   simpson = f(0.d0)+f(1.d0)
!   x=0.d0
!   do i = 1, n-1
!      x = i * dx
!      if (mod(i,2) == 1) then
!          simpson = simpson + 4 * f(x)
!      else 
!          simpson = simpson + 2 * f(x)
!      end if 
!   end do 
!   simpson = simpson * dx/3.d0
   simpson = 0.d0
   x = 0.d0
   do i = 0, n-2, 2
      x = i * dx
      simpson = simpson + (f(x) + f(x+dx)*4.d0 + f(x + 2*dx)) /3 
   end do 
   simpson = simpson *dx
   print *, simpson

   open (unit=3, file="simp.dat", status="unknown", position="append")
   write (3,*) "simpson:",n, simpson,     abs(simpson     - 1.d0)
   close (3)

!
!  go back and read another number of sub-intervals
!
   go to 1

end subroutine numerical_integration
!------------------------------------------------------------------------------
real(8) function f(x)
!
   implicit none
   real(8), intent (in) :: x
   f = cos(x)
end function
!------------------------------------------------------------------------------
