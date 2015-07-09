program harmonic
!------------------------------------------------------------------------------
!
!  Solve the quantum harmonic oscillator by forward integration
!
!    H phi = e phi   :     -1/2 phi''  + (1/2 x**2 -e) phi = 0
!
!  The symmetry of the potential [ V(-x) = V(x) ] is exploited by integrating
!  the differential equation  only in the [0,xmax] interval and reconstructing 
!  the solution in the [-xmax,0] interval by symmetry.
!
!  Eigenvalue search using "shooting method" or by bisection based on 
!  the required number of nodes (crossing of the horizontal axis).
!
   implicit none

!  Maximum number of points in the mesh
   integer, parameter :: mshx = 2000
   real(8), parameter :: pi = 3.14159265358979323846d0, eps12 = 1.d-12

   real(8) :: x(0:mshx), y(0:mshx), p(0:mshx), vpot(0:mshx), f(0:mshx)
   real(8) :: xmax, dx, ddx, xmcl, norm
   integer :: nodes, hnodes, ncross, iteration
   integer :: mesh, i, icl
   real(8) :: e, eup, elw
   logical :: iterate
   character*80 :: fileout
!
!  Input and grid initialization
!
   print '(a,$)','Maximum x (typical value: 10) ? '
   read*, xmax
   print '(a,i5,a,$)', 'Number of mesh points (max:',mshx,') ? '
   read*,mesh
   if (mesh > mshx) stop
   dx =  xmax/mesh 
   ddx = dx*dx
!
!  Potential definition (symmetric w.r.t x=0)
!
   do i = 0, mesh
      x(i) = dfloat(i) * dx
      vpot(i) = 0.5d0 * x(i)*x(i)
   enddo

   print '(a,$)','Output file name = '
   read '(a)',fileout
   open(1,file=fileout,form="formatted",status="new")

!
!  Entry for eigenvalue search loop 
!
999  continue
!
!  Input number of nodes (<0 to terminate)
!
   print '(a,$)', 'Number of nodes (-1 = stop) ? '
   read*,nodes
   if (nodes < 0) then
      close(1)
      stop 
   endif
!
!  Define eigenvalue bounding limits
!
   eup = vpot(mesh)
   elw = eup
   do i=0,mesh
      elw = min(elw,vpot(i))
      eup = max(eup,vpot(i))
   end do
!
!  Define energy for single shot evaluation
!
   print '(a,$)','Trial energy (0 = search by bisection) ? '
   read*,e
   if ( e ==  0.d0 ) then
!  Automatic eigenvalue search mode
      e = 0.5d0 * (elw + eup)
      iterate = .true.
   else
!  Single trial energy mode
      iterate = .false.
   endif

   iteration = 0
!
!  Entry for solution at a given energy e
!
1  continue
   iteration = iteration + 1
!
! define the 
!
   f(0) = ddx * 2.d0 * (vpot(0)-e)
   do i=1,mesh
!
!     f <= 0  classicaly allowed region
!     f > 0   classicaly forbidden region
!
      f(i) = ddx * 2.d0 * (vpot(i)-e)
!
!     no f(i) can vanish, otherwise the change of sign 
!     is not detected. Next line is a trick to avoid that.
!
      if( abs(f(i)) < eps12 ) f(i) = eps12    ! *** From Peng: This is very clever!
!
!     remember the index of the last sign change.
!     in the assumption that the potential is growuimg for large x this maks the!     limit of the classically allowed region
!
      if( f(i) .ne. sign(f(i),f(i-1)) ) icl=i
   end do

   if (icl >= mesh-2) then
      print*,'Last sign change is too far.'
      stop
   endif
!
   do i=0,mesh
      y(i) = 0.d0
   enddo
!
!  definition of the wfc in the first two points.
!
   hnodes = nodes/2
   if (2*hnodes.eq.nodes) then
      y(0) = 1.d0
      y(1) = y(0) + 0.5d0 * f(0)*y(0)
   else
      y(0) = 0.d0 
      y(1) = dx
   endif
!
!  outward integration and node counting
!
   ncross=0
   do i=1,mesh-1
      y(i+1) = (2.d0+f(i))*y(i)-y(i-1)
      if ( y(i) .ne. sign(y(i),y(i+1)) ) ncross=ncross+1
   enddo
!
!  verify the number of nodes
!
   print '(2i4,f14.8)',iteration, ncross, e
   if ( iterate ) then
      if (ncross.gt.hnodes) then
!        Too many sign chages =>  The energy is too high.
         eup = e
      else 
!        Right number of or too few sign changes => The energy is too low.
         elw = e
      endif
!     new trial energy value
      e = 0.5d0 * (eup+elw)
!     convergence criterium
      if (eup-elw .gt. eps12) go to 1
   endif
!
!  convergence has been acheived, or iteration was not required.
!  *** BEWARE ***  y is NOT normalized yet

!
!  Renormalization so that
!  int |psi|^2 dx = 1:   *** WE ARE NOT GOING TO DO THAT NOW ***
!  The problem is the divergence at large x values due to the fact that 
!  the outword integration is numerically unstable when V > e .

!  NORMALIZATION:
!
!  The following would be an integration by trapezoidal rule... 
!  ***assuming*** that y(xmax+dx)=0   AND  a symmetric interval [-xmax,xmax] 
!  It should be fine AFTER the problem of divergence has been solved
!
!--renormalization starts here
!
!  norm = 0.5d0 * y(0)*y(0)
!  do i=1,mesh
!     norm = norm + y(i)*y(i)
!  enddo
!  norm = dx * norm 
!
!  norm = 2.d0 * norm ! THE FACTOR OF 2 BECAUSE WE INTEGRATE IN [-xmax,xmax] 
!
!  norm = sqrt(norm)
!  do i=0,mesh
!     y(i) = y(i) / norm
!  enddo

!--renormalization ends here
!

!
!  classical probability density for energy e
!  (for comparison):
!
   xmcl = sqrt(2.d0*e)
   do i=icl,mesh
      p(i) = 0.d0
   enddo
   do i=0,icl-1
      p(i) = 1.d0/sqrt(xmcl**2 - x(i)**2)/pi
   enddo

!  WRITE SOLUTION ON FILE
   write (1,'(a)') "# x(i)     phi(i)          |phi(i)|**2   p_classical(i)    vpot(i)"
!
!  x <  0    ! this part is generated using the known symmetry for x->-x
!
   do i=mesh,1,-1
      write (1,'(f7.3,3e16.8,f12.6)') -x(i),(-1)**nodes*y(i),y(i)*y(i),p(i),vpot(i)
   enddo
!
!  x > 0     ! this is the part we actually calculated
!
   do i=0,mesh
      write (1,'(f7.3,3e16.8,f12.6)') x(i),y(i),y(i)*y(i),p(i),vpot(i)
   enddo

!
! go back for another eigenvalue
!
   go to 999

end program harmonic
