! Original script David Manolopoulos (Oxford)

      program main
      use constants
      use math 
      implicit none 
      
      integer(kind=4),parameter :: n=255,iverbose=0,nmode =1
      real(kind=8)    :: h(n,n),g(n,n),g2(n,n),e(n),f(n),w(n),q(n)

      real(kind=8)    :: qmin,qmax,mass,mass_aux
      real(kind=8)    :: on2m,pot
      integer(kind=8) :: i,j,k

      real(kind=8)    :: tempK,beta,tmax,z,c,dt,t,aux
      integer(kind=8) :: it,nt,ierr 

      real(kind=8)    :: xi ,xf, dx
      integer(kind=8) :: nsize
      real(kind=8), dimension(NDIM) :: V

!     INIT
      h(:,:) = 0d0
      g(:,:) = 0d0
      g2(:,:) = 0d0
      e(:)   = 0d0
      f(:)   = 0d0
      w(:)   = 0d0
      q(:)   = 0d0

  
!     Parameters - take care, these are system dependent

      mass_aux = 2.0 *1837.289d0 ! mass H
      mass = 1.0 !MS 
!     em=0.9482044362022515*1822.8885 ! for OH
 
      qmin = -3d0
      qmax = 3d0

      !scale by mass
      qmin = qmin*dsqrt(mass_aux) 
      qmax = qmax*dsqrt(mass_aux) 


      !acf
      tempK=400.0
      beta = 1 / (tempK * 3.1668152e-06)
      tmax = 100000.0 
      nt   = NINT(tmax/10.)

!     Read grid

      if (nmode .EQ. 2) THEN 
      call LECTOR(xi,xf,dx,nsize,V)

      if (iverbose .GE. 1) THEN
      write(6,*) ''
      write(6,*) 'Grid loaded successfully'
      WRITE(6,*) 'nsize',nsize
      WRITE(6,*) 'xi',xi*au2A,'xf',xf*au2A,'dx',dx*au2A
      write(6,*) ''
      ENDIF

      ENDIF

!-------------------------------------------------------------------------------------------
 
!     Calculation
 
      on2m = 0.5d0*hbar*hbar/mass
      call sindvr (on2m,qmin,qmax,n,h,w,q)


      if (iverbose .GE. 2) THEN
        do i=1,n
        do j=1,n
        write(6,*) h(i,j)
        enddo
        enddo
     
        write(6,*)
        do i=1,n
        write(6,*) q(i)
        enddo

        write(6,*)
        do i=1,n
        write(6,*) w(i)
        enddo
      ENDIF
   

      if (nmode .EQ. 1) THEN 
       !quartic
       do j = 1,n
          call calcquartic(q(j),pot,mass_aux)
          h(j,j) = h(j,j)+pot
       enddo
      
      ELSE IF (nmode .EQ. 2) THEN  
      !GRID 
      do j = 1,n
         if ( (q(j) .LT. xi ) .OR. ( q(j) .GT. xf ) ) THEN
           write(6,*) q(j),xi,xf
           stop 'EXTRAPOLATION'
         endif   
         i = NINT( (q(j)-xi)/dx  )
         pot = V(i)
         h(j,j) = h(j,j)+pot
      enddo
      ENDIF

      !------write H before diag---------
      IF (iverbose .GE. 2) THEN
         write(6,*)
         do i=1,n
         do j=1,n
         write(6,*) h(i,j)
         enddo
         enddo
      ENDIF

     !----diagonalize------------ 
!      call symevp (h,n,n,e,ierr)
!      if (ierr .ne. 0) stop 'main 1'
       call diag(h,e,n,iverbose)

      write(6,*)
      do j = 1,n
         if (j .GT. 1) THEN
            write(6,100) j,e(j)*au2eV,(e(j)-e(j-1))*au2eV,(e(j)-e(j-1))/invcm2au 
         ELSE
            write(6,101) j,e(j)*au2eV,' (ZPE:',e(j)/invcm2au,')'
         endif
      enddo

      OPEN (5,FILE='pot.dat')

      
      do j = 1,n
         if (nmode .EQ. 1) THEN
             call calcquartic(q(j),pot,mass_aux)
         elseif (nmode .EQ. 2) THEN
             if ( (q(j) .LT. xi ) .OR. ( q(j) .GT. xf ) ) THEN
                write(6,*) q(j),xi,xf
                stop 'EXTRAPOLATION'
             endif
             i = NINT( (q(j)-xi)/dx  )
             pot = V(i)
         endif
         write(5,*) q(j)*au2A/dsqrt(mass_aux),pot*au2eV
      enddo
      CLOSE(5)

      OPEN (5,FILE='levels.dat')
      do j = 1,10
         write(5,*)
         write(5,*) qmin*au2A/dsqrt(mass_aux), e(j)*au2eV
         write(5,*) qmax*au2A/dsqrt(mass_aux), e(j)*au2eV
      enddo
      CLOSE(5)
      

100   FORMAT (I10,2E20.10,F10.2)
101   FORMAT (I10,E20.10,A6,F10.2,A2)


     !----compute acf david-------------------------------

   
      OPEN (5,FILE='acf.dat')
      z = 0.d0
      do j = 1,n
    
         f(j) = dexp(-beta*e(j))
         z = z+f(j)
         do i = 1,j
            g(i,j) = 0.d0
            do k = 1,n
               g(i,j) = g(i,j)+q(k)*h(k,i)*h(k,j)
            enddo
            if (i .lt. j) then
               !g2(i,j)  = (f(j)-f(i))*(e(i)-e(j))*g(i,j)*g(i,j)
               g2(i,j)  = (f(j)-f(i))*((e(i)-e(j)))*g(i,j)*g(i,j)
            else
               g2(j,j) = 0.d0
            endif
         enddo
      enddo

      dt = tmax/nt
      dt = dt/hbar

      do it = 0,nt
         t = it*dt
         c = 0.d0
         do j = 1,n
            do i = 1,j-1
               c = c+g2(i,j)*dcos((e(i)-e(j))*t)
            enddo
         enddo
         c = 2.d0*c/(beta*z*hbar*hbar)
         write(5,*) t,c
      enddo
     
      close(5)
!----------------coupling elements ------------------------------
      OPEN (5,FILE='coupling.dat')

      do j =1,10
      do i =1,j-1
         write(5,*) j,i,g(i,j),g2(i,j),g2(i,j)/g2(2,3)
      enddo
      enddo 
      close(5)


      endprogram

!-------------------------------------------------------------------------------------------
      !call sindvr (on2m,qmin,qmax,n,h,w,q)
      subroutine sindvr (on2m,a,b,n,t,w,x)
      IMPLICIT NONE
      integer(kind=4) ::n
      integer(kind=8) :: i,j,k,m
      real(kind=8) :: t(n,n),w(n),x(n)
      real(kind=8) :: s(2*n)
      real(kind=8) :: pi,alfa,beta,cosa,sina,cosb,sinb,temp
      real(kind=8) :: on2m,a,b,dx,wt,t0

      !implicit double precision (a-h,o-z)
 
!     ----------------------------------------------------------------- 
!     n point sine DVR in a < x < b
!     D.T.Colbert and W.H.Miller, JCP 96 (1992) 1982, eq.(A6)
!     ----------------------------------------------------------------- 
 
 
      m = n+1
      pi = acos(-1.0d0)
      alfa = 0.5d0*on2m*(pi/(b-a))**2
      beta = pi/(2*m)
      cosa = 1.0d0
      sina = 0.0d0
      cosb = cos(beta)
      sinb = sin(beta)
      dx = (b-a)/m
      wt = sqrt(dx)
      t0 = alfa*(2*m**2+1)/3.0d0

      !ALBERTO
      !write(6,*) m
      !write(6,*) pi
      !write(6,*) alfa
      !write(6,*) beta
      !write(6,*) cosa
      !write(6,*) sina
      !write(6,*) cosb
      !write(6,*) sinb
      !write(6,*) dx
      !write(6,*) wt
      !write(6,*) t0
      !write(6,*)

      do k = 1,2*n
         alfa = -alfa
         temp = cosa*sinb+sina*cosb
         cosa = cosa*cosb-sina*sinb
         sina = temp
         s(k) = alfa/sina**2
      enddo
      do j = 2,n
         do i = 1,j-1
            t(i,j) = s(j-i)-s(j+i)
            t(j,i) = t(i,j)
         enddo
      enddo
      do j = 1,n
         t(j,j) = t0-s(j+j)
         w(j) = wt
         x(j) = a+j*dx
      enddo
      return
      end 
!-------------------------------------------------------------------------------------------
      subroutine symevp (a,lda,n,d,ierr)
      implicit double precision (a-h,o-z)
      integer(kind=8) :: ierr
!
!     ----------------------------------------------------------------- 
!     This subroutine uses LAPACK DSYEV to
!     diagonalise a real symmetric matrix.
!     ----------------------------------------------------------------- 
!
      dimension a(lda,n),d(n)
      dimension work(34*n)
 
      lwork = 34*n
      call dsyev ('V','U',n,a,lda,d,work,lwork,ierr)
      return
      end

      subroutine calcquartic (r, en,mass)
      use constants
      implicit none
      double precision, intent(in) :: r
      double precision  :: a0,a2,a4
      double precision  :: mass
      double precision, intent(out) :: en

!     ----------------------------------------------------------------- 
!     Returns quartic potential
!     ----------------------------------------------------------------- 

      !Coefficients in eV, eV/A and eV/A^2 
      a0 = 0.51515678 
      a2 = -9.93809966
      a4 = 48.64470511

      !------------------------------------------

      a0 = a0  * eV2au
      a2 = a2  * eV2au * (au2A**2) 
      a4 = a4  * eV2au * (au2A**4)

      !scale by mass
      a2  = a2  /mass
      a4  = a4  /mass/mass

      en = a0 + a2 * (r**2) + a4 * (r**4) 
      en = en 
      return
      end
