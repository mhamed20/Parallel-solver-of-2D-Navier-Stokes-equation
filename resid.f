c-----------------------------------------------------------------------
      subroutine RESID( ni, nj, u, v, p, res, x, y, eps, dt )

***********************************************************************
***
*** Resid computes a forward Euler step of the incompressible
*** Navier-Stokes equations with artificial compressibility.
***
*** A second order difference method with three point stencil is used.
***
*** Input: ni, nj  - Size of grid.
***        u, v, p - Velocity and pressure at t_n
***        res     - Work space of size 3*ni*nj
***        x, y    - The grid.
***        eps     - Viscosity parameter, epsilon in the equations.
***        dt      - Time step.
***
*** Output: u, v, p - The velocity and pressure is updated to t_{n+1} by
***                   forward Euler at the points (2..ni-1,2..nj-1)
***         res(1),res(2),res(3) - The first three elements of the
***         work space vector contains the maximum norm of the residual
***         for the three equations respectively.
***
*** NOTE: FORTRAN stores arrays column-wise, i.e., the first index
***       changes fastest. (ind = i + ni*(j-1) )
***
***********************************************************************

      implicit none


      integer ni, nj, i, j
      real*8  u(ni,nj), v(ni,nj), p(ni,nj), eps, dt, x(ni,nj), y(ni,nj)
      real*8  res(3,ni,nj), m11h, m12h, m21h, m22h, cx, cy, flx, jloc
      real*8  neps, mresu, mresv, mresp

      neps = 0.1*eps
      do j=1,nj
         do i=1,ni
            res(1,i,j) = 0
            res(2,i,j) = 0
            res(3,i,j) = 0
         enddo
      enddo
       

      do j=2,nj-1
         do i=2,ni

            m11h = 0.25d0*(y(i,j+1)-y(i,j-1)+y(i-1,j+1)-y(i-1,j-1))
            m12h =-0.25d0*(x(i,j+1)-x(i,j-1)+x(i-1,j+1)-x(i-1,j-1))
            m21h = -(y(i,j)-y(i-1,j))
            m22h =  (x(i,j)-x(i-1,j))

            jloc = m11h*m22h-m12h*m21h

            flx = (u(i,j)*u(i,j)+u(i-1,j)*u(i-1,j) +p(i,j)+p(i-1,j))
     *             *m11h*0.5d0 +
     *            (u(i,j)*v(i,j)+u(i-1,j)*v(i-1,j))*m12h*0.5d0

            cx = (m11h*(u(i,j)-u(i-1,j))
     *         +  m21h*0.25d0*(
     *            u(i,j+1)-u(i,j-1)+u(i-1,j+1)-u(i-1,j-1)))/jloc

            cy = (m12h*(u(i,j)-u(i-1,j))
     *         +   m22h*0.25d0*(
     *            u(i,j+1)-u(i,j-1)+u(i-1,j+1)-u(i-1,j-1)))/jloc

            flx = ( -flx + eps*( cx*m11h + cy*m12h ) )
c     *              + neps*(p(i,j)-p(i-1,j))
         
            res(1,i,j)   = res(1,i,j)   - flx
            res(1,i-1,j) = res(1,i-1,j) + flx

            flx = (u(i,j)*v(i,j)+u(i-1,j)*v(i-1,j))
     *             *m11h*0.5d0 +
     *            (v(i,j)*v(i,j)+v(i-1,j)*v(i-1,j)+
     *                           p(i,j)+p(i-1,j))*m12h*0.5d0

            cx = (m11h*(v(i,j)-v(i-1,j))
     *         +  m21h*0.25d0*(
     *            v(i,j+1)-v(i,j-1)+v(i-1,j+1)-v(i-1,j-1)))/jloc

            cy = (m12h*(v(i,j)-v(i-1,j))
     *         +   m22h*0.25d0*(
     *            v(i,j+1)-v(i,j-1)+v(i-1,j+1)-v(i-1,j-1)))/jloc

            flx = ( -flx + eps*( cx*m11h + cy*m12h ) )
c     *              + neps*(p(i,j)-p(i-1,j))

            res(2,i,j)   = res(2,i,j)   - flx
            res(2,i-1,j) = res(2,i-1,j) + flx

            flx = -( (u(i,j)+u(i-1,j))*m11h*0.5d0 +
     *               (v(i,j)+v(i-1,j))*m12h*0.5d0 )

            cx = (m11h*(p(i,j)-p(i-1,j))
     *         +  m21h*0.25d0*(
     *            p(i,j+1)-p(i,j-1)+p(i-1,j+1)-p(i-1,j-1)))/jloc

            cy = (m12h*(p(i,j)-p(i-1,j))
     *         +   m22h*0.25d0*(
     *            p(i,j+1)-p(i,j-1)+p(i-1,j+1)-p(i-1,j-1)))/jloc

            flx =  flx + neps*( cx*m11h + cy*m12h ) 

            res(3,i,j)   = res(3,i,j) - flx
            res(3,i-1,j) = res(3,i-1,j) + flx

         enddo
      enddo

      do j=2,nj
         do i=2,ni-1

            m21h =-0.25d0*(y(i+1,j)-y(i-1,j)+y(i+1,j-1)-y(i-1,j-1))
            m22h = 0.25d0*(x(i+1,j)-x(i-1,j)+x(i+1,j-1)-x(i-1,j-1))
            m11h =  y(i,j)-y(i,j-1)
            m12h =-(x(i,j)-x(i,j-1))
            jloc = m11h*m22h-m12h*m21h

            flx = (u(i,j)*u(i,j)+u(i,j-1)*u(i,j-1) + p(i,j)+p(i,j-1))*
     *             m21h*0.5d0+
     *            (u(i,j)*v(i,j)+u(i,j-1)*v(i,j-1))*m22h*0.5d0

            cx = (m21h*(u(i,j)-u(i,j-1))
     *         +  m11h*0.25d0*(
     *            u(i+1,j)-u(i-1,j)+u(i+1,j-1)-u(i-1,j-1)))/jloc

            cy = (m22h*(u(i,j)-u(i,j-1))
     *         +  m12h*0.25d0*(
     *            u(i+1,j)-u(i-1,j)+u(i+1,j-1)-u(i-1,j-1)))/jloc

            flx = ( -flx + eps*( cx*m21h + cy*m22h ) )
         
            res(1,i,j)   = res(1,i,j)   - flx
            res(1,i,j-1) = res(1,i,j-1) + flx

            flx = (u(i,j)*v(i,j)+u(i,j-1)*v(i,j-1))*m21h*0.5d0+
     *            (v(i,j)*v(i,j)+v(i,j-1)*v(i,j-1)+p(i,j)+p(i,j-1))
     *             *m22h*0.5d0

            cx = (m21h*(v(i,j)-v(i,j-1))
     *         +  m11h*0.25d0*(
     *            v(i+1,j)-v(i-1,j)+v(i+1,j-1)-v(i-1,j-1)))/jloc

            cy = (m22h*(v(i,j)-v(i,j-1))
     *         +  m12h*0.25d0*(
     *            v(i+1,j)-v(i-1,j)+v(i+1,j-1)-v(i-1,j-1)))/jloc

            flx = ( -flx + eps*( cx*m21h + cy*m22h ) )           

            res(2,i,j)   = res(2,i,j)   - flx
            res(2,i,j-1) = res(2,i,j-1) + flx

            flx = -( (u(i,j)+u(i,j-1))*m21h*0.5d0 +
     *               (v(i,j)+v(i,j-1))*m22h*0.5d0 )

            cx = (m21h*(p(i,j)-p(i,j-1))
     *         +  m11h*0.25d0*(
     *            p(i+1,j)-p(i-1,j)+p(i+1,j-1)-p(i-1,j-1)))/jloc

            cy = (m22h*(p(i,j)-p(i,j-1))
     *         +  m12h*0.25d0*(
     *            p(i+1,j)-p(i-1,j)+p(i+1,j-1)-p(i-1,j-1)))/jloc

            flx = ( flx + neps*( cx*m21h + cy*m22h ) )

            res(3,i,j)   = res(3,i,j)   - flx
            res(3,i,j-1) = res(3,i,j-1) + flx

         enddo
      enddo

      mresu = -1
      mresv = -1
      mresp = -1

      do j=2,nj-1
         do i=2,ni-1 
            jloc = 0.25d0*( (x(i+1,j)-x(i-1,j))*(y(i,j+1)-y(i,j-1))-
     *                      (x(i,j+1)-x(i,j-1))*(y(i+1,j)-y(i-1,j)))       
            u(i,j) = u(i,j) + dt/jloc*res(1,i,j)
            v(i,j) = v(i,j) + dt/jloc*res(2,i,j)
            p(i,j) = p(i,j) + dt/jloc*res(3,i,j)
            mresu = MAX( mresu, ABS(res(1,i,j)/jloc ))
            mresv = MAX( mresv, ABS(res(2,i,j)/jloc ))
            mresp = MAX( mresp, ABS(res(3,i,j)/jloc ))
         enddo
      enddo
      res(1,1,1) = mresu
      res(2,1,1) = mresv
      res(3,1,1) = mresp

      return
      end

