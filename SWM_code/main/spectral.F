      subroutine spinit
 
c  this routine must be called before using any of the other routines.
c  (especially before "ginit" which is part of diag.F)
 
#include "param.h"
#include "spec.h"
      common/nc/ work(1000)
      print *,'Calling spinit. Must be called before ginit.'
 
      a = 6.371e6
      call parmtr
      call cffti(ix,work)
      print*,'in spinit a=',a  
      return
      end
*******************************************************************
      subroutine gaussl(x,w,m)

      double precision z,z1,p1,p2,p3,pp,eps
      dimension x(m),w(m)
      parameter (eps=3.d-14)

c    a slightly modified version of a program in numerical recipes
c
c   input:
c      m    = number of gaussian latitudes between pole and equator
c   output:
c      x(m) = sin(gaussian latitudes) (from north pole to equator)
c      w(m) = weights in gaussian quadrature -- cos(theta)d(theta)
c             (sum should equal 1.0)

      n=2*m

      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))

    1   continue

        p1=1.d0
        p2=0.d0

        do 11 j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
   11   continue

        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp

        if(abs(z-z1).gt.eps) go to 1

        x(i)=z
        w(i)=2.d0/((1.d0-z*z)*pp*pp)

   12 continue

      return
      end

c*******************************************************************
      subroutine ncwg2(xm1,xm2,x1,x2)

#include "param.h"

      dimension x1(ix),x2(ix)
      dimension z(2,ix),xm1(2,nwaves),xm2(2,nwaves)
      common/nc/ work(1000)

c   inverse of ncgw2

      xm1(2,1) = 0.0
      xm2(2,1) = 0.0
      do 1 m=1,nwaves
        z(1,m)=xm1(1,m)-xm2(2,m)
        z(2,m)=xm1(2,m)+xm2(1,m)
    1 continue

      mm1=nwaves+1
      n2=ix+2

      do 2 m=mm1,ix
        z(1,m)=xm1(1,n2-m)+xm2(2,n2-m)
        z(2,m)=-xm1(2,n2-m)+xm2(1,n2-m)
    2 continue

      call cfftb(ix,z,work)

      do 3 i=1,ix
        x1(i)=z(1,i)
        x2(i)=z(2,i)
    3 continue

      return
      end
***************************************************************
      subroutine ncgw2(x1,x2,xm1,xm2)

#include "param.h"

      dimension x1(ix),x2(ix)
      dimension z(2,ix),xm1(2,nwaves),xm2(2,nwaves)
      common/nc/work(1000)

c   computes fourier transform of two real series, x1 and x2

      anorm=1./float(ix)
      bnorm=0.5/float(ix)

      do 1 i=1,ix
        z(1,i)=x1(i)
        z(2,i)=x2(i)
    1 continue

      call cfftf(ix,z,work)

      n2=ix+2
      xm1(1,1)=z(1,1)*anorm
      xm1(2,1)=0.
      xm2(1,1)=z(2,1)*anorm
      xm2(2,1)=0.

      do 2 m=2,nwaves
        xm1(1,m)=bnorm*(z(1,m)+z(1,n2-m))
        xm1(2,m)=bnorm*(z(2,m)-z(2,n2-m))
        xm2(1,m)=bnorm*(z(2,m)+z(2,n2-m))
        xm2(2,m)=bnorm*(-z(1,m)+z(1,n2-m))
    2 continue

      return
      end
******************************************************************
      subroutine parmtr

#include "param.h"
#include "spec.h"

      dimension repsi(mxp,jxp),consq(mxp)
      dimension paly(mxp,jx),alp(mxp,jx)
      dimension wt1(il)
      dimension emm(mxp),ell(mxp,jxp),el2p(mx,jxp),epsi(mxp,jxp)

c  computes quantities needed by other spectral routines
c  (radius of earth, a, must be predefined in common/gauss)

      call gaussl(sia,wt,iy)

      am1 = 1./a
      am2=1./(a*a)

      do 1 l=1,iy
        cosqr = 1.0-sia(l)**2
        coa(l)=sqrt(cosqr)
        coam(l)=1./coa(l)
        wt1(l)=wt(l)/cosqr
    1 continue

      do 2 m=1,mx
      do 2 n=1,jx
        mm(m)=isc*(m-1)
        ll(m,n)=mm(m)+n-1
        l2(m,n)=ll(m,n)*(ll(m,n)+1)
        el2(m,n)=float(l2(m,n))*am2
        el4(m,n)=el2(m,n)*el2(m,n)
    2 continue

      elm2(1,1)=0.
      do 3 m=2,mx
      do 3 n=1,jx
        elm2(m,n)=1./el2(m,n)
    3 continue
      do 4 n=2,jx
        elm2(1,n)=1./el2(1,n)
    4 continue

c  for internal use in this subroutine, it is useful to have el2=l*(l+1)
c    to jxp=jx+1

      do 5 m=1,mx
      do 5 n=1,jxp
        llp=mm(m)+n-1
        l2p=llp*(llp+1)
        el2p(m,n)=float(l2p)*am2
    5 continue

c    in a sector model retaining m= 0,3,6,9,12,15, then mx = 6, mxp = 16
c    all m values up to mxp are needed by recursion relation for legendr
c    epsi is a useful quantity when differentiating or generating the po

      do 10 m=1,mxp
      do 10 n=1,jxp
        emm(m)=float(m-1)
        ell(m,n)=float(n+m-2)
        emm2=emm(m)**2
        ell2=ell(m,n)**2
        if(n.eq.jxp) then
          epsi(m,n)=0.0
        else if(n.eq.1.and.m.eq.1) then
          epsi(m,n)=0.0
        else
          epsi(m,n)=sqrt((ell2-emm2)/(4.*ell2-1.))
        endif
        repsi(m,n)=0.0
        if(epsi(m,n).gt.0.) repsi(m,n)=1./epsi(m,n)
   10 continue

c  compute quantities required by subroutines grad and uvspec

      do 15 m=1,mx
      do 15 n=1,jx
        m1=mm(m)
        m2=m1+1
        el1=float(ll(m,n))
        if(n.eq.1) then
          gradx(m)=float(m1)/a
          uvdx(m,1)=-a/float(m1+1)
          uvdym(m,1)=0.0
        else
          uvdx(m,n)=-a*float(m1)/(el1*(el1+1))
          gradym(m,n)=(el1-1.)*epsi(m2,n)/a
          uvdym(m,n)=-a*epsi(m2,n)/el1
        endif
        gradyp(m,n)=(el1+2.)*epsi(m2,n+1)/a
        uvdyp(m,n)=-a*epsi(m2,n+1)/(el1+1.)
   15 continue

c  generate associated legendre polynomials, normalized so that the inte
c   over the sphere equals unity

      sqrhlf=sqrt(.5)
      do 20 m=2,mxp
        consq(m) = sqrt(.5*(2.*emm(m)+1.)/emm(m))
   20 continue

      do 30 j=1,iy

        y = coa(j)
        x = sia(j)

c  define (0,0) polynomial

        alp(1,1) = sqrhlf

c  use recursive formula for n=1 (m=l) diagonal

        do 31 m=2,mxp
          alp(m,1) = consq(m)*y*alp(m-1,1)
   31   continue

c  derive n=2 polynomials from n=1

        do 32 m=1,mxp
          alp(m,2)=(x*alp(m,1))*repsi(m,2)
   32   continue

c  use recursive formula to obtain n=3 to n=jx

        do 33 n=3,jx
        do 33 m=1,mxp
          alp(m,n)=(x*alp(m,n-1)-epsi(m,n-1)*alp(m,n-2))*repsi(m,n)
   33   continue

c  zero polynomials with absolute values smaller than 10e-30

        do 34 n=1,jx
        do 34 m=1,mxp
          if(abs(alp(m,n)) .ge. 1.0e-30) then
            paly(m,n)=alp(m,n)
          else
            paly(m,n)=0.0
          endif
   34   continue

c   pick off the required plynomials

        do 35 n=1,jx
        do 35 m=1,mx
          mm2=mm(m)+1
          pol(m,n,j)=paly(mm2,n)
   35   continue

c  differentiate, adding weights

        wght=am1*wt1(j)

        do 40 m=1,mx
          mm2=mm(m)+1
          polydx(m,1,j)= emm(mm2)*pol(m,1,j)*wght
          polydy(m,1,j)=-emm(mm2)*pol(m,2,j)*wght*
     *     epsi(mm2,2)
          polydx(m,jx,j)=0.0
          polydy(m,jx,j)=0.0
   40   continue

        do 41 n=2,jtrun1
        do 41 m=1,mx
          mm2=mm(m)+1
          polydx(m,n,j)=emm(mm2)*pol(m,n,j)*wght
          polydy(m,n,j)=
     *    ((ell(mm2,n)+1.)*epsi(mm2,n  )*pol(m,n-1,j)
     *    - ell(mm2,n)    *epsi(mm2,n+1)*pol(m,n+1,j))
     *    *wght
   41   continue

   30 continue

      return
      end
****************************************************************
      subroutine lap(strm,vorm)

#include "param.h"
#include "spec.h"

      complex vorm(mx,jx),strm(mx,jx)

c  computes the spectral laplacian, vorm, of the spectral field strm

      do 1 m=1,mx
      do 1 n=1,jx
        vorm(m,n)=-strm(m,n)*el2(m,n)
    1 continue

      return
      end
****************************************************************
      subroutine lap2(strm,vorm)

#include "param.h"
#include "spec.h"

      complex vorm(mx,jx),strm(mx,jx)

c  computes spectral laplacian**2

      do 1 m=1,mx
      do 1 n=1,jx
        vorm(m,n)=strm(m,n)*el4(m,n)
    1 continue

      return
      end
*******************************************************************
      subroutine invlap(vorm,strm)

#include "param.h"
#include "spec.h"

      complex vorm(mx,jx),strm(mx,jx)

c   computes spectral inverse laplacian, strm, of spectral field vorm

      do 1 m=1,mx
      do 1 n=1,jx
        strm(m,n)=-vorm(m,n)*elm2(m,n)
    1 continue

      return
      end
*********************************************************************
      subroutine grad(psi,psdx,psdy)

#include "param.h"
#include "spec.h"

      dimension psi(2,mx,jx),psdx(2,mx,jx),psdy(2,mx,jx)

c  in spectral domain, computes x and y spectral components,
c   multiplied by cos(latitude),
c   of the gradient of the spectral field psi

      do 1 n=1,jx
      do 1 m=1,mx
        psdx(2,m,n)=gradx(m)*psi(1,m,n)
        psdx(1,m,n)=-gradx(m)*psi(2,m,n)
    1 continue

      do 2 k=1,2
      do 2 m=1,mx
        psdy(k,m,1)=gradyp(m,1)*psi(k,m,2)
        psdy(k,m,jx)=-gradym(m,jx)*psi(k,m,jtrun1)
    2 continue

      do 3 k=1,2
      do 3 n=2,jtrun1
      do 3 m=1,mx
        psdy(k,m,n)=-gradym(m,n)*psi(k,m,n-1)+gradyp(m,n)*psi(k,m,n+1)
    3 continue

      return
      end
******************************************************************
      subroutine uvspec(vorm,divm,ucosm,vcosm)

#include "param.h"
#include "spec.h"

      dimension vorm(2,mx,jx),divm(2,mx,jx),ucosm(2,mx,jx),
     * vcosm(2,mx,jx),zc(2,mx,jx),zp(2,mx,jx)

c  computes spectral u and v (multiplied by cos(latitude)) from spectral
c  vorticity vorm and divergence divm

      do 1 n=1,jx
      do 1 m=1,mx
        zp(2,m,n)= uvdx(m,n)*vorm(1,m,n)
        zp(1,m,n)=-uvdx(m,n)*vorm(2,m,n)
        zc(2,m,n)= uvdx(m,n)*divm(1,m,n)
        zc(1,m,n)=-uvdx(m,n)*divm(2,m,n)
    1 continue

      do 3 k=1,2
      do 3 m=1,mx
        ucosm(k,m,1)=zc(k,m,1)-uvdyp(m,1)*vorm(k,m,2)
        ucosm(k,m,jx)=uvdym(m,jx)*vorm(k,m,jtrun1)
        vcosm(k,m,1)=zp(k,m,1)+uvdyp(m,1)*divm(k,m,2)
        vcosm(k,m,jx)=-uvdym(m,jx)*divm(k,m,jtrun1)
    3 continue

      do 4 k=1,2
      do 4 n=2,jtrun1
      do 4 m=1,mx
        vcosm(k,m,n)=-uvdym(m,n)*divm(k,m,n-1)+uvdyp(m,n)*
     *   divm(k,m,n+1)+zp(k,m,n)
        ucosm(k,m,n)= uvdym(m,n)*vorm(k,m,n-1)-uvdyp(m,n)*
     *   vorm(k,m,n+1)+zc(k,m,n)
    4 continue

      return
      end
*********************************************************************
      subroutine gridd(vorm,vorg,kcos)

#include "param.h"
#include "spec.h"

      dimension vorg(ix,il),var(ix,2)
      dimension vorm(2,mx,jx),polvm(2,mx,jx)
      dimension fmvar(2,nwaves,2),fmvar1(2,mx,2)

c  trasforms spectral field vorm to grid field vorg
c if  kcos =2, vorg is divided by cos(lat)

      do 1 k=1,2
      do 1 m=1,nwaves
        fmvar(k,m,1)=0.
        fmvar(k,m,2)=0.
    1 continue

      do 10 j=1,iy

        do 20 k=1,2
        do 20 n=1,jx
        do 20 m=1,mx
          polvm(k,m,n)=vorm(k,m,n)*pol(m,n,j)
   20   continue

        do 21 k=1,2
        do 21 m=1,mx
          fmvar1(k,m,1)=0.
          fmvar1(k,m,2)=0.
   21   continue

        do 22 k=1,2
        do 22 n=1,jx,2
          do 23 m=1,mx
            fmvar1(k,m,1)=fmvar1(k,m,1)+polvm(k,m,n)
   23     continue
   22   continue

        do 24 k=1,2
        do 24 n=2,jx,2
          do 25 m=1,mx
            fmvar1(k,m,2)=fmvar1(k,m,2)+polvm(k,m,n)
   25     continue
   24   continue

        do 26 k=1,2
        do 26 m=1,mx
          fmvar(k,m,1)=fmvar1(k,m,1)+fmvar1(k,m,2)
          fmvar(k,m,2)=fmvar1(k,m,1)-fmvar1(k,m,2)
   26   continue

        call ncwg2(fmvar(1,1,1),fmvar(1,1,2),var(1,1),var(1,2))

        j1=il+1-j
        if(kcos.eq.1) then
          do 27 i=1,ix
            vorg(i,j)=var(i,2)
            vorg(i,j1)=var(i,1)
   27     continue
        else
          do 28 i=1,ix
            vorg(i,j)=var(i,2)*coam(j)
            vorg(i,j1)=var(i,1)*coam(j)
   28     continue
        endif

   10 continue

      return
      end
*******************************************************************
      subroutine gridx(varm,vorg,kcos)

#include "param.h"
#include "spec.h"

      dimension vorg(ix,il),var(ix,2)
      dimension varm(2,mx,il),fmvar(2,nwaves,2)

c  transforms varm (spectral in x, grid in y) to grid field vorg
c  if kcos =2, vorg is divided by cos(lat)

      do 1 k=1,2
      do 1 m=1,nwaves
        fmvar(k,m,1)=0.
        fmvar(k,m,2)=0.
    1 continue

      do 10 j=1,iy

        j1=il+1-j

        do 20 k=1,2
        do 20 m=1,mx
          fmvar(k,m,1)=varm(k,m,j)
          fmvar(k,m,2)=varm(k,m,j1)
   20   continue

        call ncwg2(fmvar(1,1,1),fmvar(1,1,2),var(1,1),var(1,2))

        if(kcos.eq.1) then
          do 21 i=1,ix
            vorg(i,j)=var(i,1)
            vorg(i,j1)=var(i,2)
   21     continue
        else
          do 22 i=1,ix
            vorg(i,j)=var(i,1)*coam(j)
            vorg(i,j1)=var(i,2)*coam(j)
   22     continue
        endif

   10 continue

      return
      end
*******************************************************************
      subroutine gridy(vorm,varm)

#include "param.h"
#include "spec.h"

      dimension varm(2,mx,il),vorm(2,mx,jx),polvm(2,mx,jx)
      dimension fmvar1(2,mx,2)

c   transforms spectral field vorm to varm (grid in y, spectral in x)

      do 100 j=1,iy

        j1=il+1-j

        do 20 k=1,2
        do 20 n=1,jx
        do 20 m=1,mx
        polvm(k,m,n)=vorm(k,m,n)*pol(m,n,j)
   20   continue

        do 21 k=1,2
        do 21 m=1,mx
        fmvar1(k,m,1)=0.
        fmvar1(k,m,2)=0.
   21   continue

        do 22 k=1,2
        do 22 n=1,jx,2
        do 23 m=1,mx
        fmvar1(k,m,1)=fmvar1(k,m,1)+polvm(k,m,n)
   23   continue
   22   continue

        do 24 k=1,2
        do 24 n=2,jx,2
        do 25 m=1,mx
        fmvar1(k,m,2)=fmvar1(k,m,2)+polvm(k,m,n)
   25   continue
   24   continue

        do 26 k=1,2
        do 26 m=1,mx
        varm(k,m,j1)=fmvar1(k,m,1)+fmvar1(k,m,2)
        varm(k,m,j )=fmvar1(k,m,1)-fmvar1(k,m,2)
   26   continue

	varm(2,1,j1) = 0.0
	varm(2,1,j ) = 0.0
  100 continue

      return
      end
*********************************************************************
      subroutine spec(vorg,vorm)

#include "param.h"
#include "spec.h"

      dimension vorg(ix,il),gt(ix,2)
      dimension vorm(2,mx,jx),scifm(2,mx),dcifm(2,mx),cifm(2,nwaves,2),
     * rpvm(2,mx,jx)

c  transforms grid field vorg to spectral field vorm

      do 1 k=1,2
      do 1 n=1,jx
      do 1 m=1,mx
        vorm(k,m,n)=0.0
    1 continue

      do 10 j=1,iy

        j1=il+1-j

        do 20 i=1,ix
          gt(i,1)=vorg(i,j1)*wt(j)
          gt(i,2)=vorg(i,j)*wt(j)
   20   continue

        call ncgw2(gt(1,1),gt(1,2),cifm(1,1,1),cifm(1,1,2))

        do 21 k=1,2
        do 21 m=1,mx
          scifm(k,m)=cifm(k,m,1)+cifm(k,m,2)
          dcifm(k,m)=cifm(k,m,1)-cifm(k,m,2)
   21   continue

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n) = pol(m,n,j)*scifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n)=pol(m,n,j)*dcifm(k,m)

        do 3 k=1,2
        do 3 n = 1,jtrun1
        do 3 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n)+ rpvm(k,m,n)
    3   continue

   10 continue

      do 30 n=1,jx
   30 vorm(2,1,n) = 0.0
      return
      end
**********************************************************************
      subroutine specx(vorg,varm)

#include "param.h"
#include "spec.h"

      dimension vorg(ix,il),gt(ix,2)
      dimension varm(2,mx,il),vm(2,nwaves,2)

c   transforms grid field vorg to varm (grid in y, spectral in x)

      do 1 j=1,iy

        j1=il+1-j

        do 7 i=1,ix
          gt(i,1)=vorg(i,j)
          gt(i,2)=vorg(i,j1)
    7   continue

        call ncgw2(gt(1,1),gt(1,2),vm(1,1,1),vm(1,1,2))

        do 8 k=1,2
        do 8 m=1,mx
          varm(k,m,j ) = vm(k,m,1)
          varm(k,m,j1) = vm(k,m,2)
    8   continue

      varm(2,1,j ) = 0.0
      varm(2,1,j1) = 0.0
    1 continue

      return
      end
******************************************************************
      subroutine specy(varm,vorm)

#include "param.h"
#include "spec.h"

      dimension varm(2,mx,il),vorm(2,mx,jx),scifm(2,mx),dcifm(2,mx),
     * rpvm(2,mx,jx)

c  transforms varm (grid in y, spectral in x) into spectral field vorm

      do 1 k=1,2
      do 1 n=1,jx
      do 1 m=1,mx
        vorm(k,m,n)=0.0
    1 continue

      do 10 j=1,iy

        j1=il+1-j

        do 11 k=1,2
        do 11 m=1,mx
          scifm(k,m)=(varm(k,m,j1)+varm(k,m,j))*wt(j)
          dcifm(k,m)=(varm(k,m,j1)-varm(k,m,j))*wt(j)
   11   continue

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n)=pol(m,n,j)*scifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n)=pol(m,n,j)*dcifm(k,m)

        do 12 k=1,2
        do 12 n = 1,jtrun1
        do 12 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n)+ rpvm(k,m,n)
   12   continue

   10 continue

      return
      end
****************************************************************
      subroutine dspec(ug,vg,divm,ivd,kcos)

#include "param.h"
#include "spec.h"

      dimension ug(ix,il),vg(ix,il)
      dimension divm(2,mx,jx)
      dimension gucos(ix,2),gvcos(ix,2)
      dimension rpvm(2,mx,jx),szaifm(2,mx),
     * dzaifm(2,mx),sbifm(2,mx),dbifm(2,mx),
     * aifm(2,nwaves,2),bifm(2,nwaves,2)

c  computes spectral divergence (ivd =1) or vorticity (ivd=2)
c  from grid u and v (k=2), or grid u x cos(lat), v x cos(lat) (k=1)

      do 6 k=1,2
      do 6 n=1,jx
      do 6 m=1,mx
        divm(k,m,n)=0.
    6 continue

      do 1 j=1,iy

        j1=il+1-j

        if(ivd.eq.1) then
          if(kcos.eq.1) then
            do 7 i=1,ix
              gucos(i,1)=ug(i,j1)
              gucos(i,2)=ug(i,j)
              gvcos(i,1)=vg(i,j1)
              gvcos(i,2)=vg(i,j)
    7       continue
          else
            do 77 i=1,ix
              gucos(i,1)=ug(i,j1)*coa(j)
              gucos(i,2)=ug(i,j)*coa(j)
              gvcos(i,1)=vg(i,j1)*coa(j)
              gvcos(i,2)=vg(i,j)*coa(j)
   77       continue
          endif
        else
          if(kcos.eq.1) then
            do 8 i=1,ix
              gucos(i,1)=vg(i,j1)
              gucos(i,2)=vg(i,j)
              gvcos(i,1)=-ug(i,j1)
              gvcos(i,2)=-ug(i,j)
    8       continue
          else
            do 88 i=1,ix
              gucos(i,1)=vg(i,j1)*coa(j)
              gucos(i,2)=vg(i,j)*coa(j)
              gvcos(i,1)=-ug(i,j1)*coa(j)
              gvcos(i,2)=-ug(i,j)*coa(j)
   88       continue
          endif
      endif

        call ncgw2(gucos(1,1),gucos(1,2),aifm(1,1,1),aifm(1,1,2))
        call ncgw2(gvcos(1,1),gvcos(1,2),bifm(1,1,1),bifm(1,1,2))

        do 5 k=1,2
        do 5 m=1,mx
          sbifm(k,m)=bifm(k,m,1)+bifm(k,m,2)
          dbifm(k,m)=bifm(k,m,1)-bifm(k,m,2)
    5   continue
        do 55 m=1,mx
          szaifm(1,m)=-(aifm(2,m,1)+aifm(2,m,2))
          szaifm(2,m)=aifm(1,m,1)+aifm(1,m,2)
          dzaifm(1,m)=-(aifm(2,m,1)-aifm(2,m,2))
          dzaifm(2,m)=aifm(1,m,1)-aifm(1,m,2)
   55   continue

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n) = polydx(m,n,j)*szaifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n) = polydx(m,n,j)*dzaifm(k,m)

        do 3 k=1,2
        do 3 n = 1,jtrun1
        do 3 m = 1,mx
          divm(k,m,n) = divm(k,m,n) + rpvm(k,m,n)
    3   continue

        do 122 k=1,2
        do 122 n = 1,jtrun1,2
        do 122 m = 1,mx
  122     rpvm(k,m,n) = polydy(m,n,j)*dbifm(k,m)
        do 123 k=1,2
        do 123 n=2,jtrun1,2
        do 123 m=1,mx
  123     rpvm(k,m,n) = polydy(m,n,j)*sbifm(k,m)

        do 103 k=1,2
        do 103 n = 1,jtrun1
        do 103 m = 1,mx
          divm(k,m,n) = divm(k,m,n) - rpvm(k,m,n)
  103   continue

    1 continue

      return
      end
*********************************************************************
      subroutine dspy(um,vm,divm,ivd)

#include "param.h"
#include "spec.h"

      dimension um(2,mx,il),vm(2,mx,il),divm(2,mx,jx)
      dimension rpvm(2,mx,jx),szaifm(2,mx),
     * dzaifm(2,mx),sbifm(2,mx),dbifm(2,mx)

c   starts with u and v fields that are spectral in x and grid in y
c    and produces fully spectral divergence (ivd=1) or vorticity (2)

      do 6 k=1,2
      do 6 n=1,jx
      do 6 m=1,mx
        divm(k,m,n)=0.0
    6 continue

      do 1 j=1,iy

        j1=il+1-j

        if(ivd.eq.1) then
          do 5 k=1,2
          do 5 m=1,mx
            sbifm(k,m)=vm(k,m,j1)+vm(k,m,j)
            dbifm(k,m)=vm(k,m,j1)-vm(k,m,j)
    5     continue
          do 55 m=1,mx
            szaifm(1,m)=-(um(2,m,j1)+um(2,m,j))
            szaifm(2,m)=um(1,m,j1)+um(1,m,j)
            dzaifm(1,m)=-(um(2,m,j1)-um(2,m,j))
            dzaifm(2,m)=um(1,m,j1)-um(1,m,j)
   55     continue
        else
          do 51 k=1,2
          do 51 m=1,mx
            sbifm(k,m)=-(um(k,m,j1)+um(k,m,j))
            dbifm(k,m)=-(um(k,m,j1)-um(k,m,j))
   51     continue
          do 551 m=1,mx
            szaifm(1,m)=-(vm(2,m,j1)+vm(2,m,j))
            szaifm(2,m)=vm(1,m,j1)+vm(1,m,j)
            dzaifm(1,m)=-(vm(2,m,j1)-vm(2,m,j))
            dzaifm(2,m)=vm(1,m,j1)-vm(1,m,j)
  551     continue
        endif

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n) = polydx(m,n,j)*szaifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n) = polydx(m,n,j)*dzaifm(k,m)

        do 3 k=1,2
        do 3 n = 1,jtrun1
        do 3 m = 1,mx
          divm(k,m,n) = divm(k,m,n) + rpvm(k,m,n)
    3   continue

        do 122 k=1,2
        do 122 n = 1,jtrun1,2
        do 122 m = 1,mx
  122     rpvm(k,m,n) = polydy(m,n,j)*dbifm(k,m)
        do 123 k=1,2
        do 123 n=2,jtrun1,2
        do 123 m=1,mx
  123     rpvm(k,m,n) = polydy(m,n,j)*sbifm(k,m)

        do 103 k=1,2
        do 103 n = 1,jtrun1
        do 103 m = 1,mx
          divm(k,m,n) = divm(k,m,n) - rpvm(k,m,n)
  103   continue

    1 continue

      return
      end
*********************************************************************
      subroutine vdspec(ug,vg,vorm,divm,kcos)

#include "param.h"
#include "spec.h"

      dimension ug(ix,il),vg(ix,il)
      dimension vorm(2,mx,jx),divm(2,mx,jx)
      dimension gucos(ix,2),gvcos(ix,2)
      dimension rpvm(2,mx,jx),saifm(2,mx),daifm(2,mx),szaifm(2,mx),
     * dzaifm(2,mx),sbifm(2,mx),dbifm(2,mx),szbifm(2,mx),
     * dzbifm(2,mx),aifm(2,nwaves,2),bifm(2,nwaves,2)

c  computes spectral vorticity vorm and divergence divm from
c  grid u and v (k=2) or grid u x cos(lat), v x cos(lat) (k=1)

      do 6 k=1,2
      do 6 n=1,jx
      do 6 m=1,mx
        vorm(k,m,n)=0.
        divm(k,m,n)=0.
    6 continue

      do 1 j=1,iy

        j1=il+1-j

        if(kcos.eq.1) then
          do 7 i=1,ix
            gucos(i,1)=ug(i,j1)
            gucos(i,2)=ug(i,j)
            gvcos(i,1)=vg(i,j1)
            gvcos(i,2)=vg(i,j)
    7     continue
        else
          do 77 i=1,ix
            gucos(i,1)=ug(i,j1)*coa(j)
            gucos(i,2)=ug(i,j)*coa(j)
            gvcos(i,1)=vg(i,j1)*coa(j)
            gvcos(i,2)=vg(i,j)*coa(j)
   77     continue
        endif

        call ncgw2(gucos(1,1),gucos(1,2),aifm(1,1,1),aifm(1,1,2))
        call ncgw2(gvcos(1,1),gvcos(1,2),bifm(1,1,1),bifm(1,1,2))

        do 5 k=1,2
        do 5 m=1,mx
          saifm(k,m)=aifm(k,m,1)+aifm(k,m,2)
          daifm(k,m)=aifm(k,m,1)-aifm(k,m,2)
          sbifm(k,m)=bifm(k,m,1)+bifm(k,m,2)
          dbifm(k,m)=bifm(k,m,1)-bifm(k,m,2)
    5   continue
        do 55 m=1,mx
          szaifm(1,m)=-(aifm(2,m,1)+aifm(2,m,2))
          szaifm(2,m)=aifm(1,m,1)+aifm(1,m,2)
          dzaifm(1,m)=-(aifm(2,m,1)-aifm(2,m,2))
          dzaifm(2,m)=aifm(1,m,1)-aifm(1,m,2)
          szbifm(1,m)=-(bifm(2,m,1)+bifm(2,m,2))
          szbifm(2,m)=bifm(1,m,1)+bifm(1,m,2)
          dzbifm(1,m)=-(bifm(2,m,1)-bifm(2,m,2))
          dzbifm(2,m)=bifm(1,m,1)-bifm(1,m,2)
   55   continue

        do 222 k=1,2
        do 222 n = 1,jtrun1,2
        do 222 m = 1,mx
  222     rpvm(k,m,n) = polydx(m,n,j)*szbifm(k,m)
        do 223 k=1,2
        do 223 n=2,jtrun1,2
        do 223 m=1,mx
  223     rpvm(k,m,n) = polydx(m,n,j)*dzbifm(k,m)

        do 203 k=1,2
        do 203 n = 1,jtrun1
        do 203 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n) + rpvm(k,m,n)
  203   continue

        do 322 k=1,2
        do 322 n = 1,jtrun1,2
        do 322 m = 1,mx
  322     rpvm(k,m,n) = polydy(m,n,j)*daifm(k,m)
        do 323 k=1,2
        do 323 n=2,jtrun1,2
        do 323 m=1,mx
  323     rpvm(k,m,n) = polydy(m,n,j)*saifm(k,m)

        do 303 k=1,2
        do 303 n = 1,jtrun1
        do 303 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n) + rpvm(k,m,n)
  303   continue

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n) = polydx(m,n,j)*szaifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n) = polydx(m,n,j)*dzaifm(k,m)

        do 3 k=1,2
        do 3 n = 1,jtrun1
        do 3 m = 1,mx
          divm(k,m,n) = divm(k,m,n) + rpvm(k,m,n)
    3   continue

        do 122 k=1,2
        do 122 n = 1,jtrun1,2
        do 122 m = 1,mx
  122     rpvm(k,m,n) = polydy(m,n,j)*dbifm(k,m)
        do 123 k=1,2
        do 123 n=2,jtrun1,2
        do 123 m=1,mx
  123     rpvm(k,m,n) = polydy(m,n,j)*sbifm(k,m)

        do 103 k=1,2
        do 103 n = 1,jtrun1
        do 103 m = 1,mx
          divm(k,m,n) = divm(k,m,n) - rpvm(k,m,n)
  103   continue

    1 continue

      return
      end

*********************************************************************
      subroutine vdspy(um,vm,vorm,divm)

#include "param.h"
#include "spec.h"

      dimension um(2,mx,il),vm(2,mx,il),vorm(2,mx,jx),
     * divm(2,mx,jx)
      dimension rpvm(2,mx,jx),saifm(2,mx),szaifm(2,mx),
     * daifm(2,mx),dzaifm(2,mx),sbifm(2,mx),szbifm(2,mx),
     * dbifm(2,mx),dzbifm(2,mx)

c    starts with u and v fields that are spectral in x and grid in y
c  and computes fully spectal vorticity and divergence

      do 66 k=1,2
      do 66 n=1,jx
      do 66 m=1,mx
        vorm(k,m,n)=0.0
        divm(k,m,n)=0.0
   66 continue

      do 1 j=1,iy

        j1=il+1-j

        do 50 k=1,2
        do 50 m=1,mx
          saifm(k,m)=um(k,m,j1)+um(k,m,j)
          daifm(k,m)=um(k,m,j1)-um(k,m,j)
          sbifm(k,m)=vm(k,m,j1)+vm(k,m,j)
          dbifm(k,m)=vm(k,m,j1)-vm(k,m,j)
   50   continue
        do 55 m=1,mx
          szaifm(1,m)=-(um(2,m,j1)+um(2,m,j))
          szaifm(2,m)=um(1,m,j1)+um(1,m,j)
          dzaifm(1,m)=-(um(2,m,j1)-um(2,m,j))
          dzaifm(2,m)=um(1,m,j1)-um(1,m,j)
          szbifm(1,m)=-(vm(2,m,j1)+vm(2,m,j))
          szbifm(2,m)=vm(1,m,j1)+vm(1,m,j)
          dzbifm(1,m)=-(vm(2,m,j1)-vm(2,m,j))
          dzbifm(2,m)=vm(1,m,j1)-vm(1,m,j)
   55   continue

        do 2 k=1,2
        do 2 n = 1,jtrun1,2
        do 2 m = 1,mx
    2     rpvm(k,m,n) = polydx(m,n,j)*szbifm(k,m)
        do 3 k=1,2
        do 3 n=2,jtrun1,2
        do 3 m=1,mx
    3     rpvm(k,m,n) = polydx(m,n,j)*dzbifm(k,m)

        do 4 k=1,2
        do 4 n = 1,jtrun1
        do 4 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n) + rpvm(k,m,n)
    4   continue

        do 5 k=1,2
        do 5 n = 1,jtrun1,2
        do 5 m = 1,mx
    5     rpvm(k,m,n) = polydy(m,n,j)*daifm(k,m)
        do 6 k=1,2
        do 6 n=2,jtrun1,2
        do 6 m=1,mx
    6     rpvm(k,m,n) = polydy(m,n,j)*saifm(k,m)

        do 7 k=1,2
        do 7 n = 1,jtrun1
        do 7 m = 1,mx
          vorm(k,m,n) = vorm(k,m,n) + rpvm(k,m,n)
    7   continue

        do 22 k=1,2
        do 22 n = 1,jtrun1,2
        do 22 m = 1,mx
   22     rpvm(k,m,n) = polydx(m,n,j)*szaifm(k,m)
        do 23 k=1,2
        do 23 n=2,jtrun1,2
        do 23 m=1,mx
   23     rpvm(k,m,n) = polydx(m,n,j)*dzaifm(k,m)

        do 24 k=1,2
        do 24 n = 1,jtrun1
        do 24 m = 1,mx
          divm(k,m,n) = divm(k,m,n) + rpvm(k,m,n)
   24   continue

        do 122 k=1,2
        do 122 n = 1,jtrun1,2
        do 122 m = 1,mx
  122     rpvm(k,m,n) = polydy(m,n,j)*dbifm(k,m)
        do 123 k=1,2
        do 123 n=2,jtrun1,2
        do 123 m=1,mx
  123     rpvm(k,m,n) = polydy(m,n,j)*sbifm(k,m)

        do 103 k=1,2
        do 103 n = 1,jtrun1
        do 103 m = 1,mx
          divm(k,m,n) = divm(k,m,n) - rpvm(k,m,n)
  103   continue

    1 continue

      return
      end
*********************************************************************
      subroutine global(vor,vor1,av)

#include "param.h"
#include "spec.h"

      dimension vor(ix,il),vor1(ix,il)

c computes global mean (av) and departure from mean (vor1)

      av=0.

      do 1 i=1,ix
      do 1 j=1,iy
        jj=il+1-j
        av=av+(vor(i,j)+vor(i,jj))*wt(j)
    1 continue

      av=av/(2.*float(ix))

      do 2 i=1,ix
      do 2 j=1,il
        vor1(i,j)=vor(i,j)-av
    2 continue

      return
      end
**************************************************************
      subroutine zon(gg,gz,ge)

#include "param.h"

      dimension gg(ix,il),ge(ix,il),gz(il)

c   computes zonal mean (gz) and departure from zonal mean (ge)
      do 1 j=1,il
        gz(j)=0.
        do 2 i=1,ix
          gz(j)=gz(j)+gg(i,j)
    2   continue
        gz(j)=gz(j)/float(ix)
        do 3 i=1,ix
          ge(i,j)=gg(i,j)-gz(j)
    3   continue
    1 continue

      return
      end
