      subroutine ginit

c  this subroutine computes constants in common block in const1.h
c  it must be called before using any of the other routines.
c  EXCEPT for "spinit"...call "spinit" first in order to have
c  elements in common /gauss/ defined.

#include "param.h"   
#include "const1.h"  
c      include 'gauss.h'    

      dimension hsg9(10), hsg14(15)

      data hsg9/ .00001,.050,.140,.270,.430,.600,.760,.900,.980,1.00/
      data hsg14/.00001,.0300,.0707,.1311,.2102,.3036,.4062,.5138,  
     .            .6226,.7284,.8255,.9066,.9640,.9933,1.00/
      common /gauss/ sia(iy),wt(iy),a 

      print *,'Calling ginit'
      print *,'SIA in ginit=',sia(1), sia(iy)
      print *,'in ginit, a=',a
      akap=2./7.
      g=9.81  
      rgas=287.

      pi=4.*atan(1.) 
      omega=7.292e-05  
      rad=pi/180.0  

      do 4 k=1,kxp
      if(kx.eq.9) then
	hsg(k) = hsg9(k)
      else if(kx.eq.14) then
	hsg(k) = hsg14(k)
      end if
    4 continue

      do 6 k=1,kx
        dhs(k)=hsg(k+1)-hsg(k)
        fsg(k)=0.5*(hsg(k+1)+hsg(k))
    6 continue

      do k=1,kx
	DHSR(K)=0.5/DHS(K) 
	FSGR(K)=AKAP/(2.*FSG(K))
      enddo  

      do 5 j=1,iy
        jj=il+1-j
        radang(j)=-asin(sia(j))
        radang(jj)=asin(sia(j))
    5 continue
      do 66 j=1,il
        theta(j)=radang(j)*180./pi
        gcos(j)=cos(radang(j))
        gsin(j)=sin(radang(j))
        coriol(j)=2.*omega*gsin(j)
   66 continue

      return
      end
