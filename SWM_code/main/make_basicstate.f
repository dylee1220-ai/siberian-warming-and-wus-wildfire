
c-------------------------------------------------------------------
c-this program is used to make the basic state for the anomaly GCM
c-July 14, 1997.  Linhai Yu (yu@physics.uiuc.edu, 217-333-8242)
c-------------------------------------------------------------------
c input: 3-D uwind, vwind, air temperature and 2-D surface pressure
c output: basic state file for the stationary wave model
c-------------------------------------------------------------------
c to compile:
c gfortran -O make_BS.f spectral.F fftnew.F ginit.F
c-------------------------------------------------------------------
	include "param.h"  
	include "const1.h"  
c-------------------------------------------------------------

      real imiss
      parameter(imiss=-999)  

      real ps1(ix,il),puv(ix,il,kx)
      real vv(ix,il,kx) 
      real qdg(ix,il,kx)
      real ubar(il,kx),vbar(il,kx),tbar(il,kx),psbar(il)  
      complex psm(mx,jx),ppxm(mx,jx),ppym(mx,jx),vvvm(mx,jx,kx),
     *        dddm(mx,jx,kx)
      character*88 fname1,fname2,fname3,fname4,fname5        

C  THE FOLLOWING TWO SUBROUTINES MUST BE CALLED FIRST
      call spinit
      call ginit

c--read in your ps, u, v, t-----------------------------------
c---note that ps should be in logp----------

c---INPUT: read in u,v,t,qdg and ps 
c...here ps is log(surface pressure)

      fname1='../temp/U_R30L14' 
      fname2='../temp/V_R30L14'
      fname3='../temp/T_R30L14'
      fname4='../temp/PS_R30'

c      fname1='u_R30L14_merra.1std' 
c      fname2='v_R30L14_merra.1std'
c      fname3='t_R30L14_merra.1std'
c      fname4='ps_R30_merra.1std'

      open(1,file=
     & "./"
     & //fname1,  
     & form='unformatted',access='direct',recl=ix*il*kx*4)
      open(2,file=
     & "./"
     & //fname2,  
     & form='unformatted',access='direct',recl=ix*il*kx*4)
      open(3,file=
     & "./"
     & //fname3,  
     & form='unformatted',access='direct',recl=ix*il*kx*4)
      open(4,file=
     & "./"
     & //fname4,  
     & form='unformatted',access='direct',recl=ix*il*4)

      read(1,rec=1) uu  
      read(2,rec=1) vv   
      read(3,rec=1) tt 	!in Kelvin  

      read(4,rec=1) ps1 !hpa 

      do j=1,il 
      do i=1,ix 
	      ps1(i,j)=ps1(i,j)*0.01
      enddo 
      enddo 

      close(1)
      close(2)  
      close(3)  
      close(4)  

      do j=1,il 
      do i=1,ix  
      if(abs(ps1(i,j)-imiss).le.100) print*,'MissValue!!!'  
      if(ps1(i,j).le.0.0) ps1(i,j)=1  
      do k=1,kx  
	 if(abs(uu(i,j,k)-imiss).le.100) print*,'MissValue!!!'  
	 if(abs(vv(i,j,k)-imiss).le.100) print*,'MissValue!!!'  
	 if(abs(tt(i,j,k)-imiss).le.100) print*,'MissValue!!!'  
      enddo  
      enddo 
      enddo  

      do j=1,il 
      do i=1,ix 
	 ps1(i,j)=alog(ps1(i,j))  
      enddo 
      enddo  

c get zonal mean 

      if(1.eq.1) then 	!zonal mean basicstate 

      do k=1,kx 
      do j=1,il
         ubar(j,k)=0.0 
         vbar(j,k)=0.0 
         tbar(j,k)=0.0 
         do i=1,ix 
            ubar(j,k)=ubar(j,k)+uu(i,j,k)/float(ix) 
            vbar(j,k)=vbar(j,k)+vv(i,j,k)/float(ix) 
            tbar(j,k)=tbar(j,k)+tt(i,j,k)/float(ix) 
         enddo   
c Lim        do i=1,ix 
c            uu(i,j,k)=ubar(j,k)  
c            vv(i,j,k)=vbar(j,k)  
c            tt(i,j,k)=tbar(j,k)  
c Lim         enddo
       enddo 
       enddo 

       do j=1,il
          psbar(j)=0.0 
          do i=1,ix 
             psbar(j)=psbar(j)+ps1(i,j)/float(ix) 
          enddo
c Lim         do i=1,ix 
c             ps1(i,j)=psbar(j) 
c Lim         enddo
       enddo 

       endif		!end zonal mean BS if 

      print*,'1st step' 

c---calculate the zonal-asymmetric part of the basic state---

      call spec(ps1,psm)
      call grad(psm,ppxm,ppym)
      call gridd(ppxm,ppx,2)
      call gridd(ppym,ppy,2)

c-----calculate ddd and vvv----------------------

      do k=1,kx
        call vdspec(uu(1,1,k),vv(1,1,k),vvvm(1,1,k),dddm(1,1,k),2)
        call gridd(vvvm(1,1,k),vvv(1,1,k),1)
        call gridd(dddm(1,1,k),ddd(1,1,k),1)
      end do
       
       print*,'2nd step' 

c----calculate uu1, vv1, ddd1--------------------
 
      DO 6 K=1,KX
        DHS(K)=HSG(K+1)-HSG(K)
    6 CONTINUE

      DO 12 I=1,IX
      DO 12 J=1,IL
        Uu1(I,J)=0.
        Vv1(I,J)=0.
        Ddd1(I,J)=0.
   12 CONTINUE

      DO 13 K=1,KX
        DO 113 I=1,IX
        DO 113 J=1,IL
          Uu1(I,J)=uU1(I,J)+uU(I,J,K)*DHS(K)
          Vv1(I,J)=Vv1(I,J)+Vv(I,J,K)*DHS(K)
          Ddd1(I,J)=Ddd1(I,J)+Ddd(I,J,K)*DHS(K)
  113   CONTINUE
   13 CONTINUE

c----calculate ss: sigmadot. If you already have this field,
c----you can modify this program to direct read and write it
c---to the basic state but still in the position of ss and the
c---following can be neglected. 

      DO 16 I=1,IX
      DO 16 J=1,IL
        Ss(I,J,1)=0.
   16 CONTINUE

      DO 250 K=1,KX
        DO 251 I=1,IX
        DO 251 J=1,IL
          PUV(I,J,K)=(Uu(I,J,K)-Uu1(I,J))*PpX(I,J)+
     *     (Vv(I,J,K)-Vv1(I,J))*PpY(I,J)
  251   CONTINUE
  250 CONTINUE

      DO 17 K=1,KX
        DO 117 I=1,IX
        DO 117 J=1,IL
          Ss(I,J,K+1)=Ss(I,J,K)-DHS(K)*(PUV(I,J,K)+
     *     Ddd(I,J,K)-Ddd1(I,J))
  117   CONTINUE
   17 CONTINUE
      
      do 19 i=1,ix
      do 19 j=1,il
      ss(i,j,kxp)=0.
  19  continue
      print*,'3rd step'

c----now we have ss---------

c---- OUTPUT
       fname5='../temp/basicstate_3D'
c       fname5='basicstate_merra.Z_BS.2std'
       open(21,file=
     &  "./"
     & //fname5,   
     &  form='unformatted')

       write(21) uu,vv,tt,ppx,ppy,ss,ddd,vvv,uu1,vv1,ddd1

       close(21)

      stop
      end  
