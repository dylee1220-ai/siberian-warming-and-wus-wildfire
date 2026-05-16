#!/bin/bash
for var in U V T
do
cat > verintp.f << EOF

c this code performs vertical interpolation:
c from pressure levels to sigma levels

      program verintp_tran   

c      parameter (nx=576,ny=361,nv=42)
      parameter (nx=576,ny=361,nv=19)
      parameter (nz=14)  
      integer yrbegin,yrend 
      parameter (yrbegin=1979,yrend=1979) 
c      parameter (nmons=12,mend=12) 
      parameter (nmons=1,mend=1) 
      real imiss 
c      parameter (imiss=-999) 
      parameter (imiss=1.e+15) 

      real psg(nx,ny),psg0(nx,ny)   
      real heat(nx,ny,nv),heatsig(nx,ny,nz)    
      real heat1d(nv),heat1d1(nv)    
      real sig(nv),sig1(nv),p(nv)

c for others 
c      data p/1000., 975., 950., 925., 900., 875., 850., 825.,800.,775., 
c     &	750., 725., 700., 650., 600., 550.,500., 450., 400., 350.,300., 
c     &	250., 200., 150., 100., 70., 50., 40., 30., 20., 10., 7.,5.,4.,
c     &	 3., 2., 1., 0.7, 0.5, 0.4, 0.3, 0.1/
      data p/1000., 950., 900., 850., 800., 
     &	750., 700., 600., 500., 400., 300., 
     &	200., 100., 50., 30., 10., 5., 3., 1./

      character ufm*11/'unformatted'/
      character*99 fname,fname1,fname2   
      real sigma(nz)

      real qmh(15)/.000,.030,.0707,.1311,.2102,.3036,.4062,.5138,
     &  .6226,.7284,.8255,.9066,.9640,.9933,1.000/

      do k=1,nv 
	p(k)=p(k)*100.00 
      enddo 

      do k=1,nz
	sigma(k)=0.5*(qmh(k)+qmh(k+1))
      end do

      rgas=287.
      g=9.81
      gamma=6.5E-3
      fact=g/(rgas*gamma)	

C...INPUT  
       fname1='${var}-clm.bdat' 
       open(1,file=
     &	"../input/"//fname1,    
     &      form=ufm,access='direct',recl=nx*ny*nv*4)

       fname2='PS-clm.bdat'
       open(2,file=
     &	"../input/"//fname2,  
     &      form=ufm,access='direct',recl=nx*ny*4)
      
C...OUTPUT  
       open(3,file='../temp/${var}_L14_intp',
     &      form=ufm,access='direct',recl=nx*ny*nz*4)  

C...READ IN TRANSIENTS 

	irec=0 
	do iyr=yrbegin,yrend 
		print*,'iyr=',iyr 
	monend = nmons 
	if(iyr.eq.yrend) monend = mend 
	do imon=1,monend 
	irec=irec+1 

	read(1,rec=irec) heat  
	read(2,rec=irec) psg    

c for ps: from_-179.375 -> from_0.625
	do j=1,ny 
	do i=1,nx/2
c        psg0(i,j)=psg(i+nx/2,j)
c        psg0(i+nx/2,j)=psg(i,j) 
        psg0(i+nx/2,j)=psg(i+nx/2,j)  ! heat and psg (the same writing order in zonal direction)
        psg0(i,j)=psg(i,j) 
	enddo 
	enddo 

        do j=1,ny  
        do i=1,nx

          kk=0
	  do k=1,nv
	     if(heat(i,j,k).ne.imiss) then
	        kk=kk+1
	        heat1d1(kk)=heat(i,j,k)
	        sig1(kk)=p(k)/psg0(i,j)
 	     endif
	  enddo    

C...reverse in vertical direction  
	  do k=1,kk  
	     heat1d(kk+1-k)=heat1d1(k)  
	     sig(kk+1-k)=sig1(k)  
          enddo  

	  do k=1,nz
	     nvv=kk  
	     heatsig(i,j,k)=fintrp(sig,heat1d,nvv,sigma(k),1)
          enddo  

       end do
       end do

       write(3,rec=irec) heatsig  
       enddo    !imon
       enddo	!iyr
       close(3) 

       stop
       end

      REAL FUNCTION FINTRP(X,F,NTABS,X0,N)
      INTEGER I,IMAX,IMIN,IUPP,J,N,NH,NP,NTABS
      REAL F(NTABS),P,X(NTABS),X0
      NP = N + 1
      NH = NP / 2
      IF (NP .LE. NTABS) GO TO 1
      WRITE (*,1110) NP,NTABS
      PRINT 1110,NP,NTABS
 1110 FORMAT(1X,'THE NUMBER OF POINTS REQUESTED FOR THE ',
     * 'INTERPOLATION, ',I5,',',/,1X,'IS GREATER THAN THE ',
     *'NUMBER OF POINTS AT WHICH INPUT IS SPECIFIED, ',I5)
      STOP
    1 IF ((X0 .LE. X(1)) .OR. (X0 .GE. X(NTABS))) GO TO 3
      DO 10  I = 2,NTABS
      IUPP = I
      IF (X(IUPP) .GE. X0) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IMIN = MAX0(1,IUPP-NH)
      IMAX = IMIN + N
      IF (IMAX .LE. NTABS) GO TO 2
      IMAX = NTABS
      IMIN = NTABS - N
    2 CONTINUE
      GO TO 6
    3 IF (X0 .LE. X(1)) GO TO 4
      IMIN = NTABS - N
      IMAX = NTABS
      GO TO 6
    4 IMIN = 1
      IMAX = NP
    6 CONTINUE
      FINTRP = 0.0
      DO 40  I = IMIN,IMAX
      P = F(I)
      DO 30  J = IMIN,IMAX
      IF (J .EQ. I) GO TO 30
      P = P * ((X0 - X(J)) / (X(I) - X(J)))
   30 CONTINUE
      FINTRP = FINTRP + P
   40 CONTINUE
      RETURN
      END
EOF

pgf90 verintp.f
./a.out
rm verintp.f

done
