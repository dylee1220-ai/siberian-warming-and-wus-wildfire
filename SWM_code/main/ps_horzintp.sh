#!/bin/bash

cat > horintp.f << EOF

	program combine 

c	parameter(nx=576,ny=181,nz=14) 
c	parameter(nx=576,ny=361,nz=14)   !  for u, v, t
	parameter(nx=576,ny=361,nz=1)   !  for ps
	parameter(ix=96,il=80) 
c	parameter(nmon=12)  
	parameter(nmon=1)  

	integer yrbegin,yrend,nyr 
	parameter(yrbegin=1979,yrend=1979) 
	parameter(nyr=yrend+1-yrbegin) 

	real lon1,lat1,dx1,dy1 
	real lon2,dx2 
c	parameter(lon1=0.0,lat1=-90.0) 
c	parameter(dx1=1.25,dy1=1) 
	parameter(lon1=0.625,lat1=-89.375) 
	parameter(dx1=1.25,dy1=1.25) 
	parameter(lon2=0.0) 
	parameter(dx2=3.75) 

	real alon1(nx),alat1(ny) 
	real alon2(ix),alat2(il),gslat(il),gs(il/2)    
 
	real a0(nx,ny,nz),a1(nx,ny,nz),a2(ix,il,nz)  

	character*99 fname0,fname1,fname2 
	logical to_write 

c for horizontal interpolation
        pi=4.*atan(1.)
        rad=pi/180.

c input resolution 
        do i=1,nx 
                alon1(i)=lon1+(i-1)*dx1
        enddo 
        do j=1,ny 
                alat1(j)=lat1+(j-1)*dy1 
        enddo 

C output resolution: the Gaussian latitude for R30
        call PGAUSSL(il/2,gs)
        do j=1,il 
           if(j.le.il/2) gslat(j)=-gs(j)/rad
           if(j.gt.il/2) gslat(j)=gs(il-j+1)/rad
        end do

        do i=1,ix 
                alon2(i)=lon2+(i-1)*dx2
        enddo 
        do j=1,il 
                alat2(j)=gslat(j)  
        enddo 

c input file  

c	open(1,file='${var}_L14_intp', 
c       open(1,file='T_L14_576x361_merra.7911', 
	open(1,file=
     *'../input/PS-clm.bdat', 
     *form='unformatted',access='direct',recl=nx*ny*nz*4) 

c output file  

        open(2,file='../temp/PS_R30',
c	open(2,file='T_R30L14_merra.7911',
c	open(2,file='PS_R30_merra.7911',
     *  form='unformatted',access='direct',recl=ix*il*nz*4) 

c	print*,'fname1,fname2' 
c	print*,fname1,fname2 

	irec=0 

	if(1979.lt.yrbegin) then 
		do iyr=1979,yrbegin-1 
		print*,' iyr=',iyr
		do imon=1,nmon 
		   irec=irec+1 
		   read(1,rec=irec) a1
		enddo 
		enddo 
	endif

	iirec=0
	do iyr=yrbegin,yrend 
		monend=nmon
		print*,' iyr=',iyr 
	do imon=1,monend  
	   irec=irec+1 
	   iirec=iirec+1
	   read(1,rec=irec) a0

	   do k=1,nz 
	   do j=1,ny 
	   do i=1,nx/2 
		a1(i,j,k)=a0(i+nx/2,j,k)
		a1(i+nx/2,j,k)=a0(i,j,k)
	   enddo 
	   enddo 
	   enddo 

	   call intpR30(ny,nx,ix,il,nz,a1,a2,alon1,alat1,alat2)  

  	   write(2,rec=iirec) a2  
	enddo
	enddo

	close(1) 
	close(2) 

	stop 
	end 


      SUBROUTINE intpR30(lat,lon,nx,ny,nz,array,R30,alon,alat,gslat)
      real array(lon,lat,nz),R30(nx,ny,nz),intpx(nx,lat)
      real alon(lon),alat(lat),gslat(ny),a1d(lat)

      do k=1,nz
C       Interpolate in zonal direction
        do j=1,lat
          do i=1,nx
            clon=(i-1)*360./nx
            intpx(i,j)=fintrp(alon,array(1,j,k),lon,clon,1)
          end do
        end do
C       Interpolate in meridional direction
        do i=1,nx
          do j=1,lat
            a1d(j)=intpx(i,j)
          end do
          do j=1,ny
            R30(i,j,k)=fintrp(alat,a1d,lat,gslat(j),1)
          end do
        end do
      end do
      return
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


      SUBROUTINE PGAUSSL(M,X)
      DOUBLE PRECISION Z,Z1,P1,P2,P3,PP,EPS
      DIMENSION X(M)
      PARAMETER (EPS=3.D-14)
c   this is a slightly modified version of a program in Numerical Recipes
c
c   input:
c      m    = number of gaussian latitudes between pole and equator
c   output:
c      x(m) = sin(gaussian latitude)

      N=2*M

      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))

    1   CONTINUE

        P1=1.D0
        P2=0.D0

        DO 11 J=1,N
          P3=P2
          P2=P1
          P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   11   CONTINUE

        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP

        IF(ABS(Z-Z1).GT.EPS) GO TO 1

        X(I)=asin(Z)

   12 CONTINUE

      RETURN
      END

EOF

pgf90 horintp.f
./a.out
rm horintp.f
