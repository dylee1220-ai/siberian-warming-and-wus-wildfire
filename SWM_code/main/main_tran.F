
c|_________________________________________________________________|
c|_________________________________________________________________|
c|  NAME:         The Anomaly GCM                                  | 
c|  DESCRIPTION:  This model comprises of four part: 1)main.f,     |
c|      2) mat.f, 3)spectral.f, and 4) *.h files.                  |
c|       mat.f and spectral.f are consisted of utility subroutines |
c|     which are used in the calculation of matrix arithematics    |
c|     and spectral transformations. These two routines are        |
c|     adapted from GFDL dry GCM.                                  |
c|       *.h files declare many variables, defines constants and   |
c|     common blocks.                                              |
c|       main.f is the core part of the model, it is consisted of  |
c|     MAIN and other subroutines which initializes, integrates,   |
c|     monitors the integration, and writes out the output.  This  |
c|     part is the 'physical' part of the Anomaly GCM.             |
c|  MODEL INPUT and OUTPUT: The model inputs parameters, constants,|
c|     basic states, external forcings, orography, etc. Output     |
c|     include restart data, model normal output data, and         |
c|     monitoring data.  See comments in the following routines.   |
c|  MODEL CONFIGURATION:  R30 spectral horizontal. 14 sigma level. |
c|  COMPILE:   f77 -O (or other arguments) main.f mat.f spectral.f |
c|  DATE:       The model was developed in the 1995-1996.  Current |
c|     version was finaled in July 1996 and commented in June 1997 |
c|  CONTACT: Linhai Yu,  yu@physics.uiuc.edu (ph: 217-328-3153)      | 
c|_________________________________________________________________|
c|_________________________________________________________________|

c        main.f
c_______________________________________________________________
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The Anomaly GCM-- MAIN                                   c 
c     Other subroutines are called in this part of the codes   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INCLUDE 'param.h'
      INCLUDE 'const1.h'

c------ PARMTR intializes Legendre transforms
c------ FFTRIG initializes the NCAR FFTS

      CALL spinit   
      CALL ginit   

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+--- The following part specifies all the constants,   +++
c+--- model parameters, basic states, orography,        +++
c+--- external forcings, and initial values             +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ccc-- ALPH=0--forward; ALPH=1--backward; ALPH=0.5--centered
ccc-- ROB = Robert filter strength

      ALPH=0.5
      ROB=0.05

ccc-- read in steps per day, total integration steps, and---------
ccc-- step interval for the output.  If you use 48 steps per day--
ccc-- then the time step for the model will be 30 mins.  If you --
ccc-- also choose 48 as the IPRINT value, then you will get daily-
ccc-- output.-----------------------------------------------------
      PRINT *, '# STEPS/DAY, # STEPS, # STEPS/ANALY ='
CCC      READ*,NSTEPS,JSTEPS,IPRINT
c...run 50 days

      NSTEPS = 48
      JSTEPS = NSTEPS*60
      IPRINT = 48

      DELT=86400./FLOAT(NSTEPS)
      DELT5=DELT*.5
      DELT2=2.*DELT
      NSTEP=NSTEPS

	open(30,file='outfile1',form='unformatted')
        open(31,file='outfile2',form='formatted')
        open(80,file='outfile3',form='unformatted')

ccc----CONST: pre-computes various constants of interest, including:
ccc--------   gas and physics constants, model levels, dampings,
ccc--------   basic states, and topography, etc.-------------------
      CALL CONST

ccc----EXTENFORCING:  assign the values of external forcings, such as
ccc----diabatic heating, transient forcings for vorticity, divergence,
ccc----temperature, and surface pressure equations.------------------
      CALL EXTERNFORCING

ccc----IMPINT: initializes constants for the implicit gravity
ccc---------   wave computation.  It is assumed that implicit steps
ccc---------   are of length DELT2 and use the forward/backward
ccc---------   parameter ALPH.-------------------------------------  
      CALL IMPINT(DELT2,ALPH)

ccc----INITIALVALUE: gives the initial values----------------- 
      CALL INITIALVALUE

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+------------This part is: THE INTEGRATION PROCESS----+
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++

ccc-- In the following integration process, two subroutines
ccc-- are employed.  STEP is the subroutine that is making the 
ccc-- time integration. And subroutine ANALY is used to calculate
ccc-- and record the monitoring information (into path 31) and 
ccc-- the model output (into path 30).------------------------

ccc-- Start with a second-order accurate forward step----

      CALL ANALY(1)
      ISTEP=0
      CALL STEP(1,1,2,DELT5,0.)
      CALL STEP(1,2,2,DELT,0.)

ccc-- Then the semi-implicit leapfrog integration--------

      DO 10 I=1,JSTEPS

      if(mod(I,10).eq.0) print*,'step=  ', I
      ISTEP=I
      IF(MOD(I,25).NE.0) THEN
        IF(MOD(I,2).NE.0) THEN
	  J1=1
	  J2=2
	  J3=1
          CALL STEP(J1,J2,J3,DELT2,ALPH)
          IF(MOD(I,IPRINT).EQ.0) CALL ANALY(J3)
        ELSE
	  J1=2
	  J2=1
	  J3=2
          CALL STEP(J1,J2,J3,DELT2,ALPH)
          IF(MOD(I,IPRINT).EQ.0) CALL ANALY(J3)
        ENDIF
      ELSE
	CALL STEP(J3,J3,J2,DELT5,0.)
	CALL STEP(J3,J2,J2,DELT,0.)
        IF(MOD(I,IPRINT).EQ.0) CALL ANALY(J2)
      END IF
   10 CONTINUE
    
      CLOSE (30)
      CLOSE (31)

ccc-- Save a restart file.-----------------------------
ccc-- The restart file is writen with the final time model output of
ccc-- vorticity, divergence, temperature, and surface pressure in the
ccc-- spectral form. These data may be used in a later integration
ccc-- as the initialization.---------------------------

      WRITE (80) VOR,DIV,T,PS  
      CLOSE (80)

ccc--end of MAIN------------ 

      STOP
      END


c____________________________________________________________________      
*********************************************************************
*    CONST: pre-computes various constants of interest, including:  *
*           gas and physics constants, model levels, dampings,        *
*           basic states, and topography, etc.                      *
*********************************************************************
      SUBROUTINE CONST

      INCLUDE 'param.h'
      INCLUDE 'gauss.h'
      INCLUDE 'const1.h'
      real zstarbar(il)
      complex ZERO  
      real DUU(14),DVV(14)  

C....12-month round off average of Damping Coef
      DATA DUU/35.0, 35.0, 35.0, 35.0, 35.0, 35.0, 35.0, 
     &    35.0, 29.0, 23.0, 8.6, 1.0, .5, .3/

      DATA DVV/25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0,
     &    23.0, 14.0, 6.5, 6.0, 1.0, .3, .15/

c--- gas and physics constants------------- 
      ZERO=(0.,0.)
c      PI=4.*ATAN(1.)
c      OMEGA=2.*PI/86400.
c      AKAP=2./7.
c      G=9.81
      PREF=1.0E05
c      RGAS=AKAP*1004.
      DELH=60.

c-- linear dampings to vorticity and divergence (Rayleigh
c-- Friction DRAG) and to temperature(Newtonian Cooling RTAU)
c-- to the lowest two layers in the current model configeration.
c-- If you want to modify the layers that the dampings are 
c-- applied or the relative strength of the dampings  applied
c-- to the different layers, you need to modify the following.

      DO K=1,KX
	RTAU(K)=1.0/(15.0*86400.0)  
        DRAGU(K)=1.0/(DUU(k)*86400.0)
	DRAGV(K)=1.0/(DVV(K)*86400.0)  
      END DO

c-- Additional zonalmean dampings: aaa is applied to the 
c-- vor and div (zonalmean ug and vg), and bbb is applied to 
c-- zonalmean temperature.---------------------------------
      aaa=1./(86400.*3.)
      bbb=1./(86400.*3.)

c-----READ IN THE BASIC STATE-----------------------------
c----uu: zonal velocity.----vxv: meridional velocity.-----
c----tt: temperature.---ppx and ppy: pressure gradients.--
c----ss: vertical velocity in sigma coordinate.-----------
c----ddd: divergence.----vvv: vorticity.------------------
c----uu1, vv1, and ddd1 are respectively the vertical 
c----means of uu, vxv, and ddd.---------------------------
c---- ppx ppy uu1 vv1 ddd1 are (ix,il).-------------------
c---- ss is (ix,il,kxp) and all others are (ix,il,kx).----

c---- read in basic state

      open(1,file=
     *	'basicstate_3D',
c     *	'basicstate_merra.Z_BS.-1std-NAO',
c     *	'basicstate_merra.Z_BS.-2std',
     *	form='unformatted',convert='little_endian') 
      read(1) uu,vxv,tt,ppx,ppy,ss,ddd,vvv,uu1,vv1,ddd1
      close(1)

c-----ZEROS------------------------
      DO 2 J=1,IL
      DO 2 I=1,IX
        ZEROR(I,J)=0.
    2 CONTINUE
      DO 3 N=1,JX
      DO 3 M=1,MX
        ZEROC(M,N)=ZERO
    3 CONTINUE

c----Biharmonic Diffusion--------------------------
c----VNU id the biharmonic diffusion coefficient---
c----DMP is strength of biharmonic damping---------
      VNU=1.E17
      DO 10 M=1,MX
      DO 10 N=1,JX
        MM=ISC*(M-1)
        ELA=FLOAT(MM+N-1)
        EL2A=ELA*(ELA+1.)/(A*A)
        DMP(M,N)=VNU*EL2A*EL2A
   10 CONTINUE

c-----topography---------------------------------------------
c---In the current model configuration, orography is 
c---assigned as zero, which means no topography is included--
c---If you want to consider the orograph effect, you can ---
c---include it here.----------------------------------------
c      DO 20 I=1,IX
c      DO 20 J=1,IL
c      ZSTAR(I,J)=0.
c   20 CONTINUE

	do j=1,il
	do i=1,ix
	  zstar(i,j)=0.0  
	enddo
	enddo

      CALL SPEC(ZSTAR,ZS)

c----constants used by GEOP--------------------------------
      DO 30 K=1,KX
      XGEOP1(K)=RGAS*ALOG(HSG(K+1)/FSG(K))
      IF(K.NE.KX) XGEOP2(K+1)=RGAS*ALOG(FSG(K+1)/HSG(K+1))
   30 CONTINUE

      RETURN
      END

c_________________________________________________________________
c*****************************************************************
c*   INITIALVALUE: gives the initial values to the model         *
c*                 variables.                                    *
c*****************************************************************
      SUBROUTINE INITIALVALUE

      INCLUDE 'param.h'
      INCLUDE 'const1.h'
      COMPLEX ZERO

      ZERO=(0.,0.)

c----ask if you want to initialize the model from a restart
c----file or not.  If your command is 1, then the model will 
c--- be initialized with the spectral vorticity, divergence,
c--- temperature, and sufface pressure from the restart file
c--- of path 80.  If 2, then VOR, DIV, T, PS will first be 
c--- initialized as ZEROs with trivial perturbations.--------

c      PRINT*, 'RESTART? 1=Y, 2=N'
c      READ*,ISTART
      ISTART = 2

      IF(ISTART.EQ.1) THEN
        READ (80) VOR,DIV,T,PS
        CLOSE (80)
      ELSE
        DO  N=1,JX
        DO  M=1,MX
          PS(M,N,1)=0.0
        DO  K=1,KX
          VOR(M,N,K,1)=0.0
          DIV(M,N,K,1)=0.0
          T(M,N,K,1)=0.0
        END DO
        END DO
        END DO
  
        VOR(2,4,9,1)=CMPLX(1.E-07,0.)
        VOR(2,4,8,1)=CMPLX(1.E-07,0.)
        VOR(2,4,7,1)=CMPLX(1.E-07,0.)

        DO K=1,KX
        DO N=1,JX 
          VOR(1,N,K,1)=0.
          DIV(1,N,K,1)=0.
          T(1,N,K,1)=0.
        END DO 
        END DO 

      ENDIF

      RETURN
      END

c____________________________________________________________________
c********************************************************************
c*   EXTENFORCING:  assign the values of external forcings, such as *
c*         diabatic heating, transient forcings for vorticity,      *
c*         vergence,temperature, and surface pressure equations.    *
c********************************************************************
      SUBROUTINE EXTERNFORCING
      INCLUDE 'param.h'
      INCLUDE 'const1.h'

c--An idealized diabatic heating source is applied in current model
c--configuration with ecliptic horizontal and sine vertical distribution.
c--Other kinds of idealized or realistic heating forcing can also
c--be applied which may require the modification of the following part of
c--codes.---------------------------------

      real heatrate,field(ix,il),vert(kx)  
      real heat1(ix,il,kx),heat2(ix,il,kx)  
      real tranv1(ix,il,kx),tranv2(ix,il,kx) 
      real trand1(ix,il,kx),trand2(ix,il,kx) 
      real trant1(ix,il,kx),trant2(ix,il,kx) 
      real heatgrid(ix,il,kx),heatbar(il,kx)
      real tranv(ix,il,kx),trand(ix,il,kx),trant(ix,il,kx)
      real tranvbar(il,kx),trandbar(il,kx),trantbar(il,kx) 
      real alon(ix),alat(il) 
      real heatgrid0(ix,il,kx),heatorig(ix,il,kx) 

      do i=1,ix 
         alon(i)=(i-1)*3.75 
      enddo 
      do j=1,il 
         alat(j)=-88.875+(j-1)*2.25 
      enddo 

c construct idealized transient vorticity forcing
       lon1 =320-10/2.0
       lon2 =320+10/2.0
       lat1 =40-10/2.0
       lat2 =40+10/2.0

        open(88,file='tranv.vert_profile',
     *  form='unformatted',access='direct',recl=4)
        do k=1,kx
                read(88,rec=k) vert(k)
        enddo
        close(88)

        do k=1,kx
        do j=1,il
        do i=1,ix
           tranv(i,j,k)=0.0
           if(alon(i).ge.lon1.and.alon(i).le.lon2.and.
     *        alat(j).ge.lat1.and.alat(j).le.lat2) then
                tranv(i,j,k)=(4.0)*
     *          (sin(3.1415926*(alat(j)-lat1)/10)**2)*
     *          (sin(3.1415926*(alon(i)-lon1)/10)**2)*
     *          vert(kx+1-k)
           endif

        enddo
        enddo
        enddo
c end construct idealized transient vorticity forcing


       do k=1,14
       do j=1,il
         tranvbar(j,k)=0.0
         do i=1,ix
           tranvbar(j,k)=tranvbar(j,k)+tranv(i,j,k)/float(ix)
         enddo
       enddo
       enddo

       do k=1,14
       do j=1,il
       do i=1,ix
         heatgrid(i,j,k)=0.0
       enddo
       enddo
       enddo

	do k=1,14
	  CALL spec(heatgrid(1,1,k), heat(1,1,k))
	enddo

	do k=1,kx 
	do j=1,il 
	do i=1,ix 
                tranv(i,j,k)=tranv(i,j,k)-tranvbar(j,k)
		trand(i,j,k)=0.0 
		trant(i,j,k)=0.0 
	enddo 
	enddo 
	enddo 

	do k=1,kx
	  CALL spec(tranv(1,1,k),vforc(1,1,k))
	  CALL spec(trand(1,1,k),dforc(1,1,k))
	  CALL spec(trant(1,1,k),tforc(1,1,k))
	enddo

      do 20 i=1,mx
      do 20 j=1,jx
        psforc(i,j)=(0.,0.)
  20	continue

      return
      end




c____________________________________________________________________
c********************************************************************
c*    STEP is the subroutine that contains the the model equations. *
c*    STEP calculates the tendency terms and makes the integration. *
c********************************************************************
      SUBROUTINE STEP(J1,J2,J3,DT,ALPH)

      INCLUDE 'param.h'
      include 'gauss.h'  
      PARAMETER (IXIL=IX*IL)
      INCLUDE 'const1.h'

      real uggggg(il,kx),vggggg(il,kx),tggggg(il,kx)

      DIMENSION XCON(MX,JX),DIVGM(IX,IL,KX),DM(IX,IL),
     * SIGM(IX,IL,KXP),TGG(IX,IL,KX),PUV(IX,IL,KX),PUV2(IX,IL,KX)
     * ,PUV1(IX,IL,KX),PUV3(IX,IL,KX),PUV7(IX,IL,KX)
     * ,DTDLNQ(IX,IL,KX),PSG(IX,IL),PSG4(IX,IL)

 
      COMPLEX ZERO

      ZERO=(0.,0.)
      NDAY=ISTEP/NSTEP

      DO J=1,IL
      DO I=1,IX
      DO K=2,KX
	 DTDLNQ(I,J,K)=(TT(I,J,K)-TT(I,J,K-1))/
     *                (ALOG(FSG(K))-ALOG(FSG(K-1)))
      END DO
      DTDLNQ(I,J,1)=DTDLNQ(I,J,2)
      END DO
      END DO



      J4=J1
      IF(ALPH.EQ.0) J4=J2

      DO 10 K=1,KX

        CALL UVSPEC(VOR(1,1,K,J2),DIV(1,1,K,J2),
     *   DUMC(1,1,1),DUMC(1,1,2))

        CALL GRIDD(DUMC(1,1,1),UG(1,1,K),2)
        CALL GRIDD(DUMC(1,1,2),VG(1,1,K),2)
        CALL GRIDD(VOR(1,1,K,J2),VORG(1,1,K),1)
        CALL GRIDD(DIV(1,1,K,J2),DIVG(1,1,K),1)
        CALL GRIDD(DIV(1,1,K,J4),DIVGM(1,1,K),1)
        CALL GRIDD(T(1,1,K,J2),TG(1,1,K),1)
   10 CONTINUE 

c--------hydrostatic balance------------------------
      CALL GEOP(J4)
c---------------------------------------------------

      DO 12 I=1,IX
      DO 12 J=1,IL
        UMEAN(I,J)=0.
        VMEAN(I,J)=0.
        DMEAN(I,J)=0.
        DM(I,J)=0.
   12 CONTINUE

      DO 13 K=1,KX
        DO 113 I=1,IX
        DO 113 J=1,IL
          UMEAN(I,J)=UMEAN(I,J)+UG(I,J,K)*DHS(K)
          VMEAN(I,J)=VMEAN(I,J)+VG(I,J,K)*DHS(K)
          DMEAN(I,J)=DMEAN(I,J)+DIVG(I,J,K)*DHS(K)
          DM(I,J)=DM(I,J)+DIVGM(I,J,K)*DHS(K)
  113   CONTINUE

   13 CONTINUE

c-----------------pressure tendency equation--------------

      PS(1,1,J2)=ZERO


      CALL GRAD(PS(1,1,J2),DUMC(1,1,2),DUMC(1,1,3))
      CALL GRIDD(DUMC(1,1,2),PX,2)
      CALL GRIDD(DUMC(1,1,3),PY,2)

      CALL GRIDD(PS(1,1,J2),PSG,1)
      CALL LAP(PS(1,1,J2),DUMC(1,1,2))
      CALL LAP(DUMC(1,1,2),DUMC(1,1,3))
      CALL GRIDD(DUMC(1,1,3),PSG4,1)


      DO 14 I=1,IX
      DO 14 J=1,IL
        DUMR(I,J,1)=UMEAN(I,J)*PX(I,J)+VMEAN(I,J)*PY(I,J)
        DUMR(I,J,1)=DUMR(I,J,1)+UMEAN(I,J)*PPX(I,J)+VMEAN(I,J)*PPY(I,J)
        DUMR(I,J,1)=DUMR(I,J,1)+UU1(I,J)*PX(I,J)+VV1(I,J)*PY(I,J) 
        DUMR(I,J,2)=-DUMR(I,J,1)-DM(I,J)
   14 CONTINUE

      CALL SPEC(DUMR(1,1,2),PSDT)
c----add the transient surface pressure forcing---     
      do i=1,mx
      do j=1,jx
      psdt(i,j)=psdt(i,j)+psforc(i,j)
      end do
      end do

      PSDT(1,1)=ZERO

c---compute sigma-dot "vertical" velocity in two ways, 
c----    using divergence from j2 and j1 time steps

      DO 16 I=1,IX
      DO 16 J=1,IL
        SIGDT(I,J,1)=0.
        SIGM(I,J,1)=0.
   16 CONTINUE

      DO 250 K=1,KX
        DO 251 I=1,IX
        DO 251 J=1,IL
          PUV1(I,J,K)=(UG(I,J,K)-UMEAN(I,J))*PX(I,J)+
     *     (VG(I,J,K)-VMEAN(I,J))*PY(I,J)
          PUV2(I,J,K)=(UU(I,J,K)-UU1(I,J))*PX(I,J)+
     *     (VXV(I,J,K)-VV1(I,J))*PY(I,J)
          PUV3(I,J,K)=(UG(I,J,K)-UMEAN(I,J))*PPX(I,J)+
     *     (VG(I,J,K)-VMEAN(I,J))*PPY(I,J)
          PUV(I,J,K)=PUV1(I,J,K)+PUV2(I,J,K)+PUV3(I,J,K)
          PUV7(I,J,K)=(UU(I,J,K)-UU1(I,J))*PPX(I,J)+
     *     (VXV(I,J,K)-VV1(I,J))*PPY(I,J)
  251   CONTINUE

  250 CONTINUE

      DO 17 K=1,KX
        DO 117 I=1,IX
        DO 117 J=1,IL
          SIGDT(I,J,K+1)=SIGDT(I,J,K)-DHS(K)*(PUV(I,J,K)+
     *     DIVG(I,J,K)-DMEAN(I,J))
          SIGM(I,J,K+1)=SIGM(I,J,K)-DHS(K)*(PUV(I,J,K)+
     *     DIVGM(I,J,K)-DM(I,J))
  117   CONTINUE
   17 CONTINUE
c------------------------------------------------
      DO 119 K=1,KX
      DO 118 I=1,IX
      DO 118 J=1,IL
        TGG(I,J,K)=TG(I,J,K)-TREF(K)
  118 CONTINUE
  119 CONTINUE
c------------------

      DO 18 K=1,KX
c-- verical advection ------
        IF (K.EQ.1) THEN
          DO 19 J=1,IL
          DO 19 I=1,IX
            DUMR(I,J,1)=SIGDT(I,J,2)*(UG(I,J,2)-UG(I,J,1))*DHSR(K)
     *                 +SIGDT(I,J,2)*(UU(I,J,2)-UU(I,J,1))*DHSR(K)
     *                 +SS(I,J,2)*(UG(I,J,2)-UG(I,J,1))*DHSR(K)
       
            DUMR(I,J,2)=SIGDT(I,J,2)*(VG(I,J,2)-VG(I,J,1))*DHSR(K)
     *                 +SIGDT(I,J,2)*(VXV(I,J,2)-VXV(I,J,1))*DHSR(K)
     *                 +SS(I,J,2)*(VG(I,J,2)-VG(I,J,1))*DHSR(K)

            DUMR(I,J,3)=SIGDT(I,J,2)*(TGG(I,J,2)-TGG(I,J,1)+
     *                   TT(I,J,2)-TT(I,J,1))
     *             +SS(I,J,2)*(TG(I,J,2)-TG(I,J,1))
     *             +SIGM(I,J,2)*(TREF(2)-TREF(1))
            DUMR(I,J,3)=DUMR(I,J,3)*DHSR(K)
   19     CONTINUE
        ELSE IF(K.EQ.KX) THEN
          DO 199 J=1,IL
          DO 199 I=1,IX
            DUMR(I,J,1)=SIGDT(I,J,KX)*(UG(I,J,KX)-UG(I,J,KXM))*DHSR(K)
     *                 +SIGDT(I,J,KX)*(UU(I,J,KX)-UU(I,J,KXM))*DHSR(K)
     *                 +SS(I,J,KX)*(UG(I,J,KX)-UG(I,J,KXM))*DHSR(K)      
            DUMR(I,J,2)=SIGDT(I,J,KX)*(VG(I,J,KX)-VG(I,J,KXM))*DHSR(K)
     *                 +SIGDT(I,J,KX)*(VXV(I,J,KX)-VXV(I,J,KXM))*DHSR(K)
     *                 +SS(I,J,KX)*(VG(I,J,KX)-VG(I,J,KXM))*DHSR(K)

            DUMR(I,J,3)=SIGDT(I,J,KX)*(TGG(I,J,KX)-TGG(I,J,KXM)
     *                                +TT(I,J,KX)-TT(I,J,KXM))
     *             +SS(I,J,KX)*(TG(I,J,KX)-TG(I,J,KXM))      
     *             +SIGM(I,J,KX)*(TREF(KX)-TREF(KXM))
            DUMR(I,J,3)=DUMR(I,J,3)*DHSR(K)
  199     CONTINUE

        ELSE
          DO 1999 J=1,IL
          DO 1999 I=1,IX
            DUMR(I,J,1)=(SIGDT(I,J,K+1)*(UG(I,J,K+1)-UG(I,J,K))
     *           +SIGDT(I,J,K)*(UG(I,J,K)-UG(I,J,K-1)))*DHSR(K)
     *           + (SIGDT(I,J,K+1)*(UU(I,J,K+1)-UU(I,J,K))
     *           +SIGDT(I,J,K)*(UU(I,J,K)-UU(I,J,K-1)))*DHSR(K)
     *           + (SS(I,J,K+1)*(UG(I,J,K+1)-UG(I,J,K))
     *           +SS(I,J,K)*(UG(I,J,K)-UG(I,J,K-1)))*DHSR(K)
            DUMR(I,J,2)=(SIGDT(I,J,K+1)*(VG(I,J,K+1)-VG(I,J,K))
     *           +SIGDT(I,J,K)*(VG(I,J,K)-VG(I,J,K-1)))*DHSR(K)
     *           + (SIGDT(I,J,K+1)*(VXV(I,J,K+1)-VXV(I,J,K))
     *           +SIGDT(I,J,K)*(VXV(I,J,K)-VXV(I,J,K-1)))*DHSR(K)
     *           + (SS(I,J,K+1)*(VG(I,J,K+1)-VG(I,J,K))
     *           +SS(I,J,K)*(VG(I,J,K)-VG(I,J,K-1)))*DHSR(K)

            DUMR(I,J,3)=SIGDT(I,J,K+1)*(TGG(I,J,K+1)-TGG(I,J,K))
     *           +SIGDT(I,J,K)*(TGG(I,J,K)-TGG(I,J,K-1))
     *           +SIGDT(I,J,K+1)*(TT(I,J,K+1)-TT(I,J,K))
     *           +SIGDT(I,J,K)*(TT(I,J,K)-TT(I,J,K-1))
     *           +SS(I,J,K+1)*(TG(I,J,K+1)-TG(I,J,K))
     *           +SS(I,J,K)*(TG(I,J,K)-TG(I,J,K-1))
     *           +SIGM(I,J,K+1)*(TREF(K+1)-TREF(K))
     *           +SIGM(I,J,K)*(TREF(K)-TREF(K-1))
            DUMR(I,J,3)=DUMR(I,J,3)*DHSR(K)
 1999     CONTINUE
        ENDIF

c--------calculate the zonalmean of VG and UG-------------
c--uggggg(j,k) and vggggg(j,k) are respectively the ------
c--zonalmean parts of UG and VG--------------------------
      do j=1,il
      uggggg(j,k)=0.0
      vggggg(j,k)=0.0
      end do

      do i=1,ix
      do j=1,il
      uggggg(j,k)=uggggg(j,k)+UG(I,J,K)/float(ix)
      vggggg(j,k)=vggggg(j,k)+VG(I,J,K)/float(ix)
      end do
      end do
c----------------

        DO 20 J=1,IL
        DO 20 I=1,IX
          DUMR(I,J,4)=VG(I,J,K)*VORG(I,J,K)
     *     -RGAS*TGG(I,J,K)*PX(I,J)-DUMR(I,J,1)
     *      -DRAGU(K)*UG(i,J,K)
     *      -aaa*uggggg(j,k)
          DUMR(I,J,4)=Vxv(I,J,K)*VORG(I,J,K)
     *     -RGAS*TG(I,J,K)*PPX(I,J)+DUMR(I,J,4)
          DUMR(I,J,4)=VG(I,J,K)*(VVV(I,J,K)+CORIOL(j))
     *     -RGAS*TT(I,J,K)*PX(I,J)+DUMR(I,J,4)
     
          DUMR(I,J,5)=-UG(I,J,K)*VORG(I,J,K)
     *     -RGAS*TGG(I,J,K)*PY(I,J)-DUMR(I,J,2)
     *      -DRAGV(K)*VG(I,J,K)
     *      -aaa*vggggg(j,k)
          DUMR(I,J,5)=-UU(I,J,K)*VORG(I,J,K)
     *     -RGAS*TG(I,J,K)*PPY(I,J)+DUMR(I,J,5)
          DUMR(I,J,5)=-UG(I,J,K)*(VVV(I,J,K)+coriol(j))
     *     -RGAS*TT(I,J,K)*PY(I,J)+DUMR(I,J,5)

          DUMR(I,J,6)=0.5*(UG(I,J,K)*UG(I,J,K)+
     *                  VG(I,J,K)*VG(I,J,K))
     *        +UG(I,J,K)*UU(I,J,K)+VG(I,J,K)*VXV(I,J,K)
          DUMR(I,J,7)=UG(I,J,K)*(TGG(I,J,K)+TT(I,J,K))
     *                +UU(I,J,K)*TG(I,J,K)
          DUMR(I,J,8)=VG(I,J,K)*(TGG(I,J,K)+TT(I,J,K))
     *                +VXV(I,J,K)*TG(I,J,K)
          DUMR(I,J,9)=(TGG(I,J,K)+TT(I,J,K))*DIVG(I,J,K)
     *              +TG(I,J,K)*DDD(I,J,K)-DUMR(I,J,3)
     *      +FSGR(K)*(TT(I,J,K)+TGG(I,J,K))*(SIGDT(I,J,K+1)
     *                                        +SIGDT(I,J,K))
     *      +FSGR(k)*TG(I,J,K)*(SS(I,J,K+1)+SS(I,J,K)) 
     *      +FSGR(K)*TREF(K)*(SIGM(I,J,K+1)+SIGM(I,J,K))
     *      +AKAP*(
     *            TG(I,J,K)*PUV1(I,J,K)
     *           +TG(I,J,K)*PUV2(I,J,K)
     *           +TG(I,J,K)*PUV3(I,J,K)            
     *           +TT(I,J,K)*PUV1(I,J,K)
     *           +TT(I,J,K)*PUV2(I,J,K)
     *           +TT(I,J,K)*PUV3(I,J,K) 
     *           +TG(I,J,K)*PUV7(I,J,K)
     *           -TG(I,J,K)*DDD1(I,J)
     *           -(TT(I,J,K)+TGG(I,J,K))*DMEAN(I,J)
     *           -TREF(K)*DM(I,J))
   20   CONTINUE

c----------calculate the zonalmean of TG----------
c----------tggggg(j,k) is the zonalmean part of TG(i,j,k)--
      do j=1,il
      tggggg(j,k)=0.0
      end do

      do i=1,ix
      do j=1,il
      tggggg(j,k)=tggggg(j,k)+TG(i,j,k)/float(ix)
      end do
      end do
c-----------------------------------
       do 220 j=1,il 
       DO 220 I=1,IX
c       DUMR(I,J,9)=DUMR(I,J,9)-TG(I,J,K)*RTAU(K)-tggggg(j,k)*bbb
       DUMR(I,J,9)=DUMR(I,J,9)-RTAU(K)*
     *    (TG(I,J,K)-DTDLNQ(I,J,K)*PSG(I,J))-tggggg(j,k)*bbb
     *    +VNU*DTDLNQ(I,J,K)*PSG4(I,J)

  220  CONTINUE

c----calculate divdt and vordt---        
        CALL VDSPEC(DUMR(1,1,4),DUMR(1,1,5),VORDT(1,1,K),
     *   DIVDT(1,1,K),2)

        CALL SPEC(DUMR(1,1,6),DUMC(1,1,3))
        DO 221 N=1,JX
        DO 221 M=1,MX
          DUMC(M,N,3)=DUMC(M,N,3)+PHI(M,N,K)
     *     +RGAS*TREF(K)*PS(M,N,J4)
  221 CONTINUE

        CALL LAP(DUMC(1,1,3),DUMC(1,1,4))

        DO 222 N=1,JX
        DO 222 M=1,MX
          VORDT(M,N,K)=VORDT(M,N,K)
     *                 +vforc(m,n,k)
          DIVDT(M,N,K)=DIVDT(M,N,K)-DUMC(M,N,4)
     *                 +dforc(m,n,k)  
  222 CONTINUE
c-----calculate tdt---------

        CALL VDSPEC(DUMR(1,1,7),DUMR(1,1,8),DUMC(1,1,1),DUMC(1,1,2),2)
        CALL SPEC(DUMR(1,1,9),DUMC(1,1,3))

        DO 22 N=1,JX
        DO 22 M=1,MX
          TDT(M,N,K)=DUMC(M,N,3)-DUMC(M,N,2)
     *              +heat(m,n,k)
     *              +tforc(m,n,k)
   22   CONTINUE
c----------------------------
   18 CONTINUE

c----------Applying implicit correction ------------
      IF(ALPH.GT.0.) then
         CALL IMPLIC
      ENDIF  

c  damping
      DO 339 N=1,JX
      DO 339 M=1,MX
          XCON(M,N)=1./(1.+DMP(M,N)*DT)
  339 CONTINUE

      DO 318 K=1,KX
        DO 21 N=1,JX
        DO 21 M=1,MX
          VORDT(M,N,K)=(VORDT(M,N,K)-DMP(M,N)*VOR(M,N,K,J1)
     *      )*XCON(M,N)
          DIVDT(M,N,K)=(DIVDT(M,N,K)-DMP(M,N)*DIV(M,N,K,J1)
     *      )*XCON(M,N)
          TDT(M,N,K)=(TDT(M,N,K)-DMP(M,N)*T(M,N,K,J1)
     *      )*XCON(M,N)
   21   CONTINUE

  318 CONTINUE

c--------- Applying Robert filter--------

c      print*,'in STEP: ROB=',ROB  

      ROBB=1.-2.*ROB
      DO 89 K=1,KX
	DO 88 N=1,JX
        DO 88 M=1,MX
            DUMC(M,N,1)=VOR(M,N,K,J1)
            DUMC(M,N,2)=DIV(M,N,K,J1)
            DUMC(M,N,3)=T(M,N,K,J1)
            VOR(M,N,K,J3)=VOR(M,N,K,J1)+DT*VORDT(M,N,K)
            DIV(M,N,K,J3)=DIV(M,N,K,J1)+DT*DIVDT(M,N,K)
            T(M,N,K,J3)=T(M,N,K,J1)+DT*TDT(M,N,K)
            VOR(M,N,K,J2)=ROB*(DUMC(M,N,1)+VOR(M,N,K,J3))
     *       +ROBB*VOR(M,N,K,J2)
            DIV(M,N,K,J2)=ROB*(DUMC(M,N,2)+DIV(M,N,K,J3))
     *       +ROBB*DIV(M,N,K,J2)
            T(M,N,K,J2)=ROB*(DUMC(M,N,3)+T(M,N,K,J3))
     *       +ROBB*T(M,N,K,J2)
   88   CONTINUE
   89 CONTINUE

      DO 90 N=1,JX
      DO 90 M=1,MX
        DUMC(M,N,1)=PS(M,N,J1)
        PS(M,N,J3)=PS(M,N,J1)+DT*PSDT(M,N)
        PS(M,N,J2)=ROB*(DUMC(M,N,1)+PS(M,N,J3))
     *   +ROBB*PS(M,N,J2)
   90 CONTINUE
c----------------
      RETURN
      END   


c__________________________________________________________________
c******************************************************************
c*  GEOP: this subroutine is called by STEP to compute spectral   *
c*        geopotential from spectral temperature T and spectral   *
c*        topography ZS with the hydrostatic balance equation.    *
c******************************************************************  
      SUBROUTINE GEOP(JJ)

      INCLUDE 'param.h'
      INCLUDE 'const1.h'

c      print*,'in GEOP G=',G  
c      print*,'in GEOP KXM=',KXM  

      DO 1 N=1,JX
      DO 1 M=1,MX
        PHI(M,N,KX)=G*ZS(M,N)+XGEOP1(KX)*T(M,N,KX,JJ)
    1 CONTINUE

      DO 2 KK=1,KXM
        K=KX-KK
        DO 3 N=1,JX
        DO 3 M=1,MX
          PHI(M,N,K)=PHI(M,N,K+1)+XGEOP2(K+1)*T(M,N,K+1,JJ)
     *       +XGEOP1(K)*T(M,N,K,JJ)
    3   CONTINUE
    2 CONTINUE

      RETURN
      END

c___________________________________________________________________
c*******************************************************************
c*  subroutine ANALY is used to calculate (analyse) and record     *
c*  several monitoring variables (into path 31) and model output   *
c*  (into path 30).                                                *
c*******************************************************************
      SUBROUTINE ANALY(JJ)

      INCLUDE 'param.h'
      INCLUDE 'const1.h'

      COMPLEX TEMP(MX,JX)
      COMPLEX VOR1(MX,JX,KX),DIV1(MX,JX,KX),T1(MX,JX,KX),PS1(MX,JX)
      DIMENSION EK(KX)

c------------ The integration is monitored by recording several
c------------ variables into a record file of path (31).
c------------ The record file name ( path 31) is defined in MAIN.  

      WRITE (31,*) ISTEP     
 
cccc _____ calculate and write EK__kinetic energy___to path 31__

      DO 1 K=1,KX
        CALL INVLAP(VOR(1,1,K,JJ),TEMP)
        EK(K)=0.
        DO 2 M=2,MX
        DO 2 N=1,JX
          EK(K)=EK(K)-TEMP(M,N)*CONJG(VOR(M,N,K,JJ))
    2   CONTINUE
    1 CONTINUE

      WRITE(31,2000) (EK(K),K=1,KX)
c      print*,'EK in ANALY: EK=', EK  
 2000 FORMAT(10(1X,E12.5))
cccc __calculate and write regional VOR,DIV,T,and PS__to path 31__


      CALL UVSPEC(VOR(1,1,5,JJ),DIV(1,1,5,JJ),
     * DUMC(1,1,1),DUMC(1,1,2))

      CALL GRIDD(DUMC(1,1,1),DUMR(1,1,1),2)
      CALL GRIDD(DUMC(1,1,2),DUMR(1,1,2),2)
      CALL GRIDD(T(1,1,10,JJ),DUMR(1,1,3),1)

      WRITE(31,1000) (DUMR(1,J,1),J=1,40,4)
      WRITE(31,1000) (DUMR(1,J,2),J=1,40,4)
      WRITE(31,1000) (DUMR(1,J,3),J=1,40,4)
 1000 FORMAT(1X,10(1X,E12.5))

c------------ The following is the recording process of the model  
c------------ output.  The output is the spectral form vorticity,
c------------ divergence, temperature, and surface pressure.
c------------ The model output is writen into a file of path 30.
c------------ The out put file name is defined in MAIN

      DO M=1,MX
      DO N=1,JX
      DO K=1,KX
      VOR1(M,N,K)=VOR(M,N,K,JJ)
      DIV1(M,N,K)=DIV(M,N,K,JJ)
      T1(M,N,K)=T(M,N,K,JJ)
      END DO
      PS1(M,N)=PS(M,N,JJ)
      END DO
      END DO

c---add by hwang1
        CALL GEOP(JJ)
c---end add by hwang1

      WRITE(30) VOR1,DIV1,T1,PHI,PS1
c---------------------------------
      RETURN
      END

c_____________________________________________________________
c*************************************************************
c* IMPLIC is called in subroutine STEP to calculate the      *
c* implicit  correction.                                     *
c*************************************************************
      SUBROUTINE IMPLIC

      INCLUDE 'param.h'
      PARAMETER (MXJXKX=MX*JX*KX)
      INCLUDE 'const1.h'
      INCLUDE 'const2.h'

      COMPLEX YE(MX,JX,KX),YF(MX,JX,KX),ZERO

      ZERO=(0.,0.)

      DO 1 K=1,KX
      DO 1 N=1,JX
      DO 1 M=1,MX
        YE(M,N,K)=ZERO
    1 CONTINUE

      DO 2 K1=1,KX
      DO 2 K=1,KX
        DO 11 N=1,JX
        DO 11 M=1,MX
          YE(M,N,K)=YE(M,N,K)+XD(K,K1)*TDT(M,N,K1)
   11   CONTINUE
    2 CONTINUE

      DO 21 K=1,KX
      DO 22 N=1,JX
      DO 22 M=1,MX
        YE(M,N,K)=YE(M,N,K)+TREF1(K)*PSDT(M,N)
   22 CONTINUE
   21 CONTINUE

      DO 4 K=1,KX
        DO 44 N=1,JX
        DO 44 M=1,MX
          YF(M,N,K)=DIVDT(M,N,K)+ELZ(M,N)*YE(M,N,K)
   44   CONTINUE
    4 CONTINUE

      DO 5 K=1,KX
      DO 5 N=1,JX
      DO 5 M=1,MX
        DIVDT(M,N,K)=ZERO
    5 CONTINUE

      DO 6 N=1,JX
      DO 6 M=1,MX
        MM=ISC*(M-1)+1
        LL=MM+N-2
        IF(LL.NE.0) THEN
          DO 66 K1=1,KX
          DO 66 K=1,KX
            DIVDT(M,N,K)=DIVDT(M,N,K)+XJ(K,K1,LL)*YF(M,N,K1)
   66     CONTINUE
        ENDIF
    6 CONTINUE

      DO 46 K=1,KX
        DO 47 N=1,JX
        DO 47 M=1,MX
          PSDT(M,N)=PSDT(M,N)-DIVDT(M,N,K)*DHSX(K)
   47   CONTINUE
   46 CONTINUE

      DO 7 K=1,KX
      DO 7 K1=1,KX
        DO 48 N=1,JX
        DO 48 M=1,MX
          TDT(M,N,K)=TDT(M,N,K)+XC(K,K1)*DIVDT(M,N,K1)
   48   CONTINUE
    7 CONTINUE

      RETURN
      END

c_______________________________________________________________
c***************************************************************
c*  IMPINT initializes constants for the implicit gravity      *
c*  wave computation.  DELT2 is the length of Implicit step    *
c*  ALPH is the forward/backward paprameter.                   *
c***************************************************************
      SUBROUTINE IMPINT(DT,ALPH)

      INCLUDE 'param.h'
      INCLUDE 'gauss.h'
      INCLUDE 'const1.h'
      INCLUDE 'const2.h'

      DIMENSION DSUM(KX),INDX(KX),YA(KX,KX)

c      print*,'---***---'  
c      print*,'in IMPINT a=',a  

      XI=DT*ALPH
      XXI = XI/(A*A)

      DO 201 K=1,KX
        DHSX(K)=XI*DHS(K)
  201 CONTINUE

      do 60 n=1,jx
      do 60 m=1,mx      
        MN=M+N
        MM=ISC*(M-1)+1
        LL=MM+N-2
        elz(m,n)=FLOAT(LL)*FLOAT(LL+1)*XXI
   60 continue

c  define reference physicsphere, function of height only

      DO 101 K=1,KX
        TREF(K)=255.
        TREF1(K)=RGAS*TREF(K)
  101 CONTINUE

c   T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')

      DO 1 K=1,KX
      DO 1 K1=1,KXM
        XA(K,K1)=0.
    1 CONTINUE

      DO 2 K=1,KX
      DO 2 K1=1,KX
        YA(K,K1)=-AKAP*TREF(K)*DHS(K1)
    2 CONTINUE

      DO 3 K=2,KX
        XA(K,K-1)=0.5*(AKAP*TREF(K)/FSG(K)
     *      -(TREF(K)-TREF(K-1))/DHS(K))
    3 CONTINUE

      DO 4 K=1,KXM
        XA(K,K)=0.5*(AKAP*TREF(K)/FSG(K)
     *    -(TREF(K+1)-TREF(K))/DHS(K))
    4 CONTINUE

c     SIG(K)=XB(K,K')*D(K')

      DSUM(1)=DHS(1)
      DO 6 K=2,KX
        DSUM(K)=DSUM(K-1)+DHS(K)
    6 CONTINUE

      DO 5 K=1,KXM
      DO 5 K1=1,KX
        XB(K,K1)=DHS(K1)*DSUM(K)
        IF(K1.LE.K) XB(K,K1)=XB(K,K1)-DHS(K1)
    5 CONTINUE

c     T(K)=TEX(K)+XC(K,K')*D(K')

      DO 7 K=1,KX
      DO 7 K1=1,KX
        XC(K,K1)=YA(K,K1)
        DO 77 K2=1,KXM
          XC(K,K1)=XC(K,K1)+XA(K,K2)*XB(K2,K1)
   77   CONTINUE
    7 CONTINUE

c      P(K)=XD(K,K')*T(K') 

      DO 8 K=1,KX
      DO 8 K1=1,KX
        XD(K,K1)=0.
    8 CONTINUE

      DO 91 K=1,KX
      DO 91 K1=K+1,KX
        XD(K,K1)=RGAS*ALOG(HSG(K1+1)/HSG(K1))
   91 CONTINUE
      DO 92 K=1,KX
        XD(K,K)=RGAS*ALOG(HSG(K+1)/FSG(K))
   92 CONTINUE

c      P(K)=YE(K)+XE(K,K')*D(K')

      DO 10 K=1,KX
      DO 10 K1=1,KX
        XE(K,K1)=0.
      DO 10 K2=1,KX
        XE(K,K1)=XE(K,K1)+XD(K,K2)*XC(K2,K1)
   10 CONTINUE

      DO 11 L=1,LMAX
        XXX=(FLOAT(L)*FLOAT(L+1))/(A*A)
        DO 12 K=1,KX
        DO 12 K1=1,KX
          XF(K,K1,L)=XI*XI*XXX*(RGAS*TREF(K)*DHS(K1)-XE(K,K1))
   12   CONTINUE
        DO 13 K=1,KX
          XF(K,K,L)=XF(K,K,L)+1.
   13   CONTINUE
   11 CONTINUE

      DO 30 L=1,LMAX
        CALL INV(XF(1,1,L),XJ(1,1,L),INDX,KX)
   30 CONTINUE

      DO 50 K=1,KX
      DO 50 K1=1,KX
        XC(K,K1)=XC(K,K1)*XI
   50 CONTINUE

      RETURN
      END
