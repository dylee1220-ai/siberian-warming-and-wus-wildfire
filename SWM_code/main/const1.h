
      COMMON /ACOM/ THETA(IL),CORIOL(IL),GCOS(IL),GSIN(IL),RADANG(IL),
     * HSG(KXP),DHS(KX),FSG(KX),XGEOP1(KX),XGEOP2(KX),
     * DHSR(KX),FSGR(KX),
     * TREF(KX),
     * DMP(MX,JX),VNU,DRAG(KX),DRAGU(KX),DRAGV(KX),TZDAMP(KX),TDAMP(KX),
     * RADEQ(IL,KX),UEQ(IL,KX),VEQ(IL,KX),RTAU(KX),RTAU1(KX),
     * PI,OMEGA,AKAP,G,RGAS,DELT,ROB,
     * ISTEP,NSTEP
     * ,GS(IL/2),RAD    

      COMMON /BCOM/ VOR(MX,JX,KX,2),DIV(MX,JX,KX,2),T(MX,JX,KX,2),
     * PS(MX,JX,2),PHI(MX,JX,KX),PHI1(MX,JX,KXP),ZS(MX,JX),
     * VORDT(MX,JX,KX),DIVDT(MX,JX,KX),TDT(MX,JX,KX),PSDT(MX,JX)
     *,VORI(MX,JX,KX),DIVI(MX,JX,KX),TI(MX,JX,KX),PSI(MX,JX) 
      COMPLEX VOR,DIV,T,PS,PHI,PHI1,ZS,VORDT,DIVDT,TDT,PSDT
      complex vori,divi,ti,psi  

      COMMON /CCOM/ DUMC(MX,JX,6),ZEROC(MX,JX),
     * DUMR(IX,IL,20),ZEROR(IX,IL)
      COMPLEX DUMC,ZEROC

      COMMON /DCOM/ UG(IX,IL,KX),VG(IX,IL,KX),TG(IX,IL,KX),
     * PX(IX,IL),PY(IX,IL),
     * VORG(IX,IL,KX),DIVG(IX,IL,KX),UMEAN(IX,IL),
     * VMEAN(IX,IL),DMEAN(IX,IL),SIGDT(IX,IL,KXP),PHIG(IX,IL,KX),
     * PSTAR(IX,IL),ZSTAR(IX,IL),SIGG(ix,il,kx),sigst(ix,il,kx)
      
      COMMON /FCOM/ VFORC(MX,JX,KX),DFORC(MX,JX,KX),TFORC(MX,JX,KX),
     *  PSFORC(MX,JX),HEAT(MX,JX,KX),VCLIM(MX,JX,KX),DCLIM(MX,JX,KX),
     *  TCLIM(MX,JX,KX),PSCLIM(MX,JX),  
     *  nlinvm(MX,JX,KX),nlindm(MX,JX,KX),nlintm(MX,JX,KX),
     *  nlinsm(MX,JX,KX),psstm(MX,JX)    
      complex vforc,dforc,tforc,psforc,heat,vclim,dclim,tclim,psclim
      complex nlinvm,nlindm,nlintm,nlinsm,psstm  

      common /basicstate/ uu(ix,il,kx),vxv(ix,il,kx),tt(ix,il,kx),
     *  ppx(ix,il), ppy(ix,il), ss(ix,il,kxp), DDD(ix,il,kx),
     *  vvv(ix,il,kx), uu1(ix,il),vv1(ix,il),ddd1(ix,il)
 
      common /zonalmeandamp/ aaa, bbb, aaaa, bbbb      
      real aaa, bbb, aaaa, bbbb   
