c#ifdef r15
c      PARAMETER ( JTRUN=15, MTRUN=15, IX=48, IY=20 )   
c#endif
c#ifdef r30
      PARAMETER ( JTRUN=30, MTRUN=30, IX=96, IY=40 )   
c#endif
      PARAMETER ( ISC=1)
      PARAMETER ( JX=JTRUN+2, MX=MTRUN+1 )
      PARAMETER ( IL=2*IY, JTRUN1=JTRUN+1 )
      PARAMETER ( MXP=ISC*MTRUN+1, JXP=JX+1 )
      PARAMETER ( NWAVES=1+IX/2)
      PARAMETER ( LMAX=MXP+JX-2)  
c#ifdef k9
c      PARAMETER (KX=9)
c#endif
c#ifdef k14
      PARAMETER (KX=14)
c#endif
c#ifdef k20
c      PARAMETER (KX=20)
c#endif
      PARAMETER (KXM=KX-1,KXP=KX+1)
