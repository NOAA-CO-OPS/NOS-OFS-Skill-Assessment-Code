C      This subroutine reads a tidal time series and find the MHW, MLW,
C      Using SVD, or find events for nontidal time series.
C      Time in hours

      SUBROUTINE EXTREMES(TTMP,HORI,NMAX1,IPR1,DELHR1,DELAMP1,
     1  DELPCT1,IOPTA1,THIGHS0,HHIGHS0,IDX0,NSMAX0,TCUT,DELTT,
     2  ZHALL,THALL,NUM_H,ZLALL,TLALL,NUM_L,NTYPE)

      PARAMETER (NMX=735*240,NMX2=NMX/10)
      REAL*8 DELTT
      COMMON/B1/H(NMX),T(NMX),NMAX,NPERHR,
     1 THIGHS(NMX2),HHIGHS(NMX2),IHH(NMX2),ILL(NMX2),IDX(NMX2),
     2 NSMAX,DELHR,DELAMP,DELPCT,IOPTA,HMAX,HMIN

      DIMENSION HH(NMX),TH1(NMX2),ZH1(NMX2),TH2(NMX2),
     1   ZH2(NMX2),TTMP(NMX),HTMP(NMX),HORI(NMX),
     2   ZHALL(NMX2),THALL(NMX2),ZLALL(NMX2),TLALL(NMX2),
     3   THIGHS0(NMX2),HHIGHS0(NMX2),IDX0(NMX2)
      CHARACTER FNAME*80
      DATA LU,LU2/10,20/

      WRITE(*,*) 'PROGRAM EXTREMES.F'
      NMAX = NMAX1
      IF(NMAX .GT. NMX) THEN
        WRITE(*,*) 'Data number exceeds dimension of array',nmax,nmx
        WRITE(*,*) 'Reset array dimension nmx in extremes.'
        STOP
      ENDIF
      DO I = 1, NMAX
        T(I) = TTMP(I)
        H(I) = HORI(I)
      ENDDO
      DELHR = DELHR1
      DELAMP = DELAMP1
      DELPCT = DELPCT1
      IOPTA = IOPTA1
C      WRITE(*,*) 'delhr,delamp,delpct=',delhr,delamp,delpct,iopta
      IF(NMAX .EQ. 0) STOP
      NPERHR = INT(1.0/DELTT+0.001)
C      WRITE(*,*) ' nperhr= ',nperhr,'  delt= ',deltt

      HMIN =  1.0E+5
      HMAX = -1.0E+5
      DO N = 1, NMAX
        HMIN = AMIN1(H(N),HMIN)
        HMAX = AMAX1(H(N),HMAX)
      ENDDO
      WRITE(*,"(' Time series read: hmin, hmax=',2f10.5)") HMIN,HMAX

      IF(NTYPE .EQ. 0) THEN
        CALL HILO3_EVENTS(IPR1)
      ELSE IF(NTYPE .EQ. 1) THEN
        CALL HILO3(IPR1)
      ELSE
        WRITE(*,*) 'NTYPE is not 0 or 1, you have to reset.'
        STOP
      ENDIF  

      NUM_H = 0
      NUM_L = 0
      DO N = 1, NSMAX
        IF(IDX(N) .EQ. 1) THEN
          NUM_H = NUM_H+1
          THALL(NUM_H) = THIGHS(N)
          ZHALL(NUM_H) = HHIGHS(N)
        ELSE IF(IDX(N) .EQ. -1) THEN
          NUM_L = NUM_L+1
          TLALL(NUM_L) = THIGHS(N)
          ZLALL(NUM_L) = HHIGHS(N)
        ENDIF
      ENDDO
        
      NNN = 0
      NUM_H = 0
      NUM_L = 0
      DO N = 1, NSMAX
        NM = 0
        DO I = 1, NMAX
          IF(ABS(THIGHS(N)-T(I)) .LT. 0.01) THEN
            NM = I
            GOTO 566
          ENDIF
        ENDDO
566     CONTINUE
        IF(NM .GT. 0) THEN
          IF(IDX(N).EQ. 1) THEN
            HLOCALMAX = -1.0E+5
            DO NX = MAX0(1,NM-NPERHR), MIN0(NMAX,NM+NPERHR)
              IF(HORI(NX) .GT. HLOCALMAX) THEN
                HLOCALMAX = HORI(NX)
                TIMEMAX = T(NX)
              ENDIF
            ENDDO
            NUM_H = NUM_H+1
            NNN = NNN+1
            THALL(NUM_H) = TIMEMAX
            ZHALL(NUM_H) = HLOCALMAX
            HHIGHS(NNN) = HLOCALMAX
            THIGHS(NNN) = TIMEMAX
            IDX0(NNN) = IDX(N)
          ELSE IF(IDX(N) .EQ. -1) THEN
            HLOCALMIN = 1.0E+5
            DO NX = MAX0(1,NM-NPERHR), MIN0(NMAX,NM+NPERHR)
              IF(HORI(NX) .LT. HLOCALMIN) THEN
                HLOCALMIN = HORI(NX)
                TIMEMIN = T(NX)
              ENDIF
            ENDDO
            NNN = NNN+1
            NUM_L = NUM_L+1
            TLALL(NUM_L) = TIMEMIN
            ZLALL(NUM_L) = HLOCALMIN
            HHIGHS(NNN) = HLOCALMIN
            THIGHS(NNN) = TIMEMIN
            IDX0(NNN) = IDX(N)
          ENDIF
        ELSE
          WRITE(*,*) 'not find match time of t(i) with thighs(n)'
          WRITE(*,*) 'NM=', NM
        ENDIF 
      ENDDO
      NSMAX = NNN
999   CONTINUE

      NSMAX0 = NSMAX
      DO N = 1, NSMAX
        THIGHS0(N) = THIGHS(N)
        HHIGHS0(N) = HHIGHS(N)
        IDX0(N) = IDX(N)
      ENDDO

      WRITE(*,*) 'num_h,num_l = ', NUM_H,NUM_L
      WRITE(*,*) 'END EXTREMES PROGRAM'

      RETURN
      END


CC    use SVD to get refined peaks around the point m
CC    indx = index for hi (1) or low (-1)
      subroutine getpeak(ipr,m,indx,nwin,xmax,ymax,mpk,kpr)
      parameter (ndata=66,ma=6,mp=ndata,np=ma)
      parameter (nmx=735*240,nmx2=nmx/10)
      COMMON/B1/H(NMX),T(NMX),NMAX,NPERHR,
     1 THIGHS(NMX2),HHIGHS(NMX2),IHH(NMX2),ILL(NMX2),IDX(NMX2),
     2 NSMAX,DELHR,DELAMP,DELPCT,IOPTA,HMAX,HMIN
      REAL X(NDATA),Y(NDATA),SIG(NDATA),A(MA),U(MP,NP),V(NP,NP),
     1     W(NP),YP(NDATA),AFUNC(MA)

      IF(IPR .EQ. 1) WRITE(*,"(/,' <getpeak>',/,10x,'m=',i8)") M
CC  Initial window
      NWIN1 = NWIN
      IF(IPR .EQ. 1) WRITE(*,"('   window=',i4)") NWIN1

CC  Initialize
100   CONTINUE
      N1 = MAX(1,   M-NWIN1)
      N2 = MIN(NMAX,M+NWIN1)
      NTOT = N2-N1+1
      T0 = T(N1)

CC  Store time and initial height
      NN = 0
      NM = 0
      DO N = N1, N2
        NN = NN+1
        IF(N .EQ. M) NM = NN
        X(NN) = T(N)-T0
        Y(NN) = H(N)
        N3 = NN-NWIN1-1
        SIG(NN) = 1.0
      ENDDO

CC  Get SVD fitted curve
      CALL SVDFIT(X,Y,SIG,NTOT,A,MA,U,V,W,MP,NP,CHISQ)

CC  Get the predictions
      DO N = 1, NTOT
        CALL FUNCS(X(N),AFUNC,MA)
        SUM = 0.0
        DO J = 1, MA
          SUM = SUM+A(J)*AFUNC(J)
        ENDDO    
        YP(N) = SUM
        N3 = N-NWIN1-1
        IF(IPR .EQ. 1) WRITE(*,"(' n=',i4,' win=',i4,' t=',
     1    f12.5,' yp=',f9.5)") N,N3,X(N)+T0,YP(N)
      ENDDO

CC  Find max or min
      METH = 2
      IF(IPR.EQ.1) WRITE(*,"(/,'  find max or min: nm=',i5,
     1  ' meth=',i2,' indx=',i2)") NM,METH,INDX
      IF(METH .EQ. 2) GOTO 120

CC  Method 1: look for extrema within a fixed time domain
      YMAX = -1.0E+5
      YMIN =  1.0E+5
      NDEL = NPERHR
      NPK = 0
      MPK = 0

      DO N = MAX0(1,NM-NDEL),MIN0(NTOT,NM+NDEL)
        IF(INDX .EQ. 1 .AND. N .GT. 0 .AND. N .LE. NTOT) THEN
          IF(YP(N) .GE. YP(N-1) .AND. YP(N) .GE. YP(N+1)) THEN
            NPK = N
            GOTO 125
          ENDIF
        ELSE IF(INDX .EQ. -1 .AND. N .GT. 0 .AND. N .LE. NTOT) THEN
          IF(YP(N) .LE. YP(N-1) .AND. YP(N) .LE. YP(N+1)) THEN
            NPK = N
            GOTO 125
          ENDIF
        ENDIF
      ENDDO
120   CONTINUE

      NDEL = NPERHR
      NPK = 0
      MPK = 0
      DO L = 0, NDEL, 1
        DO N = NM-L, NM+L, MAX0(1,2*L)
          IF(IPR .EQ. 1) WRITE(*,"(' l=',I2,' n=',I3,' y=',3F9.5)") 
     1      L,N,YP(N-1),YP(N),YP(N+1)
C  Zheng
C          IF(INDX .EQ. 1 .AND. N .GT. 0 .AND. N .LE. NTOT) THEN
          IF(INDX .EQ. 1 .AND. N .GT. 1 .AND. N .LE. NTOT) THEN
            IF(YP(N).GE.YP(N-1).AND.YP(N).GE.YP(N+1)) THEN
              NPK = N
              GOTO 125
            ENDIF
C  Zheng
C          ELSE IF(INDX .EQ. -1.AND.N .GT. 0.AND.N .LE. NTOT) THEN
          ELSE IF(INDX .EQ. -1 .AND. N .GT. 1 .AND. N .LE. NTOT) THEN
            IF(YP(N).LE.YP(N-1).AND.YP(N).LE.YP(N+1)) THEN
              NPK = N
              GOTO 125
            ENDIF
          ENDIF
        ENDDO
      ENDDO
125   CONTINUE

CC  Save here
      IF(NPK .GT. 0) THEN
        YMAX = YP(NPK) 
        XMAX = X(NPK)+T0
        MPK = N1+NPK-1
      ENDIF

CC  Check for convergence
      IF(NPK .EQ. 0) THEN
        IF(NWIN1 .GT. NPERHR) THEN
          NWIN1 = NWIN1/2
          WRITE(*,"(4X,'NON-CONVERGENT. REDUCE NWIN TO ',I3)") NWIN1
          GOTO 100
        ELSE
          GOTO 200
        ENDIF
      ENDIF
200   CONTINUE

CC  Write results
      IF(IPR .EQ. 1) WRITE(*,"('  xmax=',f12.5,' ymax=',f8.4,
     1 ' at npk=',i6,' mpk=',i6)") XMAX,YMAX,NPK,MPK
      IF(NPK .EQ. 0) THEN
        WRITE(*,*)'NPK=0, NO EXTREMA WAS NOT FOUND AROUND TIME=',T(M)
      ENDIF

      RETURN
      END


CC  Find peaks in 1/2-hr avgs, put SVD curve thru,
CC  Get new peaks, eliminate if too close in time or elev
      SUBROUTINE HILO3(IPR)
      PARAMETER (NMX=735*240,NMX2=NMX/10)
      COMMON/B1/H(NMX),T(NMX),NMAX,NPERHR,
     1 THIGHS(NMX2),HHIGHS(NMX2),IHH(NMX2),ILL(NMX2),IDX(NMX2),
     2 NSMAX,DELHR,DELAMP,DELPCT,IOPTA,HMAX,HMIN
      DIMENSION ZB(NMX),IB(NMX),TB(NMX),NAP(NMX),JDX(NMX)
      CHARACTER TYPE*16,SGN*1

      WRITE(*,*) 'PROGRAM HILO3'
C      WRITE(*,*) 'delhr,delamp,delpct = ',delhr,delamp,delpct,iopta

CC  Get half-hour averages
      WRITE(*,*) 'Get 30 minute averages'
      NA = 0  ! the number of the set of points being averaged
      INC = MAX(1,NPERHR/2)
      DO 110 N = 1, NMAX-INC, INC
        NA = NA+1
        ZB(NA) = 0.0
        TB(NA) = 0.0
        DO I = 1, INC
          ZB(NA) = ZB(NA)+H(N+I-1)
          TB(NA) = TB(NA)+T(N+I-1)
        ENDDO
        ZB(NA) = ZB(NA)/FLOAT(INC)
        TB(NA) = TB(NA)/FLOAT(INC)
110   CONTINUE
      NAMAX = NA

CC  Make sure there are no duplicate peaks
      DO N = 2, NAMAX
        IF(ZB(N) .EQ. ZB(N-1)) ZB(N-1) = ZB(N-1)+.001
      ENDDO

CC  <hilo3> Look for peaks in half-hour averages
      WRITE(*,*) 'Positive and negative peaks in 30 min averages'
      NA = 1
      IB(NA) = 0
      IF(IPR .EQ. 1) WRITE(*,"(' na=',i6,' tb=',f11.3,' zb=',f9.4,
     1  ' ib=',i3)") NA,TB(NA),ZB(NA),IB(NA)

      NA = NAMAX
      IB(NA) = 0
CC  Loop thru avergaes, tag as high (1) or low (-1)
      DO NA = 2, NAMAX-1
        IB(NA) = 0
        IF(ZB(NA) .GE. ZB(NA-1) .AND. ZB(NA) .GE. ZB(NA+1)) IB(NA) = 1
        IF(ZB(NA) .LE. ZB(NA-1) .AND. ZB(NA) .LE. ZB(NA+1)) IB(NA) = -1
        IF(IPR .EQ. 1) WRITE(*,"(' na=',i6,' tb=',f11.3,' zb=',f9.4,
     1    ' ib=',i3)") NA,TB(NA),ZB(NA),IB(NA)
      enddo
      NA = NAMAX
      IB(NA) = 0
      IF(IPR .EQ. 1) WRITE(*,"(' na=',i6,' tb=',f11.3,' zb=',f9.4,
     *  ' ib=',i3)") NA,TB(NA),ZB(NA),IB(NA)

CC  <hilo3> Collect only peaks of avgs, save to file for plotting
      NP = 0
      DO NA = 2, NAMAX-1
        IF(IB(NA) .NE. 0) THEN
          NP = NP+1
          NAP(NP) = NA
        ENDIF
      ENDDO
      NPMAX = NP
      WRITE(*,*) 'First peak number = ', NPMAX

CC  <hilo3> Find revised peaks by SVD
      WRITE(*,*) 'Find revised peaks using SVD'
      NP = 0
      NS = 0
      ZMAX = -1.0E5
      ZMIN =  1.0E5
      DO NX = 1, NAMAX
        ZMAX = AMAX1(ZMAX,ZB(NX))
        ZMIN = AMIN1(ZMIN,ZB(NX))
      ENDDO
      DELZMAX = (ZMAX-ZMIN)/2.0
      WRITE(*,"(' zmax,min,diff = ',3f7.3)") ZMAX,ZMIN,ZMAX-ZMIN

      DO NP = 1, NPMAX
        NA = NAP(NP)
        M = INC*(NA-1)+INC/2+1
        INDX = IB(NA)
        JPR = 0           
CC  Find max and min in 3.0 hr interval, set window
        ZMAX = -1.0E5
        ZMIN =  1.0E5
        DO NX = MAX0(1,NA-3), MIN0(NAMAX,NA+3)
          ZMAX = AMAX1(ZMAX,ZB(NX))
          ZMIN = AMIN1(ZMIN,ZB(NX))
        ENDDO
        NWIN = 2*NPERHR
        IF(ZMAX-ZMIN .GE. DELZMAX) NWIN = NWIN/2
        XMAX = 0.0
        YMAX = 0.0
        KPR = 0
        CALL GETPEAK(JPR,M,INDX,NWIN,XMAX,YMAX,MPK,KPR)

CC  Save result
        IF(MPK .GT. 0) THEN
          NS = NS+1
          THIGHS(NS) = XMAX
          HHIGHS(NS) = YMAX
          IDX(NS) = INDX
        ENDIF
      ENDDO
      NSMAX = NS
      WRITE(*,*) 'Finish getpeak'

CC  <hilo3> Sort peaks into ascending order (by time)
130   CONTINUE
      NSW = 0
      DO N = 1, NSMAX-1
        IF(THIGHS(N) .GT. THIGHS(N+1)) THEN
          T1 = THIGHS(N)
          H1 = HHIGHS(N)
          I1 = IDX(N)
          THIGHS(N) = THIGHS(N+1)
          THIGHS(N) = THIGHS(N+1)
          IDX(N) = IDX(N+1)
          HHIGHS(N+1) = H1
          THIGHS(N+1) = T1
          IDX(N+1) = I1
          NSW = NSW+1
        ENDIF
      ENDDO
      IF(NSW .GT. 0) GOTO 130

CC  If two consecutive highs or two lows, remove one 
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DO N = 2, NSMAX
        IF(IDX(N) .EQ. IDX(N-1)) THEN
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N-1 
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N-1
          IF(HHIGHS(N) .EQ. HHIGHS(N-1)) NN = N-1
          JDX(NN) = 0
        ENDIF
      ENDDO

CC  Remove doubles
      N = 0
      DO NN = 1, NSMAX
        IF(JDX(NN) .EQ. 1) THEN
          N = N+1
          THIGHS(N) = THIGHS(NN)
          HHIGHS(N) = HHIGHS(NN)
          IDX(N) = IDX(NN)
        ENDIF
      ENDDO
      NSMAX = N

CC  Reset amplitude criterion
C      WRITE(*,"(//,' set delamp. iopta=',i2)") IOPTA
      IF(IOPTA .EQ. 1) THEN
        WRITE(*,"(' delamp=',f6.4)") DELAMP
      ELSE IF(IOPTA .EQ. 2) THEN
        DELAMP = DELPCT*(HMAX-HMIN)
        WRITE(*,"(' hmin,hmax=',2f10.5,' delpct=',f6.4)") 
     1    HMIN,HMAX,DELPCT
        WRITE(*,"(' delamp=',f6.4)") DELAMP
      ELSE IF(IOPTA .EQ. 3) THEN
        HMAX = 0.0
        DO N = 2, NSMAX-1
          HMAX = HMAX+ABS(HHIGHS(N)-.5*(HHIGHS(N-1)+HHIGHS(N+1)))
        ENDDO
        HMAX = HMAX/(NSMAX-2)
        DELAMP = DELPCT*HMAX
        WRITE(*,"(' hmax=',f10.5,' delpct=',f7.4,' delamp=',f7.4)") 
     1     HMAX,DELPCT,DELAMP
      ENDIF

CC  Check for extrema too close in time or amplitude
C      WRITE(*,"(//,' extrema check: delhr=',f5.3,' delamp=',f8.5)")
C     1  DELHR,DELAMP

CC  Loop thru the highs and lows
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DELTT = 0.0
      DO 150 N = 2, NSMAX
        DELTT = THIGHS(N)-THIGHS(N-1)

CC  Check on times
        IF(DELTT .LT. DELHR) THEN
C          WRITE(*,*) 'deltt = ', DELTT
          IF(JDX(N-1) .EQ. 0) GOTO 150
C          WRITE(*,"(5x,'time prob at n=',i5)") n
          JDX(N) = 0
          WRITE(*,"(1x,'remove peaks at n=',i5)") n
          GOTO 150
        ENDIF
150   CONTINUE        

CC  Recount after eliminating spurious peaks
      NS = 0
      DO N = 1, NSMAX
        IF(JDX(N) .NE. 0) THEN
          NS = NS+1
          THIGHS(NS) = THIGHS(N)
          HHIGHS(NS) = HHIGHS(N)
          IDX(NS) = IDX(N)
        ENDIF
      ENDDO
      NSMAX = NS

CC  If two consecutive highs or two lows, remove one 
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DO N = 2, NSMAX
        IF(IDX(N) .EQ. IDX(N-1)) THEN
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N-1
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N-1
          IF(HHIGHS(N) .EQ. HHIGHS(N-1)) NN = N-1
          JDX(NN) = 0
        ENDIF
      ENDDO

CC  Remove doubles
      N = 0
      DO NN = 1, NSMAX
        IF(JDX(NN) .EQ. 1) THEN
          N = N+1
          THIGHS(N) = THIGHS(NN)
          HHIGHS(N) = HHIGHS(NN)
          IDX(N) = IDX(NN)
        ENDIF
      ENDDO
      NSMAX = N

CC  Check on amplitude differences between highs and lows > delamp
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DELHH = 0.0
      DO 160 N = 2, NSMAX
        DELHH = ABS(HHIGHS(N)-HHIGHS(N-1))
        IF(DELHH .LT. DELAMP) THEN
          WRITE(*,"(1x,'n=',i5,'  h=',2f10.4,'  delh=',f10.4)")
     1          N,HHIGHS(N-1),HHIGHS(N),DELHH
          IF(IDX(N)+IDX(N-1) .EQ. 0) THEN
            JDX(N) = 0
            WRITE(*,"(1x,'Remove peaks at n=',i5,'  h=',2f10.4)")
     1          N,THIGHS(N),HHIGHS(N)
          ENDIF
          GOTO 160
        ENDIF
160   CONTINUE
      
CC  Recount after eliminating spurious peaks
      NS=0
      DO N = 1, NSMAX
        IF(JDX(N) .NE. 0) THEN
          NS = NS+1
          THIGHS(NS) = THIGHS(N)
          HHIGHS(NS) = HHIGHS(N)
          IDX(NS) = IDX(N)
        ENDIF
      ENDDO
      NSMAX = NS

CC  If two consecutive highs or two lows, remove one 
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DO N = 2, NSMAX
        IF(IDX(N) .EQ. IDX(N-1)) THEN
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N-1 
          IF(IDX(N) .EQ. 1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .GT. HHIGHS(N-1)) NN = N
          IF(IDX(N) .EQ. -1 .AND. HHIGHS(N) .LT. HHIGHS(N-1)) NN = N-1
          IF(HHIGHS(N) .EQ. HHIGHS(N-1)) NN = N-1
          JDX(NN) = 0
        ENDIF
      ENDDO

CC  Remove doubles
      N = 0
      DO NN = 1, NSMAX
        IF(JDX(NN) .EQ. 1) THEN
          N = N+1
          THIGHS(N) = THIGHS(NN)
          HHIGHS(N) = HHIGHS(NN)
          IDX(N) = IDX(NN)
        ENDIF
      ENDDO
      NSMAX = N

CC  Print
      WRITE(*,"(' Summary of peaks passing extrema test')")
      WRITE(*,"(' (but look for double highs or lows)')")
      IBAL = 0
      DO N = 2, NSMAX
        TYPE = ' '
        IF(IDX(N) .EQ. IDX(N-1)) THEN
          IBAL = IBAL+1
          TYPE = 'same as previous'
        ENDIF
      ENDDO
      IF(IBAL. GT. 0) THEN
        WRITE(*,"(/,'  ***',i5,' double highs or lows ***')") IBAL
      ENDIF

      RETURN
      END


CC  Find peaks in 1/2-hr avgs, put SVD curve thru.
CC  Get new peaks, eliminate if too close in time or elev.
      SUBROUTINE HILO3_EVENTS(IPR)
      PARAMETER (NMX=735*240,NMX2=NMX/10)
      COMMON/B1/H(NMX),T(NMX),NMAX,NPERHR,
     1 THIGHS(NMX2),HHIGHS(NMX2),IHH(NMX2),ILL(NMX2),IDX(NMX2),
     2 NSMAX,DELHR,DELAMP,DELPCT,IOPTA,HMAX,HMIN
      CHARACTER TYPE*16,SGN*1
      DIMENSION TH1(NMX2),HH1(NMX2),IDXH(NMX2),
     1          TL1(NMX2),HL1(NMX2),IDXL(NMX2)
      DIMENSION ZB(NMX),TB(NMX),IB(NMX),NAP(NMX),JDX(NMX)

      WRITE(*,*) 'PROGRAM HILO3_EVEBTS'
C      WRITE(*,*) 'delhr,delamp,delpct=',DELHR,DELAMP,DELPCT,IOPTA

CC  Get half-hour averages
      WRITE(*,*) ' Get 30 minute averages'
      NA = 0  ! the number of the set of points being averaged
      INC = MAX(1,NPERHR/2)
      DO 110 N = 1, NMAX-INC, INC
        NA = NA+1
        ZB(NA) = 0.0
        TB(NA) = 0.0
        DO I = 1, INC
          ZB(NA) = ZB(NA)+H(N+I-1)
          TB(NA) = TB(NA)+T(N+I-1)
        ENDDO
        ZB(NA) = ZB(NA)/FLOAT(INC)
        TB(NA) = TB(NA)/FLOAT(INC)
110   CONTINUE
      NAMAX = NA

CC  Make sure there are no duplicate peaks
      DO N = 2, NAMAX
        IF(ZB(N) .EQ. ZB(N-1)) ZB(N-1) = ZB(N-1) + 0.001
      ENDDO

CC  <hilo3_events> look for peaks in half-hour averages
      WRITE(*,*) 'Positive and negative peaks in 30 min averages'
      NA = 1
      IB(NA) = 0
      IF(IPR .EQ. 1) WRITE(*,"(' NA=', I6,' TB=', F11.3,
     1  ' ZB=', F9.4,' IB=', I3)") NA,TB(NA),ZB(NA),IB(NA)

CC  Loop thru avergaes, tag as high (1) or low (-1)
      DO NA = 2, NAMAX-1
        IB(NA) = 0
        IF(ZB(NA) .GE. ZB(NA-1) .AND. ZB(NA) .GE. ZB(NA+1)) IB(NA) =  1
        IF(ZB(NA) .LE. ZB(NA-1) .AND. ZB(NA) .LE. ZB(NA+1)) IB(NA) = -1
        IF(IPR .EQ. 1) WRITE(*,"(' NA=', I6,' TB=', F11.3,
     1    ' ZB=', F9.4,' IB=', I3)") NA,TB(NA),ZB(NA),IB(NA)
      ENDDO
      NA = NAMAX
      IB(NA) = 0
      IF(IPR .EQ. 1) WRITE(*,"(' NA=', I6,' TB=', F11.3,' ZB=', F9.4,
     1  ' IB=', I3)") NA,TB(NA),ZB(NA),IB(NA)

CC  <hilo3_events> collect only peaks of avgs, save to file for plotting
      NP = 0
      DO NA = 2, NAMAX-1
        IF(IB(NA) .NE. 0) THEN
          NP = NP+1
          NAP(NP) = NA
        ENDIF
      ENDDO
      NPMAX = NP
      WRITE(*,*) 'First peak number=', NPMAX

CC  <hilo3_events> find revised peaks by SVD
      WRITE(*,*) 'Find revised peaks using SVD'
      NP = 0
      NS = 0
      ZMAX = -1.0E5
      ZMIN =  1.0E5
      DO NX = 1, NAMAX
        ZMAX = AMAX1(ZMAX,ZB(NX))
        ZMIN = AMIN1(ZMIN,ZB(NX))
      ENDDO
      DELZMAX = (ZMAX-ZMIN)/2.0
      WRITE(*,"(' ZMAX,MIN,DIFF=',3F7.3)") ZMAX,ZMIN,ZMAX-ZMIN

      DO NP = 1, NPMAX
        NA = NAP(NP)
        M = INC*(NA-1)+INC/2+1  ! mid-point in original time series t(n)
        INDX = IB(NA)
        JPR=0           

CC  Find max and min in 3.0 hr interval, set window
        ZMAX=-1.0E5
        ZMIN= 1.0E5
        DO NX = MAX0(1,NA-3), MIN0(NAMAX,NA+3)
          ZMAX = AMAX1(ZMAX,ZB(NX))
          ZMIN = AMIN1(ZMIN,ZB(NX))
        ENDDO
        NWIN = 2*NPERHR
        IF(ZMAX-ZMIN .GE. DELZMAX) NWIN = NWIN/2

CC  Put SVD curve thru
        XMAX = 0.0
        YMAX = 0.0
        KPR = 0
        CALL GETPEAK(JPR,M,INDX,NWIN,XMAX,YMAX,MPK,KPR)

CC  Save result
        IF(MPK .GT. 0) THEN
          NS = NS+1
          THIGHS(NS) = XMAX
          HHIGHS(NS) = YMAX
          IDX(NS) = INDX
        ENDIF
      ENDDO
      NSMAX = NS
      WRITE(*,*) 'Finish getpeak'

CC  <hilo3_events> put peaks into ascending order (by time)
130   CONTINUE
      NSW = 0
      DO N = 1, NSMAX-1
        IF(THIGHS(N) .GT. THIGHS(N+1)) THEN
          T1 = THIGHS(N)
          H1 = HHIGHS(N)
          I1 = IDX(N)
          THIGHS(N) = THIGHS(N+1)
          THIGHS(N) = THIGHS(N+1)
          IDX(N) = IDX(N+1)
          HHIGHS(N+1) = H1
          THIGHS(N+1) = T1
          IDX(N+1) = I1
          NSW = NSW+1
        ENDIF
      ENDDO
      IF(NSW .GT. 0) GOTO 130

CC  Reset amplitude criterion
C      WRITE(*,"(//,' set delamp. iopta = ',i2)") IOPTA
      IF(IOPTA .EQ. 1) THEN
        WRITE(*,"('   delamp = ',f6.4)") DELAMP
      ELSE IF(IOPTA .EQ. 2) THEN
        DELAMP=DELPCT*(HMAX-HMIN)
        WRITE(*,"('   hmin,hmax = ',2f10.5,' delpct = ',f6.4)") 
     *    HMIN,HMAX,DELPCT
        WRITE(*,"('   delamp = ',f6.4)") DELAMP
      ELSE IF(IOPTA .EQ. 3) THEN
        HMAX = 0.0
        DO N = 2, NSMAX-1
          HMAX = HMAX+ABS(HHIGHS(N)-0.5*(HHIGHS(N-1)+HHIGHS(N+1)))
        ENDDO
        HMAX = HMAX/(NSMAX-2)
        DELAMP = DELPCT*HMAX
        WRITE(*,"(' hmax =',f10.5,' delpct =',f7.4,' delamp =',f7.4)")
     1    HMAX,DELPCT,DELAMP
      ENDIF

CC  Loop thru the highs and lows
      DO N = 1,NSMAX
        JDX(N) = 1
      ENDDO
      DELTT = 0.0
      NS = 1
      DO 150 N = 2, NSMAX
        DELTT = THIGHS(N)-THIGHS(NS)
CC  Check on times
        IF(DELTT .GT. DELHR) THEN
          NS = NS+1
          THIGHS(NS) = THIGHS(N)
          HHIGHS(NS) = HHIGHS(N)
          IDX(NS) = IDX(N)
        ELSE IF(DELTT .LE. DELHR) THEN
          IF(((IDX(NS)+IDX(N)) .GT. 1) .AND.
     1        (HHIGHS(N) .GT. HHIGHS(NS))) THEN
            THIGHS(NS) = THIGHS(N)
            HHIGHS(NS) = HHIGHS(N)
            IDX(NS) = IDX(N)
          ELSE IF(((IDX(NS)+IDX(N)) .LT. -1) .AND.
     1        (HHIGHS(N) .LT. HHIGHS(NS))) THEN
            THIGHS(NS) = THIGHS(N)
            HHIGHS(NS) = HHIGHS(N)
            IDX(NS) = IDX(N)
          ELSE IF(((IDX(NS)+IDX(N)) .EQ. 0) .AND.
     1        (ABS(HHIGHS(N)-HHIGHS(NS)) .GT. DELAMP) .AND.
     1        (DELTT .GT. DELHR/3)) THEN    !! add on 5/5/2006 zaj
            NS = NS+1
            THIGHS(NS) = THIGHS(N)
            HHIGHS(NS) = HHIGHS(N)
            IDX(NS) = IDX(N)
          ENDIF    
        ENDIF
150   CONTINUE        
      NSMAX = NS

CC  Check on amplitude differences between highs and lows > delamp
      DO N = 1, NSMAX
        JDX(N) = 1
      ENDDO
      DELHH = 0.0
      DO 160 N = 2, NSMAX
        DELHH = ABS(HHIGHS(N)-HHIGHS(N-1))
        IF(DELHH .LT. DELAMP) THEN
          IF(IDX(N)+IDX(N-1) .EQ. 0) THEN
            JDX(N) = 0
            WRITE(*,"(5x,'Remove peaks at n=',i5,'h=',3f10.4)")
     1            N,THIGHS(N),HHIGHS(N),DELHH
          ENDIF
          GOTO 160
        ENDIF
 160  CONTINUE

CC  Recount after eliminating spurious peaks
      NS = 0
      DO N = 1, NSMAX
        IF(JDX(N) .NE. 0) THEN
          NS = NS+1
          THIGHS(NS) = THIGHS(N)
          HHIGHS(NS) = HHIGHS(N)
          IDX(NS) = IDX(N)
        ENDIF
      ENDDO
      NSMAX = NS

      WRITE(*,"(' Summary of peaks passing extrema test')")
      WRITE(*,"(' (but look for double highs or lows)')")
      IBAL = 0     ! assume no double highs or lows
      DO N=2,NSMAX
        TYPE=' '
        IF(IDX(N) .EQ. IDX(N-1)) THEN
          IBAL=IBAL+1
          TYPE='Same as previous'
        ENDIF
      ENDDO
      IF(IBAL .GT. 0) THEN
        WRITE(*,"(/,'  ***',i5,' double highs or lows ***')") IBAL
      ENDIF

      RETURN
      END


CC  Add factor parameter: 
CC  upper limit=mean + factor * SD
CC  lower limit=mean - factor * SD
CC  if factor =< 0 then using fixed value instead of sigma rule 
CC  upper limit= -factor      
CC  lower limit= factor
      SUBROUTINE EVENTS(XTMP,YTMP,NUMB,THIGHS0,HHIGHS0,IDX0,
     1   NSMAX0,HUPPER,HLOWER,FACTOR)
      DIMENSION XTMP(NUMB),YTMP(NUMB),XTMP1(NUMB),YTMP1(NUMB),
     1   THIGHS0(NSMAX0),HHIGHS0(NSMAX0),IDX0(NSMAX0),
     2   THIGHS1(NSMAX0),HHIGHS1(NSMAX0),IDX1(NSMAX0)
      HWIN = 85   !! means using 7 days window (170 hours) to calculate SD
      N2 = 0
      DO J = 1, NSMAX0
        N0 = 0
        DO I = 1, NUMB
          IF((XTMP(I) .GE. THIGHS0(J)-HWIN) .AND.
     1       (XTMP(I) .LT. THIGHS0(J)+HWIN)) THEN
            N0 = N0+1    
            XTMP1(N0) = XTMP(I)
            YTMP1(N0) = YTMP(I)
          ENDIF
        ENDDO
        IF(N0 .GT. 1) THEN
          IF(FACTOR .GT. 0.0) THEN
            SM0 = SM(YTMP1,NUMB,N0)
            SD0 = SD(YTMP1,NUMB,N0)
            HUPPER = SM0+FACTOR*SD0
            HLOWER = SM0-FACTOR*SD0
          ELSE  !! FACTOR < 0 FOR SPECIFYING A CONSTANT UPPER AND LOWER LIMIT
            SM0 = SM(YTMP1,NUMB,N0)
            HUPPER = SM0-FACTOR
            HLOWER = SM0+FACTOR
          ENDIF   
          IF(SD0 .LT. 0.03) GOTO 50  !! IGNORE THE SMALLER PEAKS 05/04/2006
        ENDIF

        IF((HHIGHS0(J) .LE. HLOWER) .AND. (IDX0(J) .LT. 0)) THEN
          N2 = N2+1
          THIGHS1(N2) = THIGHS0(J)
          HHIGHS1(N2) = HHIGHS0(J)
          IDX1(N2) = IDX0(J)
          WRITE(78,100)J,THIGHS0(J)/24.0,HHIGHS0(J),IDX0(J),
     1        HUPPER,HLOWER,SD0,N0,N2
        ELSE IF((HHIGHS0(J) .GE. HUPPER).AND.(IDX0(J) .GT. 0)) THEN 
          N2 = N2+1
          THIGHS1(N2) = THIGHS0(J)
          HHIGHS1(N2) = HHIGHS0(J)
          IDX1(N2) = IDX0(J)
          WRITE(78,100) J,THIGHS0(J)/24.0,HHIGHS0(J),IDX0(J),
     1        HUPPER,HLOWER,SD0,N0,N2
        ENDIF
50      CONTINUE
      ENDDO
      WRITE(*,*)  'NMAX,N2=', NSMAX0,N2

      NSMAX0 = N2
      DO I = 1, N2
        THIGHS0(I) = THIGHS1(I)
        HHIGHS0(I) = HHIGHS1(I)
        IDX0(I) = IDX1(I)
      ENDDO
100   FORMAT(I7,2F10.4,I3,3F10.4,2I7) 

      RETURN
      END

