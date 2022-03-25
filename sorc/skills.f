c      skills.f. read data and create NOS skill scores of water level
c                (All times in hrs)
C      skills.x < skills.in
C f90 skills.f foufil.f extremes.f svd.f equal_interval.f slack.f table.f -o skills.x
c      variables -
c         tstart = start time at beginning of data used for SA in days 
c         tfinish = end time of used data for skill assessment in days
c         delt = time interval of time series data (hr)
c         jmax = number of forecasts per day (includes 00), =4 for 4 times
C                daily forecasts at cycle 00Z, 06Z, 12Z and 18Z
c         tcut = cutoff period (hrs) Fourier filtering, 
c              =30 is for 30-hour low-pass filtering
C         Igapfill is control switch of gap filling with interpolation 
C             0: filling with missing value -999.0; 
C             1: filling with interpolation value
C         method is index of interpolation method
C             0: Cubic spline  
C             1: Singular Value Decomposition(SVD); 
C         criteria1 (in hours) uses linear interpolation if gap < criteria1    
C         criteria2 (in hours) uses spline or SVD interpolation method 
C               when criteria1 < gap < criteria2 
C               fill gaps using missing value -999.0 while gap > criteria2    
c         a(i) = astronomical tidal prediction series
c         o(i) = observation series
c         of(i) = observation series after filtering
c         s(i) = tide-only model simulation series
c         er(i) = random error series
c         t() = time of astro tide, obs, model tides (hrs)
C         hi(i)= model hindcast
c         c(i) = nowcast
c         tc(i) = time for nowcasts (hrs)
c         f(L,J,K) = forecast
c         tf(L,J,K) = time for forecasts (f and pf) (hrs)
c         pf(L,J,K) = persistance forecast (same time as fcst, tf)
C         L = index in cycle including 00z and 24Z forecasts.
C            L=241 for 24-hour forecasts with time interval=6 minutes.
C         J= number of forecast cycle, J=1:00Z; J=2:06Z; J=3:12Z; J=4:18Z 
C         K= number of forecast days. K=1:first day; K=2:second day; 
C            K=3:third day, etc.
c         p() = predicted series
c         r() = reference series
c         e() = error     series
c         tp() = time for predicted, ref, and error series
C         ISWITCH(6) = scenario switch. = 0 is off; =1, is on
C*******************************************************************
C           SCENARIO 1:  TIDAL SIMULATION ONLY  
C           SCENARIO 2   SCENARIO: HINDCAST                 
C           SCENARIO 3   SCENARIO: NOWCAST                 
C           SCENARIO 4   SCENARIO: FORECAST                
C           SCENARIO 5   COMPARISON: PERSISTENCE FORECAST          
C           SCENARIO 6   COMPARISON: ASTRONOMICAL TIDE PREDICTION 
C*******************************************************************
C     Modified by Lianyuan Zheng pn 03/10/2017

      PROGRAM SKILLS
      include 'skills.inc'

      CHARACTER*200 FILESHORT,STATION_FILE,STAID,STATIONNAME,DATADIR
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      DIMENSION XTMP(IMX),YTMP(IMX),XNEW(IMX),YNEW(IMX)
      DIMENSION ZHALL(NMX),THALL(NMX),ZLALL(NMX),TLALL(NMX)
      DIMENSION IFRST(IMX),ILAST(IMX),HHWP0(NMX),
     1  THWP0(NMX),HLWP0(NMX),TLWP0(NMX),IDXR(NMX),IDX0(NMX),
     2  THIGHSR(NMX),HHIGHSR(NMX),THIGHS0(NMX),HHIGHS0(NMX)

      WRITE(*,"(' PROGRAM skills.f ')")
      READ(5,*) IYRS,MMS,IDDS,IHHS,IMINS 
      YEARB  = IYRS*1.0
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      YEARB  = IYRS*1.0
      MONTHB = MMS*1.0
      DAYB   = IDDS*1.0
      HOURB  = IHHS*1.0 + IMINS/60.0
      JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JBASE_DATE+1.0
      TSTART = JDAY

      READ(5,*) IYRE,MME,IDDE,IHHE,IMINE
      YEARB  = IYRE*1.0
      MONTHB = MME*1.0
      DAYB   = IDDE*1.0
      HOURB  = IHHE*1.0 + IMINE/60.0
      JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JBASE_DATE+1.0
      TFINISH = JDAY + 2  !!! TWO DAYS MORE FOR FORECAST SCENARIO
      IYEAR  = IYRS
      READ(5,*) NTYPE,DELT,DELT_O,DELT_T,DELT_M,JMAX,TCUT,FACTOR
      DELT   = DELT/60.0    !! 11/30/2006  CONVERT FROM MINUTES INTO HOURS
      DELT_O = DELT_O/60.0  !! 11/30/2006  CONVERT FROM MINUTES INTO HOURS
      DELT_T = DELT_T/60.0  !! 11/30/2006  CONVERT FROM MINUTES INTO HOURS
      DELT_M = DELT_M/60.0  !! 11/30/2006  CONVERT FROM MINUTES INTO HOURS
      READ(5,*) (ISWITCH(I),I = 1, 6),ISURGE  ! =0=off, 1=on
      READ(5,*) IGAPFILL,CRITERIA1,CRITERIA2,METHOD 

C  Read in parameters for subroutine extremes
      READ(5,*) IPR1
      READ(5,*) DELHR,DELAMP,DELPCT,IOPTA
      READ(5,*) X1,X2,X11
      X00 = X1
      READ(5,*) KINDAT
      READ(5,"(A200)") STATION_FILE
      READ(5,*) NFDURATION
      READ(5,*) FINCLUDE0
      READ(5,"(A200)") DATADIR
      IF(FINCLUDE0) THEN
        NCUT0  = INT(NFDURATION/DELT_M+0.1) + 1 
        IDUMMY = INT(NFDURATION/DELT+0.1) + 1
      ELSE
        NCUT0  = INT(NFDURATION/DELT_M+0.1)
        IDUMMY = INT(NFDURATION/DELT+0.1)
      END IF

      OPEN(1,FILE = TRIM(STATION_FILE))
      NDAY  = INT(TFINISH-TSTART)  ! number of days
      IPDAY = INT(24.0/DELT+0.1)   ! number of values in 24 hours
      IMAXA = INT(0.1+FLOAT(NDAY)*24.0/DELT)+1
      IF(IMAXA .GT. IMX) THEN
        WRITE(*,*) '************************************************'
        WRITE(*,*) 'IMX= ',IMX,'  IMAXA= ',IMAXA
        WRITE(*,*) 'Length of time series exceeds Max. Value (IMX)'
        WRITE(*,*) 'Reset KMX in skills.inc, and run COMPILE.sh '
        WRITE(*,*) '************************************************'
        STOP
      ENDIF 

C  Print input parameters
      WRITE(*,"(/,' TSTART =',F6.0, /,' TFINISH =',F6.0, /,
     1 ' DELT  =',F6.3, /,' JMAX   =',I6, /,' TCUT   =',F6.1)")
     2   TSTART,TFINISH,DELT,JMAX,TCUT
      WRITE(*,'(6I2)') ISWITCH
      WRITE(*,"(' IGAPFILL=',I2,' CRITERIA1=',F5.2,' CRITERIA2=',F5.2,
     1  ' METHOD=',I2)") IGAPFILL,CRITERIA1,CRITERIA2,METHOD 
      WRITE(*,*) 'IPR1=',IPR1
      WRITE(*,"(' DELHR,DELAMP   = ',2F8.3)") DELHR,DELAMP
      WRITE(*,"(' DELPCT,IOPTA   = ',F8.3,I4)") DELPCT,IOPTA
      WRITE(*,*) 'ERROR CRITERIA X1 =',X1,'  X2 = ',X2,'  X11 = ',X11
      WRITE(*,*) 'KINDAT = ',KINDAT
      WRITE(*,*) 'STATIONFILE = ',TRIM(STATION_FILE)
      WRITE(*,"(' NDAY = ',I3,'  IPDAY = ',I5)") NDAY,IPDAY
      WRITE(*,"(' IMAXA = ',I8)") IMAXA

C  Generate table for each station, loop through all stations      
10    READ(1,*,ERR=999,END=999) STAID,FILESHORT,STATIONNAME
      WRITE(*,*)
      WRITE(*,*) 'FILESHORT= ',TRIM(FILESHORT)
      WRITE(*,*) 'STATIONNAME= ',TRIM(STATIONNAME)
      READ(1,*) ALAT, ALON, DIRFLOOD, SDEPTH
      IMAXB = INT(0.1+FLOAT(NDAY)*24.0/DELT)+1

C  Read all data
      DO I = 1, IMX
        A(I)  = -999.9
        O(I)  = -999.9
        S(I)  = -999.9
        HI(I) = -999.9
        C(I)  = -999.9
        T(I)  = -999.9
        IF(ISURGE .GT. 0) THEN
          SGO(I) = -999.9
          SGM(I) = -999.9
        ENDIF
      ENDDO
      DO I = 1, MAXDURATION 
        DO J = 1, JMAX
          DO K = 1, KMX
            F(I,J,K)=-999.9
          ENDDO
        ENDDO
      ENDDO
         
      CALL READDATA(KINDAT,FILESHORT,DATADIR)
C  Convert times to hours
      DO I = 1, IDUMMY
        DO J = 1, JMAX
          DO K = 1, IMAXF
            TF(I,J,K)=TF(I,J,K)*24.0
          ENDDO
        ENDDO
      ENDDO
      DO I = 1, IMAXB
        T(I) = T(I)*24.0
        TC(I) = T(I)
      ENDDO
      NNN0 = 0
      DO I = 1, IMAXB
        OF(I) = O(I)
        IF(O(I) .GT. -900) NNN0 = NNN0+1
      ENDDO

C  Make Tabular Output
C      WRITE(*,*) 'imaxb after read = ',IMAXB,NNN0
      IF(NNN0 .GT. 1) THEN
        IF(NTYPE .EQ. 0) THEN
          CALL TABLE_NONTIDE(FILESHORT,STATIONNAME,KINDAT)
        ELSE IF(NTYPE .EQ. 1) THEN

C  Compute tidal ranges from tidal prediction time series
          IF(KINDAT .EQ. 2 .and. X00 .LT. 0.0) THEN
            NTMP = 0
            DO I = 1, IMAXB
              IF(ABS(A(I)) .LT. 900.0) THEN
                NTMP = NTMP+1
                XNEW(NTMP) = T(I)
                YNEW(NTMP) = A(I)
              ENDIF
            ENDDO

            GAP = 1.5*DELT
            CALL CONTINUOUS(XNEW,NTMP,GAP,NSEGMENTS,IFRST,ILAST)
            WRITE(*,*) 'Nsegments=', Nsegments
            NMAXRH9 = 0
            NMAXRL9 = 0
            DO NG = 1, NSEGMENTS
              ISTART = IFRST(NG)
              IEND = ILAST(NG)
              DIF = XNEW(IEND)-XNEW(ISTART)
              IF(DIF .GT. 48.0) THEN
                NUMB = IEND-ISTART+1  
                DO I = ISTART, IEND
                  I0 = I-ISTART+1 
                  XTMP(I0) = XNEW(I)
                  YTMP(I0) = YNEW(I)
                ENDDO
                CALL EXTREMES(XTMP,YTMP,NUMB,IPR1,DELHR,DELAMP,
     1           DELPCT,IOPTA,THIGHS0,HHIGHS0,IDX0,NSMAX0,TCUT,DELT,
     2           HHWP0,THWP0,NMAXRH0,HLWP0,TLWP0,NMAXRL0,NTYPE)
 
                DO J = 1, NMAXRH0
                  NMAXRH9 = NMAXRH9+1
                  HHWR(NMAXRH9) = HHWP0(J)
                  THWR(NMAXRH9) = THWP0(J)
                ENDDO

                DO J = 1, NMAXRL0
                  NMAXRL9 = NMAXRL9+1
                  HLWR(NMAXRL9) = HLWP0(J)
                  TLWR(NMAXRL9) = TLWP0(J)
                ENDDO
              ENDIF
            ENDDO

            AVG = 0.0
            DO IH = 1, NMAXRH9
              AVG = AVG+HHWR(IH)  
            ENDDO
            HAVG = AVG/FLOAT(NMAXRH9)
            AVG = 0.0
            DO IH = 1,NMAXRL9
              AVG=AVG+HLWR(IH)  
            ENDDO
            AVG = AVG/FLOAT(NMAXRL9)
            TIDAL_RANGE = HAVG-AVG

C  Redefine the error criteria X1 according to tidal range
            X1 = ABS(TIDAL_RANGE*X00)/100.0
            IF(X1 .LT. 0.15) X1 = 0.15
            WRITE(*,500) HAVG,AVG,X1,TRIM(FILESHORT)
          ENDIF

          CALL TABLE_TIDE(FILESHORT,STATIONNAME,KINDAT)
        ENDIF
      ELSE
        WRITE(*,"(/, 'No observation at ',a40)") TRIM(STATIONNAME)
      ENDIF
      GOTO 10
500   FORMAT('MHW, MLW, and New error criteria X1= ',3F8.2,' at',a10)
999   CONTINUE       

      STOP
      END


C  Print one line to table in a standard format
C  Variables:
C  Indx = index for line:  1=SM
C                          2=SM...MDPO, WOF
C                          3=SM...MDPO
C                          4=SM...POF
      SUBROUTINE PRTLINE1(TTL,TTL2,LU,B,TP,JJMX,IMAX,ERR0,INDX,GAP)
      include 'skills.inc'
      DIMENSION B(IMX),TP(IMX),XTMP(IMX),YTMP(IMX),
     1          IFRST(IMX),ILAST(IMX)
      CHARACTER TTL*9,TTL2*15

      SM0 = -999.0
      RMSE0 = -999.0
      SD0 = -999.0
      VOF0 = -999.0
      CF0 = -999.0
      POF0 = -999.0
      TMDON = -999.0
      TMDOP = -999.0
      IF(IMAX .GT. 1) THEN
        CALL CONTINUOUS(TP,IMAX,GAP,NSEGMENTS,IFRST,ILAST)
      END IF

      DO 150 K = 1, 2
C  Choose output channel
        IF(K .EQ. 1) LL = 6
        IF(K .EQ. 2) LL = LU

        IF(INDX .EQ. 1) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            WRITE(LL,130) TTL,TTL2,IMAX,SM0
          ENDIF

        ELSE IF(INDX .EQ. 2) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            RMSE0 = RMSE(B,IMX,IMAX)
            SD0 = SD(B,IMX,IMAX)
            VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
            CF0 = CF(B,IMX,IMAX,ERR0)
            POF0 = POF(B,IMX,IMAX,2.0*ERR0)

C  Using continuous time series to calculate tmdo
            TMDON = 0.0
            TMDOP = 0.0
            IF(IMAX .GT. 1) THEN
              DO NG = 1, NSEGMENTS
                ISTART = IFRST(NG)
                IEND = ILAST(NG)
                DIF = TP(IEND)-TP(ISTART)
                IF(DIF .GT. 48.0) THEN
                  NMAX = IEND-ISTART+1  
                  DO I = ISTART, IEND
                    I0 = I-ISTART+1 
                    XTMP(I0) = TP(I)
                    YTMP(I0) = B(I)
                  ENDDO
                  TMDO0 = TMDO(YTMP,XTMP,IMX,NMAX,-2.0*ERR0)
                  TMDO1 = TMDO(YTMP,XTMP,IMX,NMAX,2.0*ERR0)
                  IF(TMDO0 .GT. TMDON) TMDON = TMDO0
                  IF(TMDO1 .GT. TMDOP) TMDOP = TMDO1
                ENDIF
              ENDDO
            ENDIF
            WRITE(LL,130) TTL,TTL2,IMAX,SM0,RMSE0,SD0,VOF0,CF0,
     1         POF0,TMDON,TMDOP,WCOF,CORR_C
          ENDIF

        ELSE IF(INDX .EQ. 3) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            RMSE0 = RMSE(B,IMX,IMAX)
            SD0 = SD(B,IMX,IMAX)
            VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
            CF0 = CF(B,IMX,IMAX,ERR0)
            POF0 = POF(B,IMX,IMAX,2.0*ERR0)

C  Using continuous time series to calculate tmdo
            TMDON = 0.0
            TMDOP = 0.0
            IF(IMAX .GT. 1) THEN
              DO NG = 1,NSEGMENTS
                ISTART = IFRST(NG)
                IEND = ILAST(NG)
                DIF = TP(IEND)-TP(ISTART)
                IF(DIF .GT. 48.0) THEN
                  NMAX = IEND-ISTART+1  
                  DO I = ISTART, IEND
                    I0 = I-ISTART+1 
                    XTMP(I0) = TP(I)
                    YTMP(I0) = B(I)
                  ENDDO
                  TMDO0 = TMDO(YTMP,XTMP,IMX,NMAX,-2.0*ERR0)
                  TMDO1 = TMDO(YTMP,XTMP,IMX,NMAX,2.0*ERR0)
                  IF(TMDO0 .GT. TMDON) TMDON = TMDO0
                  IF(TMDO1 .GT. TMDOP) TMDOP = TMDO1
                ENDIF
              ENDDO
            ENDIF
            WRITE(LL,131) TTL,TTL2,IMAX,SM0,RMSE0,SD0,VOF0,CF0,
     1         POF0,TMDON,TMDOP,CORR_C
          ENDIF

        ELSE IF(INDX .EQ. 4) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            RMSE0 = RMSE(B,IMX,IMAX)
            SD0 = SD(B,IMX,IMAX)
            VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
            CF0 = CF(B,IMX,IMAX,ERR0)
            POF0 = POF(B,IMX,IMAX,2.0*ERR0)
            WRITE(LL,130) TTL,TTL2,IMAX,SM0,RMSE0,SD0,VOF0,CF0,POF0
          ENDIF

        ELSE IF(INDX .EQ. 5) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            RMSE0 = RMSE(B,IMX,IMAX)
            SD0 = SD(B,IMX,IMAX)
            VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
            CF0 = CF(B,IMX,IMAX,ERR0)
            POF0 = POF(B,IMX,IMAX,2.0*ERR0)

CC  Using continuous time series to calculate tmdo
            TMDON = 0.0
            TMDOP = 0.0
            IF(IMAX .GT. 1) THEN
              DO NG = 1, NSEGMENTS
                ISTART = IFRST(NG)
                IEND = ILAST(NG)
                DIF = TP(IEND)-TP(ISTART)
                IF(DIF .GT. 48.0) THEN
                  NMAX = IEND-ISTART+1  
                  DO I = ISTART, IEND
                    I0 = I-ISTART+1 
                    XTMP(I0) = TP(I)
                    YTMP(I0) = B(I)
                  ENDDO
                  TMDO0 = TMDO(YTMP,XTMP,IMX,NMAX,-2.0*ERR0)
                  TMDO1 = TMDO(YTMP,XTMP,IMX,NMAX,2.0*ERR0)
                  IF(TMDO0 .GT. TMDON) TMDON = TMDO0
                  IF(TMDO1 .GT. TMDOP) TMDOP = TMDO1
                ENDIF
              ENDDO
            ENDIF
            WRITE(LL,130) TTL,TTL2,IMAX,SM0,RMSE0,SD0,VOF0,CF0,
     1         POF0,TMDON,TMDOP,WCOF,CORR_C,SKILLV
          ENDIF

        ELSE IF(INDX .EQ. 6) THEN
          IF(IMAX .GT. 0) THEN
            SM0 = SM(B,IMX,IMAX)
            RMSE0 = RMSE(B,IMX,IMAX)
            SD0 = SD(B,IMX,IMAX)
            VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
            CF0 = CF(B,IMX,IMAX,ERR0)
            POF0 = POF(B,IMX,IMAX,2.0*ERR0)

CC  Using continuous time series to calculate tmdo
            TMDON = 0.0
            TMDOP = 0.0
            IF(IMAX .GT. 1) THEN
              DO NG = 1, NSEGMENTS
                ISTART = IFRST(NG)
                IEND = ILAST(NG)
                DIF = TP(IEND)-TP(ISTART)
                IF(DIF .GT. 48.0) THEN
                  NMAX = IEND-ISTART+1  
                  DO I = ISTART, IEND
                    I0 = I-ISTART+1 
                    XTMP(I0) = TP(I)
                    YTMP(I0) = B(I)
                  ENDDO
                  TMDO0 = TMDO(YTMP,XTMP,IMX,NMAX,-2.0*ERR0)
                  TMDO1 = TMDO(YTMP,XTMP,IMX,NMAX,2.0*ERR0)
                  IF(TMDO0 .GT. TMDON) TMDON = TMDO0
                  IF(TMDO1 .GT. TMDOP) TMDOP = TMDO1
                ENDIF
              ENDDO
            ENDIF
            WRITE(LL,131) TTL,TTL2,IMAX,SM0,RMSE0,SD0,VOF0,CF0,
     1         POF0,TMDON,TMDOP,CORR_C,SKILLV
          ENDIF
        ENDIF
130     FORMAT(A9,1X,A15,I6,1X,3F7.3,3F6.1,1X,F6.1,1x,F5.1,1x,F5.1,F6.2,1x,F6.2) 
131     FORMAT(A9,1X,A15,I6,1X,3F7.3,3F6.1,1X,F6.1,1x,F5.1,6X,F6.2,F6.2) 
150   CONTINUE

      RETURN
      END


C  Given a series of tidal extrema and times, compute error 
C    time series and store in 'b' 
C  Variables:
C    indx = index for computing: 1=amplitudes, 2=times
C    hw1  = reference ampl
C    hw2  = predicted ampl
C    tw1  = reference time (hrs)
C    tw2  = predicted time
      SUBROUTINE EXTEMS(HW1,TW1,NMAX1,HW2,TW2,NMAX2,B,XTMP,NMAXB,INDX)
      include 'skills.inc'
      DIMENSION HW1(NMX),TW1(NMX),HW2(NMX),TW2(NMX),B(NMX),XTMP(IMX)

      IPR = 0
      IF(IPR .EQ. 1) THEN
        WRITE(*,"(/,' <extrems>')")
        DO N = 1, NMAX1
          WRITE(*,"(' set 1: n=',i3,' tw=',f8.3,' hw=',f6.3)")
     1       N,TW1(N),HW1(N)
        ENDDO
        DO N = 1, NMAX2
          WRITE(*,"(' set 2: n=',i3,' tw=',f8.3,' hw=',f6.3)")
     1       N,TW2(N),HW2(N)
        ENDDO
      ENDIF

C  Loop thru series
      DELTA = 3   ! changed from 15 hrs to 3 hrs on 5/5/2008 
      M = 0
      DO N = 1, NMAX1
        NC = 0
        TIMEMIN = 999.0
        DO N2 = 1, NMAX2
          DIFFF = ABS(TW1(N)-TW2(N2))
          IF(DIFFF .LT. TIMEMIN) THEN
            TIMEMIN = DIFFF
            NC = N2
          ENDIF
        ENDDO

        IF((TIMEMIN .LT. DELTA) .AND. (NC .GT. 0)) THEN 
          M = M+1
          XTMP(M) = TW1(N)
          IF(INDX .EQ. 1) B(M) = HW2(NC)-HW1(N)
          IF(INDX .EQ. 2) B(M) = TW2(NC)-TW1(N)
          IF(IPR .EQ. 1) THEN
            WRITE(*,"('  m=',i3,' tw1,2=',2f8.3,' hw1,2=',2f6.3,
     1       ' diff=',f6.3)") M,TW1(N),TW2(NC),HW1(N),HW2(NC),B(M)
          ENDIF
        ENDIF
      ENDDO
      NMAXB = M

      RETURN
      END


C  Construct observed data time series corresponding to the nowcast/forecast data
      SUBROUTINE COLLECT_OBS(OP,AP,DIROP,DIRAP,TP,TMP,IPMAX,NMAX)
      include 'skills.inc'
      DIMENSION OP(IMX),AP(IMX),TP(IMX),TMP(IMX),DIROP(IMX),DIRAP(IMX)

      WRITE(*,"(' <Collect_obs>')")
C      WRITE(*,*) 'The number of prediction data ipmax= ', IPMAX
      IPR = 0
      N = 0
      DO I = 1, IPMAX
        KK = 0
        DO II = 1, IMAXB
          IF(ABS(TP(I)-T(II)) .LT. DELT*0.1) THEN
            N = N+1
            OP(N) = OF(II)
            AP(N) = A(II)
            DIROP(N) = DIRO(II)
            DIRAP(N) = DIRA(II)
            TMP(N) = TP(I)
            KK = II
            GOTO 100
          ENDIF
        ENDDO

100     CONTINUE
        IF(KK .EQ. 0) THEN
          WRITE(*,*) 'Match point of obs not found',TP(I),T(1),T(IMAXB)
          N = N+1
          OP(N) = -999.0
          AP(N) = -999.0
          DIROP(N) = -999.0
          DIRAP(N) = -999.0
        ENDIF
      ENDDO
      NMAX = N

      IF(IPR .GE. 1) THEN
        DO N = 1, NMAX
          WRITE(*,200) N,TP(N)/24.0,OP(N),AP(N)
        ENDDO
      ENDIF
200   FORMAT(' n=',i6,' t=',f9.4,' op=',f10.3,' ap=',f10.3)

      RETURN
      END


C  Construct project forecast time series of Group 3
C  Variables-
C     indx1 = variable selection: 1=ncst, 2=fcst, 3=persistence
C     indx2 = mode of comparison: 
C             1: each time in t(i), 
C             2: each 24/(jmax-1) hrs and at projection jj
C     jj = projection time: 00, 24/(jmax-1) hr, 2*24/(jmax-1) hr, etc
C     jj = 0:00Z; jj = 1:06Z; jj = 2:12Z; jj = 3:18Z; jj = 4:24Z
      SUBROUTINE ADDPRED(INDX1,INDX2,P,DIRP,TP,IPMAX,JJ)
      include 'skills.inc'
      DIMENSION P(IMX),TP(IMX),DIRP(IMX)

      WRITE(*,"(/,' <Sub ADDPRED> ','jj=',i2,' indx1=',i1,
     1   ' indx2=',i1)") JJ,INDX1,INDX2

C  Set print index
      IPR = 0
      IF(INDX1 .EQ. 1) IPR = 2  ! ncst
      IF(INDX1 .EQ. 2) IPR = 1  ! fcst
      IF(INDX1 .EQ. 3) IPR = 0  ! persistence

C  Loop thru preds, compare with valid times
C  Modified by zaj on Dec. 16, 2004
      IF(INDX2 .EQ. 1) THEN
        INTV = 1
        I1 = 1
        J1 = 1
        J2 = JMAX
        NNN = INT(IPDAY/JMAX)
        IP = 0
        DO K = 1, IMAXF
          DO J = 1, 1  !Modified by AJ concatenate only 1st cycle's fcst
            DO II = 1, IPDAY !!!   NNN
              IF(INDX1 .EQ. 2) THEN
                VAL  = F(II,J,K)
                VAL2 = DIRF(II,J,K)
              ELSE IF(INDX1 .EQ. 3) THEN
                VAL  = PF(II,J,K)
                VAL2 = DIRPF(II,J,K)
              END IF

              IP = IP+1
              P(IP) = VAL
              DIRP(IP) = VAL2
              TP(IP) = TF(II,J,K)
            ENDDO  
          ENDDO
        ENDDO
        IPMAX = IP
      ELSE IF(INDX2 .EQ. 2) THEN  !! each 24/jmax hrs and at projection jj  
        IP = 0
        IF(FINCLUDE0) THEN
          NPROJECT_TIME = JJ*6/DELT+1 !! do 00z 06z 12z 18z 24z  AJ
        ELSE
          NPROJECT_TIME = JJ*6/DELT
        END IF

        DO 151 K = 1, IMAXF
          DO 141 J = 1, JMAX
            IF(INDX1 .GE. 2) THEN
              TT = TF(NPROJECT_TIME,J,K)
            END IF
            IF(INDX1 .EQ. 2) THEN
              VAL  = F(NPROJECT_TIME,J,K)
              VAL2 = DIRF(NPROJECT_TIME,J,K)
            ELSE IF(INDX1 .EQ. 3) THEN
              VAL  = PF(NPROJECT_TIME,J,K)
              VAL2 = DIRPF(NPROJECT_TIME,J,K)
            END IF

            IP = IP+1
            P(IP) = VAL
            DIRP(IP) = VAL2
            TP(IP) = TT
            TMATCH = TT
141       CONTINUE
151     CONTINUE
        IPMAX = IP
      ENDIF
C300   FORMAT(' ip=',I6,' t=',F9.3,' p=',4F6.3)

      RETURN
      END


C  Get high waters
      SUBROUTINE GETHI(IPR,Z,HHW,THW,NMAX,NMAXRH0)
      include 'skills.inc'
      DIMENSION HHW(NMX),THW(NMX),Z(IMX),SMO(IMX)
      CHARACTER CHARR*1

      WRITE(*,"(/,' <gethi>')")
      I1 = 1
      I2 = NMAX 
      WRITE(*,"(/,' i1,i2=',2I5)") I1,I2

C  Pick off peaks
      NM = 0
      DO I = I1+1,I2-1 
        CHARR = ' '
        IPICK = 0
        IF(I .EQ. I1 .AND. Z(I) .GT. Z(I+1)) IPICK = I
        IF(I .EQ. I2 .AND. Z(I) .GT. Z(I-1)) IPICK = I
        IF(I .GT. I1 .AND. I. LE. I2-2 .AND. (Z(I) .GT. Z(I-1) 
     1    .AND. Z(I+1) .GT. Z(I+2) .AND. Z(I) .EQ. Z(I+1))) IPICK = I
        IF(I .GT. I1 .AND. I .LT. I2 .AND. (Z(I) .GT. Z(I-1)
     1    .AND. Z(I) .GT. Z(I+1))) IPICK = I
        IF(I .EQ. I1+1 .AND. Z(I) .EQ. Z(I-1)) IPICK = I
        IF(I .EQ. I2-1 .AND. Z(I) .EQ. Z(I+1)) IPICK = I

        IF(IPICK .GT. 0) THEN
          CHARR = '*'
          NM = NM+1
          HHW(NM) = Z(I)
          THW(NM) = T(I)
        ENDIF

        IF(IPR .EQ. 1) THEN
          WRITE(*,"(' i=',I7,' t=',F12.5,'   z=',F9.4,3x,a1)")
     1      I,T(I),Z(I),CHARR
        ENDIF
      ENDDO

C  Print final status
      NMAXRH0 = NM
      WRITE(*,"(' nmaxrh=',I4)") NMAXRH0
      IF(NMAXRH0 .GT. 0) THEN
        DO N = 1,NMAXRH0
          WRITE(*,"(' n=',I3,' thw=',F8.3,' hhw=',F6.3)") 
     1       N,THW(N),HHW(N)
        ENDDO
      ENDIF

      RETURN
      END


C  Get low  waters
      SUBROUTINE GETLO(IPR,Z,HLW,TLW,NMAX,NMAXRL0)
      include 'skills.inc'
      DIMENSION HLW(NMX),TLW(NMX),Z(IMX),SMO(IMX)
      CHARACTER CHARR*1

      WRITE(*,"(/,' <getlo>')")
      I1 = 1
      I2 = NMAX
      WRITE(*,"(/,' i1,i2=',2I5)") I1,I2

C  Pick off peaks
      NM = 0
      DO I = I1+1,I2-1 
        CHARR = ' '
        IPICK = 0
        IF(I .EQ. I1 .AND. Z(I) .LT. Z(I+1)) IPICK = I
        IF(I .EQ. I2 .AND. Z(I) .LT. Z(I-1)) IPICK = I
        IF(I .GT. I1 .AND. I .LE. I2-2 .AND. (Z(I) .LT. Z(I-1)
     1    .AND. Z(I+1) .LT. Z(I+2) .AND. Z(I) .EQ. Z(I+1))) IPICK = I
        IF(I .GT. I1 .AND. I .LT. I2 .AND. (Z(I) .LT. Z(I-1)
     1    .AND. Z(I) .LT. Z(I+1))) IPICK = I
        IF(I .EQ. I1+1 .AND. Z(I) .EQ. Z(I-1)) IPICK = I
        IF(I .EQ. I2-1 .AND. Z(I) .EQ. Z(I+1)) IPICK = I

        IF(IPICK .GT. 0) THEN
          CHARR = '*'
          NM = NM+1
          HLW(NM) = Z(I)
          TLW(NM) = T(I)
        ENDIF

        IF(IPR .EQ. 1) THEN
          WRITE(*,"(' i=',I7,' t=',F12.5,' z=',F8.4,3x,a1)") 
     1      I,T(I),Z(I),CHARR
        ENDIF
      ENDDO

C  Print final status
      NMAXRL0 = NM
      WRITE(*,"('nmaxrl =',I4)") NMAXRL0
      IF(NMAXRL0 .GT. 0) THEN
        DO N = 1,NMAXRL0
          WRITE(*,"(' n=',I3,' tlw=',F8.3,' hlw=',F6.3)") 
     1       N,TLW(N),HLW(N)
        ENDDO
      ENDIF

      RETURN
      END


C  Maximize time duration of (pos or neg) outliers
      FUNCTION TMDO(Z,TZ,IMX,IMAX,ERR0)
      DIMENSION Z(IMX),TZ(IMX),M(IMX),IFRST(IMX),ILAST(IMX),TMP(IMX)

      DO I = 1, IMAX
        TMP(I) = Z(I)
        IF(Z(I) .LT. -990.0) TMP(I) = 0.0
      ENDDO

      IPR = 0
      IF(IPR .EQ. 1) THEN
        WRITE(*,"(/,' mdo. error=',f8.3)") ERR0
      END IF

C  Tag outliers with a '1'
      DO I = 1,IMAX
        M(I) = 0
        IF(ERR0 .LT. 0.0 .AND. TMP(I) .LT. ERR0) M(I) = 1     
        IF(ERR0 .GT. 0.0 .AND. TMP(I) .GT. ERR0) M(I) = 1     
        IF(IPR .EQ. 1) THEN
          WRITE(*,"(' i=',I4,' m=',I2)") I,M(I)
        END IF
      ENDDO

C  Look for start of a series of '1's.  ifrst is i of first '1'
      IS = 0
      DO I = 1,IMAX-1
        IF(I .EQ. 1 .AND. M(I) .EQ. 1) THEN
          IS = IS+1
          IFRST(IS) = I
        ENDIF
        IF(I. GE. 1 .AND. M(I) .EQ. 0 .AND. M(I+1) .EQ. 1) THEN
          IS = IS+1
          IFRST(IS) = I+1
        ENDIF
      ENDDO
      IF(IPR .EQ. 1) THEN
        WRITE(*,*) ' is=', IS
        IF(IS .GT. 0) THEN
          DO I = 1, IS
            WRITE(*,*) ' i=', I, ' ifrst=', IFRST(I)
          ENDDO
        ENDIF
      ENDIF

C  Look for end of a series of '1's.  ilast is i of last '1'
      IE = 0
      DO I = 1,IMAX-1
        IF(M(I) .EQ. 1 .AND. M(I+1) .EQ. 0) THEN
          IE = IE+1
          ILAST(IE) = I
        ENDIF
      ENDDO

      IF(M(IMAX) .EQ. 1) THEN
        IE = IE+1
        ILAST(IE) = I
      ENDIF

      IF(IPR .EQ. 1) THEN
        WRITE(*,*) ' ie=', IE
        IF(IE .GT. 0) THEN
          DO I = 1, IE
            WRITE(*,*) ' i=', I, ' ilast=', ILAST(I)
          ENDDO
        ENDIF
      ENDIF

C  Find longest duration
      TMDO = 0
      IF(IS .NE. IE) THEN
        WRITE(*,"(/,3x,'**function tmdo. starts not equal to ends**')")
        WRITE(*,"(5x,'is=',I4,' ie=',I4)") IS,IE
        IF(IS .GT. 0) THEN
          DO I = 1, IS
            WRITE(*,*)' i=', I, ' ifrst=', IFRST(I)
          ENDDO
        ENDIF
        IF(IE .GT. 0) THEN
          DO I = 1, IE
            WRITE(*,*)' i=', I, ' ilast=', ILAST(I)
          ENDDO
        ENDIF
        DO I = 1, IMAX
          WRITE(*,"(' i=',i4,' m=',i2)") I, M(I)
        ENDDO
        STOP
      ELSE IF(IS .GT. 0) THEN
        DO I = 1, IS
          TMDO = AMAX1(TMDO,TZ(ILAST(I))-TZ(IFRST(I)))
        ENDDO
      ENDIF

      RETURN
      END

      FUNCTION SKILLA(P,O,IMX,IMAX)
      DIMENSION P(IMX),O(IMX)

      NSUM = 0
      SM0  = SM(O,IMX,IMAX)
      SUM1 = 0.0
      SUM2 = 0.0
      DO I = 1,IMAX
        IF((ABS(P(I)) .LT. 900.0) .AND.(ABS(O(I)) .LT. 900.0)) THEN
          SUM1 = SUM1+(P(I)-O(I))**2
          SUM2 = SUM2+(ABS(P(I)-SM0)+ABS(O(I)-SM0))**2
          NSUM = NSUM+1
        ENDIF
      ENDDO

      IF(ABS(SUM2) .GT. 0.0) THEN
        SKILLA = 1-SUM1/SUM2
      ELSE
        SKILLA = 0.0
      ENDIF

      RETURN
      END

      FUNCTION CORRELATION(P,O,IMX,IMAX)
      DIMENSION P(IMX),O(IMX)

      NSUM = 0
      SM0  = SM(O,IMX,IMAX)
      SMP  = SM(P,IMX,IMAX)
      SUM1 = 0.0
      SUM2 = 0.0
      SUM3 = 0.0
      DO I = 1,IMAX
        IF((ABS(P(I)) .LT. 900.0) .AND.(ABS(O(I)) .LT. 900.0)) THEN
          SUM1 = SUM1+ (P(I)-SMP)*(O(I)-SMO)
          SUM2 = SUM2+(P(I)-SMP)**2
          SUM3 = SUM3+(O(I)-SM0)**2
          NSUM = NSUM+1
        ENDIF
      ENDDO
      SQRT23=SQRT(SUM2 * SUM3)
      IF (ABS(SQRT23) .GT. 0.001)THEN
         CORRELATION=SUM1/SQRT23
      ELSE
         CORRELATION=0.0
      ENDIF


      RETURN
      END


C  Worst cast outlier frequency. p=modeled, a=astro tide, o=obs
      FUNCTION WOF(P,A,O,IMX,IMAX,ERR0)
      DIMENSION P(IMX),A(IMX),O(IMX)

      NSUM = 0
      DO I = 1,IMAX
        IF(ABS(P(I)-O(I)) .GT. ERR0) THEN
          IF(P(I) .GT. A(I) .AND. O(I) .LT. A(I)) NSUM = NSUM+1
          IF(P(I) .LT. A(I) .AND. O(I) .GT. A(I)) NSUM = NSUM+1
        ENDIF
      ENDDO
      WOF = 100.0*FLOAT(NSUM)/FLOAT(IMAX)

      RETURN
      END


C  Central frequency
      FUNCTION CF(Z,IMX,IMAX,ERR0)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO

      NSUM = 0
      DO I = 1,N
        IF(TMP(I) .GT. -ERR0 .AND. TMP(I) .LT. ERR0) NSUM = NSUM+1
      ENDDO
      CF = 100.0*FLOAT(NSUM)/FLOAT(N)

      RETURN
      END


C  Positive outlier frequency
      FUNCTION POF(Z,IMX,IMAX,ERR0)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO

      NSUM = 0
      DO I = 1,N
        IF(TMP(I) .GT. ERR0) NSUM = NSUM+1
      ENDDO
      POF = 100.0*FLOAT(NSUM)/FLOAT(N)

      RETURN
      END


C  Negative outlier frequency
      FUNCTION VOF(Z,IMX,IMAX,ERR0)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO
      
      NSUM = 0
      DO I = 1,N
        IF(TMP(I) .LT. -ERR0) NSUM = NSUM+1
      ENDDO
      VOF = 100.0*FLOAT(NSUM)/FLOAT(N)

      RETURN
      END


C  Series mean
      FUNCTION SM(Z,IMX,IMAX)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO

      IF(N .LT. 1) THEN
        SM = 0.0
      ELSE
        SUM = 0.0
        DO I = 1,N
          SUM = SUM+TMP(I)
        ENDDO
        SM = SUM/N
      ENDIF

      RETURN
      END


C  Standard deviation
      FUNCTION SD(Z,IMX,IMAX)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO
      SUM = SM(TMP,IMX,N)

      IF(N .LT. 2) THEN
        SD = 0.0
      ELSE
        S2 = 0.0
        DO I = 1,N
          S2 = S2+(TMP(I)-SUM)**2
        ENDDO
        SD = SQRT(S2/(N-1))
      ENDIF

      RETURN
      END


C  Root mean squared error
      FUNCTION RMSE(Z,IMX,IMAX)
      DIMENSION Z(IMX),TMP(IMX)

      N = 0
      DO I = 1,IMAX
        IF(Z(I) .GT. -999.0) THEN
          N = N+1
          TMP(N) = Z(I)
        ENDIF
      ENDDO

      IF(N .LT. 1) THEN
        RMSE = 0.0
      ELSE
        S2 = 0.0
        DO I = 1,N
          S2 = S2+TMP(I)**2
        ENDDO
        RMSE = SQRT(S2/N)
      ENDIF

      RETURN
      END


      SUBROUTINE READDATA(KINDAT,FILESHORT,DATADIR)
      include 'skills.inc'

      DIMENSION YTMP(IMX),YNEW(IMX),TMP(2,IMX),
     1  B(IMX),P2(IMX),U(IMX),V(IMX),UTMP(IMX),VTMP(IMX),WLTMP(99),
     2  F_STATMP(2,MAXDURATION,JMX,KMX)
      REAL*8 TFTMP(MAXDURATION,JMX,KMX)
      REAL*8 TTMP,XNEW(IMX),TSTART0,TFINISH0,XTMP(IMX)
      CHARACTER*200 FILESHORT,FILEIN,DATADIR
      LOGICAL FEXIST

      PII = 3.1415936/180.0
      DO IJK = 1, IMX
        U(IJK)    = -999.0
        V(IJK)    = -999.0
        XNEW(IJK) = -999.0
        YNEW(IJK) = -999.0
        UTMP(IJK) = -999.0
        VTMP(IJK) = -999.0
      END DO

C  Reading water levels and currents in tidal regions for tide prediction only a(i)
      IF(NTYPE .NE. 0 .AND. KINDAT .LE. 2) THEN
!        FILEIN = '../data/prediction/'//TRIM(FILESHORT)//'.prd'
        FILEIN = trim(DATADIR)//'/prediction/'//TRIM(FILESHORT)//'.prd'
        INQUIRE(FILE=TRIM(FILEIN),EXIST=FEXIST)
        IF(.NOT. FEXIST) THEN
          WRITE(*,*) 'filein=',TRIM(FILEIN), ' does not exist!'
          WRITE(*,*) 'Skill assessment table is not generated.' 
          WRITE(*,*) 'Tidal prediction has to be existed in order to'
          WRITE(*,*) 'do skill assessment in tidal regions ...'
          RETURN
        ENDIF

        I = 0
        OPEN(10,FILE=TRIM(FILEIN),FORM='FORMATTED',STATUS='OLD')
        IF(KINDAT .EQ. 1) THEN
          DO IJK = 1, 200000
            READ(10,*,END=20) TTMP, IYR0, ICM, ICD, IHR, IMN,
     1        (WLTMP(N),N = 1, 2)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I    = I + 1
              T(I) = TTMP * 24.0
              U(I) = WLTMP(1)*SIN(WLTMP(2)*PII)
              V(I) = WLTMP(1)*COS(WLTMP(2)*PII)
            END IF
          END DO
20        CLOSE(10)

          IMAXT = I
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_T,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,T,U,XNEW,YNEW,IMAXT,IMAXB)
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_T,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,T,V,XNEW,YTMP,IMAXT,IMAXB)

          DO J1 = 1, IMAXB
            T(J1) = XNEW(J1)/24.0
            IF((YNEW(J1) .LT. -900.0) .OR. 
     1         (YTMP(J1) .LT. -900.0)) THEN
              A(J1)    = -999.0
              DIRA(J1) = -999.0
            ELSE   
              CALL VELDIR(YTMP(J1),YNEW(J1),DR,SP)
              IF(DR .GT. 360.0) DR = DR - 360.0
              A(J1)    = SP
              DIRA(J1) = DR
            ENDIF
          ENDDO
        ELSE IF(KINDAT .GE. 2) THEN
          DO IJK = 1, 200000
            READ(10,*,END=25) TTMP,IYR0,ICM,ICD,IHR,IMN,WLTMP(1)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I    = I + 1
              T(I) = TTMP * 24.0
              U(I) = WLTMP(1)
            END IF
          END DO
25        CLOSE(10)

          IMAXT = I
          WRITE(*,*) 'TSTART: ',TSTART
          WRITE(*,*) 'TFINISH: ',TFINISH
          WRITE(*,*) 'DELT, DLET_T:',DELT,DELT_T
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_T,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,T,U,XNEW,YNEW,IMAXT,IMAXB)
          
          DO J1 = 1, IMAXB
             T(J1) = XNEW(J1)/24.0 
             A(J1) = YNEW(J1)
          ENDDO
        ENDIF
        WRITE(*,*) 'The number of tide prediction data, imaxt=',IMAXT
      ELSE
        DO I = 1, IMAXB
          T(I) = TSTART + (I-1)*DELT/24.0
          A(I) = -999.0
          DIRA(I) = -999.0
        ENDDO
      ENDIF

C  Reading observation data o(i)
!      FILEIN = '../data/obs/'//TRIM(FILESHORT)//'.obs'
      FILEIN = trim(DATADIR)//'/obs/'//TRIM(FILESHORT)//'.obs'
C  If only conducting skill assessment for scenario of tidal 
C    simulation, using prediction to replace the obs to perform 
C    skill assessment for those stations without observations 
C    as long as having tidal predictions
      IF(ISWITCH(1) .GT. 0 .AND. ISWITCH(2) .EQ. 0 .AND.
     1   ISWITCH(3) .EQ. 0 .AND. ISWITCH(4) .EQ. 0 .AND.
     1   ISWITCH(5) .EQ. 0 ) THEN 
        FILEIN = trim(DATADIR)//'/prediction/'//TRIM(FILESHORT)//'.prd'
      ENDIF

      WRITE(*,*) 'obs = ', TRIM(FILEIN)
      INQUIRE(FILE=TRIM(FILEIN),EXIST=FEXIST)
      IF(.NOT. FEXIST) THEN
        WRITE(*,*) 'filein=',TRIM(FILEIN), ' does not exist!'
        IMAXB = 0
        RETURN
      ENDIF

      I = 0
      OPEN(10, FILE=TRIM(FILEIN),FORM='FORMATTED',STATUS='OLD')
      IF(KINDAT .EQ. 1) THEN
        DO IJK = 1, 200000
          READ(10,*,END=40) TTMP, IYR0, ICM, ICD, IHR, IMN,
     1      (WLTMP(N),N = 1, 2)
          IF(TTMP .GE. TSTART .AND. TTMP .LE. (TFINISH+2) .AND.
     1       WLTMP(1) .GT. -90.0) THEN
            I = I + 1
            XTMP(I) = TTMP
            TMP(1,I) = WLTMP(1)
            TMP(2,I) = WLTMP(2)
            U(I) = WLTMP(1)*SIN(WLTMP(2)*PII)
            V(I) = WLTMP(1)*COS(WLTMP(2)*PII)
          END IF
        END DO
40      CLOSE(10)

        IMAXO = I
        NMAX  = IMAXO
        WRITE(*,*) 'The number of observation data, IMAXO=',IMAXO
        WRITE(*,*) 'First data=', XTMP(1),TMP(1,1),TMP(2,1)
        WRITE(*,*) 'Last data=',  XTMP(NMAX),TMP(1,NMAX),TMP(2,NMAX)

        IF(TCUT .GT. 0.0) THEN
          CALL FOUFIL0(IMX,NMAX,DELT_O*60.0,TCUT,U,UTMP)
          CALL FOUFIL0(IMX,NMAX,DELT_O*60.0,TCUT,V,VTMP)
          DO I = 1, NMAX
            AAVN = VTMP(I)
            AAVE = UTMP(I)
            CALL VELDIR(AAVN,AAVE,AANGLE,AVEL)
            IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
            TMP(1,I) = AVEL
            TMP(2,I) = AANGLE
          ENDDO
        ENDIF  

      ELSE IF(KINDAT .GE. 2) THEN
        DO IJK = 1, 200000
          READ(10,*,END=45) TTMP,IYR0,ICM,ICD,IHR,IMN,WLTMP(1)
          IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1       WLTMP(1) .GT. -90.0) THEN
            I = I + 1
            XTMP(I) = TTMP
            TMP(1,I) = WLTMP(1)
            TMP(2,I) = 0.0
          END IF
        END DO
45      CLOSE(10)

        IMAXO = I
        NMAX = IMAXO
        WRITE(*,*) 'The number of observation data, IMAXO=',IMAXO
        WRITE(*,*) 'First data=',XTMP(1),TMP(1,1)
        WRITE(*,*) 'Last data=', XTMP(NMAX),TMP(1,NMAX)

        IF(TCUT .GT. 0.0) THEN
          DO I = 1,NMAX
            U(I) = TMP(1,I)
          ENDDO
          CALL FOUFIL(NMAX,DELT_O*60.0,TCUT,U,UTMP)
          DO I = 1, NMAX
            TMP(1,I) = UTMP(I)
          ENDDO
          IMAXO = NMAX
        ENDIF
      ENDIF

      NMAX = IMAXO
      IF(NMAX .LE. 1) THEN
        DO J1 = 1, IMAXB
          O(J1) = -999.0
          DIRO(J1) = -999.0
        ENDDO
      ELSE  
        DO J1 = 1,NMAX
          TC(J1) = XTMP(J1)*24.0
        ENDDO
        IF(KINDAT .EQ. 1) THEN
          DO I = 1, NMAX
            SP = TMP(1,I)
            DR = TMP(2,I)
            U(I) = SP*SIN(DR*PII)
            V(I) = SP*COS(DR*PII)
          ENDDO 
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_O,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_O,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NMAX,IMAXB)
          DO J1 = 1, IMAXB
            IF((YNEW(J1) .LT. -900.0) .OR. 
     1         (YTMP(J1) .LT. -900.0)) THEN
              O(J1) = -999.0
              DIRO(J1) = -999.0
            ELSE   
              CALL VELDIR(YTMP(J1),YNEW(J1),DR,SP)
              IF(DR .GT. 360.0) DR = DR - 360.0
              O(J1) = SP
              DIRO(J1) = DR
            ENDIF
          ENDDO
        ELSE IF(KINDAT .GE. 2) THEN
          DO I = 1, NMAX
            U(I) = TMP(1,I)
          ENDDO 
          CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_O,METHOD,
     1         CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
          DO J1 = 1, IMAXB
            O(J1)=YNEW(J1)
            IF(ISURGE .GT. 0 .AND. O(J1) .GT. -900.0 .AND.
     1         A(J1) .GT. -900.0) THEN
              SGO(J1) = O(J1)-A(J1)
            END IF
          ENDDO
        ENDIF 
      ENDIF 

CC  Reading model tidal simulation s(i)
      IF(ISWITCH(1) .GT. 0) THEN 
        FILEIN = TRIM(FILESHORT)//'_modeltides.dat'
        OPEN(10,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD')
        I = 0
        IF(KINDAT .EQ. 1) THEN
          DO IJK = 1, 200000
            READ(10,*,END=60) TTMP,IYR0,ICM,ICD,IHR,IMN,
     1        (WLTMP(N),N=1,4)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = WLTMP(2)
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            END IF
          END DO
60        CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,U,UTMP)
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,V,VTMP)
            DO I = 1,NMAX
              AAVN = VTMP(I)
              AAVE = UTMP(I)
              CALL VELDIR(AAVN,AAVE,AANGLE,AVEL)
              IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
              TMP(1,I) = AVEL
              TMP(2,I) = AANGLE
            ENDDO
            IMAXO = NMAX
          ENDIF  
        ELSE IF(KINDAT .GE. 2) THEN
          DO IJK = 1, 200000
            READ(10,*,END=65)TTMP,IYR0,ICM,ICD,IHR,IMN,WLTMP(1)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = 0.0
            END IF
          END DO
65        CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            DO I = 1,NMAX
              U(I) = TMP(1,I)
            ENDDO
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,U,UTMP)
            DO I = 1, NMAX
              TMP(1,I) = UTMP(I)
            ENDDO
            IMAXO = NMAX
          ENDIF
        ENDIF
        NMAX = IMAXO
        WRITE(*,*) 'The number of tidal simulation, IMAX=',IMAXO

        IF(NMAX .LE. 1) THEN
          DO J1 = 1, IMAXB
            S(J1) = -999.0
            DIRS(J1) = -999.0
          ENDDO
        ELSE  
          DO J1 = 1, NMAX
            TC(J1) = XTMP(J1)*24.0
          ENDDO
          IF(KINDAT .EQ. 1) THEN
            DO I = 1, NMAX
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NMAX,IMAXB)

            DO J1 = 1, IMAXB
              IF((YNEW(J1) .LT. -900.0) .OR. 
     1           (YTMP(J1) .LT. -900.0)) THEN
                S(J1) = -999.0
                DIRS(J1) = -999.0
              ELSE   
                CALL VELDIR(YTMP(J1),YNEW(J1),DR,SP)
                IF(DR .GT. 360.0) DR = DR - 360.0
                S(J1) = SP
                DIRS(J1) = DR
              ENDIF
            ENDDO
          ELSE IF(KINDAT .GE. 2) THEN
            DO I = 1, NMAX
              U(I) = TMP(1,I)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            DO J1 = 1, IMAXB
              S(J1) = YNEW(J1)
            ENDDO
          ENDIF 
        ENDIF 
      ENDIF

CC  Reading model hindcast hi(i)
      IF(ISWITCH(2) .GT. 0) THEN 
        FILEIN = TRIM(FILESHORT)//'_hindcast.dat'
        OPEN(10,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD')
        I = 0
        IF(KINDAT .EQ. 1) THEN
          DO IJK = 1, 200000
            READ(10,*,END=80) TTMP,IYR0,ICM,ICD,IHR,IMN,
     1        (WLTMP(N),N=1,4)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = WLTMP(2)
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            END IF
          END DO
80        CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,U,UTMP)
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,V,VTMP)
            DO I = 1, NMAX
              AAVN = VTMP(I)
              AAVE = UTMP(I)
              CALL VELDIR(AAVN,AAVE,AANGLE,AVEL)
              IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
              TMP(1,I) = AVEL
              TMP(2,I) = AANGLE
            ENDDO
            IMAXO = NMAX
          ENDIF  
        ELSE IF(KINDAT .GE. 2) THEN
          DO IJK = 1, 200000
            READ(10,*,END=85) TTMP,IYR0,ICM,ICD,IHR,IMN,WLTMP(1)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = 0.0
            END IF
          END DO
85        CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            DO I = 1, NMAX
              U(I) = TMP(1,I)
            ENDDO
            CALL FOUFIL(NMAX,DELT_M*60.0,TCUT,U,UTMP)
            DO I = 1, NMAX
              TMP(1,I) = UTMP(I)
            ENDDO
            IMAXO = NMAX
          ENDIF
        ENDIF
        NMAX = IMAXO 
        WRITE(*,*) 'The number of hindcast simulation, IMAX=',IMAXO

        IF(NMAX .LE. 1) THEN
          DO J1 = 1, IMAXB
            HI(J1) = -999.0
            DIRHI(J1) = -999.0
          ENDDO
        ELSE  
          DO J1 = 1, NMAX
            TC(J1) = XTMP(J1)*24.0
          ENDDO
          IF(KINDAT .EQ. 1) THEN
            DO I = 1, NMAX
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NMAX,IMAXB)
            DO J1 = 1, IMAXB
              IF((YNEW(J1) .LT. -900.0) .OR. 
     1           (YTMP(J1) .LT. -900.0)) THEN
                HI(J1) = -999.0
                DIRHI(J1) = -999.0
              ELSE   
                CALL VELDIR(YTMP(J1),YNEW(J1),DR,SP)
                IF(DR .GT. 360.0) DR = DR - 360.0
                HI(J1) = SP
                DIRHI(J1) = DR
              ENDIF
            ENDDO
          ELSE IF(KINDAT .GE. 2) THEN
            DO I = 1, NMAX
              U(I) = TMP(1,I)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            DO J1 = 1, IMAXB
              HI(J1) = YNEW(J1)
              IF(ISURGE .GT. 0 .AND. HI(J1) .GT. -900.0 .AND.
     1           S(J1) .GT. -900.0) THEN
                SGM(J1) = HI(J1)-S(J1)
              END IF
            ENDDO
          ENDIF 
        ENDIF 
      ENDIF

CC  Reading model nowcast c(i),tc(i)
      IF(ISWITCH(3) .GT. 0) THEN 
        WRITE(*,*)
        FILEIN = TRIM(FILESHORT)//'_nowcast.dat'
        WRITE(*,*) 'nowcast = ./',TRIM(FILEIN)
        OPEN(10,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD')
        I = 0
        IF(KINDAT .EQ. 1) THEN
          DO IJK = 1, 200000
            READ(10,*,END=100) TTMP,IYR0,ICM,ICD,IHR,IMN,
     1        (WLTMP(N),N=1,4)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = WLTMP(2)
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            END IF
          END DO
100       CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            CALL FOUFIL(NMAX,real(DELT_M)*60.0,TCUT,U,UTMP)
            CALL FOUFIL(NMAX,real(DELT_M)*60.0,TCUT,V,VTMP)
            DO I = 1, NMAX
              AAVN = VTMP(I)
              AAVE = UTMP(I)
              CALL VELDIR(AAVN,AAVE,AANGLE,AVEL)
              IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
              TMP(1,I) = AVEL
              TMP(2,I) = AANGLE
            ENDDO
            IMAXO = NMAX
          ENDIF  
        ELSE IF(KINDAT .GE. 2) THEN
          DO IJK = 1, 200000
            READ(10,*,END=105)TTMP,IYR0,ICM,ICD,IHR,IMN,WLTMP(1)
            IF(TTMP .GE. TSTART .AND. TTMP .LE. TFINISH .AND.
     1         WLTMP(1) .GT. -90.0) THEN
              I = I+1
              XTMP(I) = TTMP
              TMP(1,I) = WLTMP(1)
              TMP(2,I) = 0.0
            END IF
          END DO
105       CLOSE(10)

          IMAXO = I
          NMAX = IMAXO
          IF(TCUT .GT. 0.0) THEN
            DO I = 1, NMAX
              U(I) = TMP(1,I)
            ENDDO
            CALL FOUFIL(NMAX,real(DELT_M)*60.0,TCUT,U,UTMP)
            DO I = 1, NMAX
              TMP(1,I) = UTMP(I)
            ENDDO
            IMAXO = NMAX
          ENDIF
        ENDIF
        NMAX = IMAXO 
        WRITE(*,*) 'The number of nowcast simulation, NMAX=',NMAX

        IF(NMAX .LE. 1) THEN
          DO J1 = 1,IMAXB
            C(J1) = -999.0
            DIRC(J1) = -999.0
          ENDDO
        ELSE  
          DO J1 = 1,NMAX
            TC(J1) = XTMP(J1)*24.0
          ENDDO
          IF(KINDAT .EQ. 1) THEN
            DO I = 1, NMAX
              SP = TMP(1,I)
              DR = TMP(2,I)
              U(I) = SP*SIN(DR*PII)
              V(I) = SP*COS(DR*PII)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NMAX,IMAXB)
            DO J1 = 1, IMAXB
              IF((YNEW(J1) .LT. -900.0) .OR. 
     1           (YTMP(J1) .LT. -900.0)) THEN
                C(J1) = -999.0
                DIRC(J1) = -999.0
              ELSE   
                CALL VELDIR(YTMP(J1),YNEW(J1),DR,SP)
                IF(DR .GT. 360.0) DR = DR - 360.0
                C(J1) = SP
                DIRC(J1) = DR
              ENDIF
            ENDDO
          ELSE IF(KINDAT .GE. 2) THEN
            DO I = 1, NMAX
              U(I) = TMP(1,I)
            ENDDO 
            CALL EQUAL_INTERVAL1(TSTART,TFINISH,DELT,DELT_M,METHOD,
     1           CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NMAX,IMAXB)
            DO J1 = 1, IMAXB
              C(J1) = YNEW(J1)
            ENDDO
          ENDIF 
        ENDIF 
      ENDIF

CC  Reading model forecast f(L,J,K), tf(L,J,K)
      IF(ISWITCH(4) .GT. 0) THEN 
        IF(FINCLUDE0) THEN
          NCUT0 = INT(NFDURATION/DELT_M+0.1)+1 ! fcst duration > 24h
          IDUMMY = INT(NFDURATION/DELT+0.1)+1  ! fcst duration > 24h
        ELSE
          NCUT0 = INT(NFDURATION/DELT_M+0.1)   ! fcst duration > 24h
          IDUMMY = INT(NFDURATION/DELT+0.1)    ! fcst duration > 24h
        END IF
        FILEIN = TRIM(FILESHORT)//'_forecast.dat'
        WRITE(*,*) 'forecast = ./',TRIM(FILEIN)
        OPEN(10,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD')
        IF(KINDAT .EQ. 1) THEN
          DO I = 1, NDAY
            DO J = 1, JMAX
              DO N = 1, NCUT0
                READ(10,*,END=120) TFTMP(N,J,I),IYR0,ICM,ICD,
     1            IHR,IMN,SP,DR
                F_STATMP(1,N,J,I) = SP
                F_STATMP(2,N,J,I) = DR
                U(N) = SP*SIN(DR*PII)
                V(N) = SP*COS(DR*PII)
              ENDDO
            ENDDO
          ENDDO 
        ELSE IF(KINDAT .GE. 2) THEN
          DO I = 1, NDAY
            DO J = 1, JMAX
              DO N = 1, NCUT0
                READ(10,*,END=120) TFTMP(N,J,I),IYR0,ICM,ICD,
     1            IHR,IMN,WLTMP(1)
                F_STATMP(1,N,J,I) = WLTMP(1)
                F_STATMP(2,N,J,I) = 0.0
                U(N) = WLTMP(1)
              ENDDO
            ENDDO
          ENDDO 
        ENDIF 
120     CONTINUE
        CLOSE(10)
        WRITE(*,*) 'NCUT0 JMAX NDAY: ',NCUT0,JMAX,NDAY

        IMAXF = I-1
        WRITE(*,*) 'The days of forecast simulation, IMAXF=',IMAXF
	WRITE(*,*) 'DELT=',DELT,'DELT_M=',DELT_M
        IF(ABS(DELT-DELT_M) .GT. 0.01) THEN 
          DO I = 1, IMAXF
            DO J = 1, JMAX
              TSTART0 = TFTMP(1,J,I)
              TFINISH0 = TFTMP(NCUT0,J,I)  !!  zaj 06/16/2008 tstart0+1
C              WRITE(*,*) 'TSTART0=',TSTART0,'  TFINISH0=',TFINISH0

              DO N = 1, NCUT0
                TC(N)=TFTMP(N,J,I)*dble(24.0)
              ENDDO

              IF(KINDAT .EQ. 1) THEN
                DO N = 1, NCUT0
                  SP = F_STATMP(1,N,J,I)
                  DR = F_STATMP(2,N,J,I)
                  U(N) = SP*SIN(DR*PII)
                  V(N) = SP*COS(DR*PII)
                ENDDO 
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NCUT0,IDUMMY)
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NCUT0,IDUMMY)
                DO N = 1, IDUMMY 
                  CALL VELDIR(YTMP(N),YNEW(N),DR,SP)
                  IF(DR .GT. 360.0) DR = DR - 360.0
                  F_STATMP(1,N,J,I)=SP
                  F_STATMP(2,N,J,I)=DR
                  TFTMP(N,J,I)=XNEW(N)/24.0
                ENDDO
              ELSE IF(KINDAT .GE. 2) THEN
                DO N = 1, NCUT0
                  U(N) = F_STATMP(1,N,J,I)
                ENDDO
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NCUT0,IDUMMY)
                DO N = 1, IDUMMY 
                  F_STATMP(1,N,J,I) = YNEW(N)
                  F_STATMP(2,N,J,I) = 0.0
                  TFTMP(N,J,I) = XNEW(N)/24.0
                ENDDO
              ENDIF 
            ENDDO
          ENDDO
        ENDIF

C        WRITE(*,*) 'IMAXF=',IMAXF,'NDAY=',NDAY
        DO I = 1, IMAXF
          DO J = 1, JMAX
            DO N = 1, IDUMMY 
              TF(N,J,I) = TFTMP(N,J,I)
              F(N,J,I) = F_STATMP(1,N,J,I)
              DIRF(N,J,I) = F_STATMP(2,N,J,I)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
399   FORMAT(F10.5,50F10.4)

CC  Reading persistent forecast pf(L,J,K),tf(L,J,K)
      IF(ISWITCH(5) .GT. 0) THEN 
        IF(FINCLUDE0) THEN
          NCUT0 = INT(NFDURATION/DELT_M+0.1)+1 ! fcst duration > 24h
          IDUMMY = INT(NFDURATION/DELT+0.1)+1  ! fcst duration > 24h
        ELSE
          NCUT0 = INT(NFDURATION/DELT_M+0.1)   ! fcst duration > 24h
          IDUMMY = INT(NFDURATION/DELT+0.1)    ! fcst duration > 24h
        END IF

        FILEIN = TRIM(FILESHORT)//'_persistence.dat'
        OPEN(10,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD')
        IF(KINDAT .EQ. 1) THEN
          DO I = 1, NDAY
            DO J = 1, JMAX
              DO N = 1, NCUT0
                READ(10,*,END=140)TFTMP(N,J,I),IYR0,ICM,ICD,
     1            IHR,IMN,(WLTMP(N9),N9=1,4)
                F_STATMP(1,N,J,I) = WLTMP(1)
                F_STATMP(2,N,J,I) = WLTMP(2)
              ENDDO
            ENDDO
          ENDDO 
        ELSE IF(KINDAT .GE. 2) THEN
          DO I = 1, NDAY
            DO J = 1, JMAX
              DO N = 1, NCUT0
                READ(10,*,END=140) TFTMP(N,J,I),IYR0,ICM,ICD,
     1            IHR,IMN,WLTMP(1)
                F_STATMP(1,N,J,I) = WLTMP(1)
                F_STATMP(2,N,J,I) = 0.0
              ENDDO
            ENDDO
          ENDDO 
        ENDIF 
140     CONTINUE
        CLOSE(10)

        IMAXF = I-1
        IF(ABS(DELT-DELT_M) .GT. 0.01) THEN  
          DO I = 1, IMAXF
            DO J = 1, JMAX
              TSTART0 = TFTMP(1,J,I)
              TFINISH0 = TFTMP(NCUT0,J,I)
              DO N = 1, NCUT0
                TC(N) = TFTMP(N,J,I)*24.0
              ENDDO
              IF(KINDAT .EQ. 1) THEN
                DO N = 1, NCUT0  
                  SP = F_STATMP(1,N,J,I)
                  DR = F_STATMP(2,N,J,I)
                  U(N) = SP*SIN(DR*PII)
                  V(N) = SP*COS(DR*PII)
                ENDDO 
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NCUT0,IDUMMY)
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,V,XNEW,YTMP,NCUT0,IDUMMY)
                DO N = 1, IDUMMY
                  CALL VELDIR(YTMP(N),YNEW(N),DR,SP)
                  IF(DR .GT. 360.0) DR = DR - 360.0
                  F_STATMP(1,N,J,I) = SP
                  F_STATMP(2,N,J,I) = DR
                  TFTMP(N,J,I) = XNEW(N)/24.0
                ENDDO
              ELSE IF(KINDAT .GE. 2) THEN
                DO N = 1, NCUT0
                  U(N)=F_STATMP(1,N,J,I)
                ENDDO 
                CALL EQUAL_INTERVAL1(TSTART0,TFINISH0,DELT,DELT_M,0,
     1            CRITERIA1,CRITERIA2,IMX,TC,U,XNEW,YNEW,NCUT0,IDUMMY)
                DO N = 1,IDUMMY
                  F_STATMP(1,N,J,I) = YNEW(N)
                  TFTMP(N,J,I) = XNEW(N)/24.0
                ENDDO
              ENDIF 
            ENDDO
          ENDDO
        ENDIF

        DO I = 1, IMAXF
          DO J = 1, JMAX
            DO N = 1, IDUMMY
              TF(N,J,I) = TFTMP(N,J,I)
              PF(N,J,I) = F_STATMP(1,N,J,I)
              DIRPF(N,J,I) = F_STATMP(2,N,J,I)
            ENDDO
          ENDDO
        ENDDO
        WRITE(*,*) 'The days of persistence forecast, IMAXF=',IMAXF
      ENDIF

CC  Do the datum adjustment for all model simulations 
      IF(KINDAT .NE. 1) THEN
        IF(FINCLUDE0) THEN
          NCUT0 = INT(NFDURATION/DELT_M+0.1)+1
          IDUMMY = INT(NFDURATION/DELT+0.1)+1
        ELSE
          NCUT0 = INT(NFDURATION/DELT_M+0.1)
          IDUMMY = INT(NFDURATION/DELT+0.1)
        END IF

        IF(ISWITCH(1) .GT. 0) THEN
          DO I = 1, IMAXB
            IF(ABS(S(I)) .LT. 900.0) THEN
              S(I) = S(I)+DIRFLOOD
            END IF
          ENDDO
        ENDIF

        IF(ISWITCH(2) .GT. 0) THEN
          DO I = 1, IMAXB
            IF(ABS(HI(I)) .LT. 900.0) THEN
              HI(I) = HI(I)+DIRFLOOD
            END IF
          ENDDO
        ENDIF

        IF(ISWITCH(3) .GT. 0) THEN
          DO I = 1,IMAXB
            IF(ABS(C(I)) .LT. 900.0) THEN
              C(I) = C(I)+DIRFLOOD
            END IF
          ENDDO
        ENDIF

        IF(ISWITCH(4) .GT. 0) THEN
          DO I=1,IMAXF
            DO J=1,JMAX
              DO N=1,IDUMMY
                IF(ABS(F(N,J,I)) .LT. 900.0) THEN
                  F(N,J,I) = F(N,J,I)+DIRFLOOD
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
              
      FILEIN = TRIM(FILESHORT)//'_alldata.dat'
      OPEN(10,FILE=FILEIN,FORM='FORMATTED')
      IF(KINDAT .EQ. 1) THEN
        IF(ISWITCH(4) .GT. 0) THEN
          III = 0
          DO I = 1, IMAXF
            DO J = 1, 1
              DO N = 1, IPDAY
                III = III+1
                XTMP(III) = TF(N,J,I)
                XNEW(III) = F(N,J,I)
                YNEW(III) = DIRF(N,J,I)
              ENDDO
            ENDDO
          ENDDO
          DO I = 1, IMAXB
            TMP(1,I) = -999.0
            TMP(2,I) = -999.0
            DO II = 1, III
              IF(ABS(T(I)-XTMP(II)) .LT. 0.001) THEN   
                TMP(1,I) = XNEW(II)
                TMP(2,I) = YNEW(II)
                EXIT
C                GOTO 333
              ENDIF
            ENDDO
C 333        CONTINUE          
          ENDDO
        ELSE
          DO I = 1, IMAXB
            XTMP(I) = T(I)
            TMP(1,I) = -999.0
            TMP(2,I) = -999.0
          ENDDO
        ENDIF

        DO I = 1, IMAXB
          S_M = -999.9
          DIRS_M = -999.9
          IF(ISWITCH(1) .GT. 0) THEN
	    S_M = S(I)
	    DIRS_M = DIRS(I)
	  ENDIF  
          HI_M = -999.9
	  DIRHI_M = -999.9

          IF(ISWITCH(2) .GT. 0) THEN
	    HI_M = HI(I)
   	    DIRHI_M = DIRHI(I)
	  ENDIF  

          C_M = -999.9
	  DIRC_M = -999.0
          IF(ISWITCH(3) .GT. 0) THEN
	    C_M = C(I)
	    DIRC_M = DIRC(I)
	  ENDIF  
          WRITE(10,666) T(I),O(I),DIRO(I),A(I),DIRA(I),S_M,DIRS_M,
     1	     HI_M,DIRHI_M,C_M,DIRC_M,TMP(1,I),TMP(2,I)
        ENDDO
      ELSE IF(KINDAT .GE. 2) THEN
        IF(ISWITCH(4) .GT. 0) THEN
          III = 0
          DO I = 1, IMAXF
            DO J = 1, 1
              DO N = 1, IPDAY
                III = III+1
                XTMP(III) = TF(N,J,I)
                XNEW(III) = F(N,J,I)
              ENDDO
            ENDDO
          ENDDO

          DO I = 1, IMAXB
            TMP(1,I) = -999.0
            DO II = 1, III
              IF(ABS(T(I)-XTMP(II)) .LT. 0.001) THEN   
                TMP(1,I) = XNEW(II)
                GOTO 334
              ENDIF
            ENDDO
334         CONTINUE
          ENDDO
        ELSE
          DO I = 1, IMAXB
            XTMP(I) = T(I)
            TMP(1,I) = -999.0
          ENDDO
        ENDIF
        DO I = 1, IMAXB
          S_M = -999.9
          IF(ISWITCH(1) .GT. 0) THEN
            S_M = S(I)
          END IF
          HI_M = -999.9
          IF(ISWITCH(2) .GT. 0) THEN
            HI_M = HI(I)
          END IF
          C_M = -999.9
          IF(ISWITCH(3) .GT. 0) THEN
            C_M = C(I)
          END IF
          WRITE(10,667) T(I),O(I),A(I),S_M,HI_M,C_M,TMP(1,I)
	ENDDO
      ENDIF
666   FORMAT(F13.8,6(1X,F9.4,1X,F7.1))      	    
667   FORMAT(F13.8,6(1X,F9.4))      	    
      WRITE(*,*) 'The data reading is completed!'

      RETURN
      END

