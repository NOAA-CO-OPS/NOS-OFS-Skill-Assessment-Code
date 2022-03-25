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

      PROGRAM SKILLS
      include 'skills.inc'

      CHARACTER*200 FILESHORT,STATION_FILE,STAID,STATIONNAME,DATADIR
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      DIMENSION XTMP(IMX),YTMP(IMX),XNEW(IMX),YNEW(IMX)
      DIMENSION ZHALL(NMX),THALL(NMX),ZLALL(NMX),TLALL(NMX)
      CHARACTER*200 FILEIN
      FILEIN='8638901W_HS.obs'
      FILEIN='8638901W_HS_modeltides.dat.noDA'
      OPEn(1,file=FILEIN)
C  Read all data
      DO I = 1, IMX
       read(1,*,end=100)tc(i),j,j,j,j,j,a(i)
      ENDDO
100   CONTINUE
      CLOSE(1)
      NUM=I-1
      print *,'NUM=',num
      FILEIN='8638901W_HS_modeltides.dat.DA'
      OPEn(1,file=FILEIN)
C  Read all data
      DO I = 1, IMX
       read(1,*,end=110)xtmp(i),j,j,j,j,j,s(i)
       if( (xtmp(i) - tc(NUM)) .gt. 0.02)goto 110
       if (abs(tc(i)-xtmp(i)) .gt. 0.02)then
         print *,'I= ',I,tc(i),xtmp(i)
       endif
      ENDDO
110   CONTINUE
      NUM1=I-1     
      print *,'NUM1=',num1
      FILEIN='8638901W_HS.obs'
      OPEn(1,file=FILEIN)
      read(1,*,end=100)tc11,j,j,j,j,j,a111
C  Read all data
      DO I = 1, IMX
       read(1,*,end=120)xtmp(i),j,j,j,j,j,o(i)
       if( (xtmp(i) - tc(NUM)) .gt. 0.02)goto 120
      ENDDO
120   CONTINUE
      CLOSE(1)
      NUM2=I-1
      print *,'NUM2=',num2
      avg1=0.0
      avg2=0.0
      avg3=0.0
      DO I=1,NUM
        avg1=avg1+o(i)
        avg2=avg2+a(i)
        avg3=avg3+s(i)
        write(11,'(4f12.5)')tc(i),o(i),a(i),s(i)
      ENDDO
      avg1=avg1/num
      avg2=avg2/num
      avg3=avg3/num
      print *,'mean=',avg1,avg2,avg3
      STOP


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


500   FORMAT('MHW, MLW, and New error criteria X1= ',3F8.2,' at',a10)
999   CONTINUE       

      STOP
      END

!            SM0 = SM(B,IMX,IMAX)
!            RMSE0 = RMSE(B,IMX,IMAX)
!            SD0 = SD(B,IMX,IMAX)
!           VOF0 = VOF(B,IMX,IMAX,2.0*ERR0)
!            CF0 = CF(B,IMX,IMAX,ERR0)
!            POF0 = POF(B,IMX,IMAX,2.0*ERR0)




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





 





