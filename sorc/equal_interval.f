C   This subroutine converts time series with interval of delt0 into 
C     a new time series with equally interval of delta.
C   If delt0 is less than delt, cubic spline interpolation is used 
C     to series convertion.
C   If there are gaps in the original series, the gaps will be filled 
C     using linear, cubic or SVD method.
C   If gap < criteria1 hours, using linear interpolated values to 
C     fill the data gap.
C   If criteria1 < gap <criteria2 hours, using Spline or SVD 
C     interpolated values to fill the data gap.
C   If gap > criteria2 hours, using -999.0 to fill the data gap.
C   Method = 0, using Cubic Spline interpolation method.
C   Method = 1, using Singular Value Decomposition method.
C   Delta= desired time interval of new time series in hours.
C   Delt0= time interval of original time series in hours.
C   Day_begin is beginning time of new series in days.
C   Day_end is the end time of the new series in days.

      SUBROUTINE EQUAL_INTERVAL(DAY_BEGIN,DAY_END,DELTA,DELT0,METHOD,
     1  CRITERIA1,CRITERIA2,TIME,WL,TIME_NEW,WL_NEW,NUM,M_NEW)

      PARAMETER(NDATA = 1000,MA = 5,NUMT = 220000)
      REAL*8 TIME,DAY_BEGIN,DAY_END,DELTA,DELT0,TIME_NEW(NUMT)
      DIMENSION TIME(NUMT),WL(NUMT),WL_NEW(NUMT),
     1  XTMP(NDATA),YTMP(NDATA),X(NDATA),Y(NDATA),A(MA),AFUNC(MA)

      WRITE(*,*) 'numt = ', NUMT
      IF(NUM .GT. NUMT .OR. M_NEW .GT. NUMT) THEN
        WRITE(*,*) 'Data number exceeds the maximum array dimension.'
        WRITE(*,*) 'Redefine numt value in equal_interval.f.'
        STOP
      ENDIF

      WRITE(*,*)'input values =',DAY_BEGIN,DAY_END,DELTA,DELT0


      HOUR_BEGIN = DAY_BEGIN*24.0
      HOUR_END = DAY_END*24.0
      IF(NUM .LE. 1) STOP
      TIME0 = HOUR_BEGIN
      Y0 = WL(1)
      TWINDOW = 10.0   ! in hours,time window the data within this period are used 
      NWIN = INT(TWINDOW/delt0)
      IF(NWIN .LT. 10) NWIN = 10
      J0 = 1
      T1 = HOUR_BEGIN
      DO I = 1, M_NEW
        DO J = J0, NUM
          IF(ABS(T1-TIME(J)) .LE. 0.001) THEN
            TIME_NEW(I) = T1
            WL_NEW(I) = WL(J)
            J0 = J
            GOTO 40
          ELSE IF(T1 .LT. TIME(1)) THEN
            TIME_NEW(I) = T1
            WL_NEW(I) = -999.0
            J0 = 1
            GOTO 40        
          ELSE IF(T1 .GT. TIME(NUM)) THEN
            TIME_NEW(I) = T1
            WL_NEW(I) = -999.0
            J0 = NUM
            GOTO 40        
          ELSE IF(T1 .LT. TIME(J) ) THEN
            J0 = J
            GOTO 10
          ENDIF
        ENDDO
10      CONTINUE

        GAP = TIME(J0)-TIME(J0-1)
        IF((J0 .GT. 1) .AND. 
     1     (TIME(J0) .LT. TIME(1)+TWINDOW)) THEN
          IF(GAP .LE. CRITERIA1) THEN
            TIMEB = TIME(J0-1)
            TIMEE = TIME(J0)
            YB = WL(J0-1)
            YE = WL(J0)
            CALL LINEAR(TIMEB,YB,TIMEE,YE,T1,Y0)
            TIME_NEW(I) = T1
            WL_NEW(I) = Y0
          ELSE
            TIME_NEW(I) = T1
            WL_NEW(I) = -999.0
          ENDIF              
        ELSE IF((TIME(J0) .GT. TIME(NUM)-TWINDOW) .AND.
     1          (J0 .LE. NUM)) THEN
          IF(GAP .LE. CRITERIA1) THEN
            TIMEB = TIME(J0-1)
            TIMEE = TIME(J0)
            YB = WL(J0-1)
            YE = WL(J0)
            CALL LINEAR(TIMEB,YB,TIMEE,YE,T1,Y0)
            TIME_NEW(I) = T1
            WL_NEW(I) = Y0
          ELSE
            TIME_NEW(I) = T1
            WL_NEW(I) = -999.0
          ENDIF              
        ELSE
          N0 = 0
          DO N = max0(1,J0-NWIN/2), min0(NUM,J0+NWIN/2+1)
            IF((TIME(N).GE. TIME(J0-1)-TWINDOW/2) .and. 
     1         (TIME(N).LE. TIME(J0)+TWINDOW/2) ) THEN
              N0 = N0+1
              X(N0) = TIME(N)
              Y(N0) = WL(N)
            ENDIF
          ENDDO
          NFIT = N0
          IF(NFIT .LT. NWIN/2) THEN
            IF(GAP .LE. CRITERIA1) THEN
              TIMEB = TIME(J0-1)
              TIMEE = TIME(J0)
              YB = WL(J0-1)
              YE = WL(J0)
              CALL LINEAR(TIMEB,YB,TIMEE,YE,T1,Y0)
              TIME_NEW(I) = T1
              WL_NEW(I) = Y0
            ELSE IF(GAP .GT. CRITERIA1) THEN
              TIME_NEW(I) = T1
              WL_NEW(I) = -999.0
            ENDIF
          ELSE IF(NFIT .GE. NWIN/2) THEN 
            IF(GAP .LE. CRITERIA1) THEN
              TIMEB = TIME(J0-1)
              TIMEE = TIME(J0)
              YB = WL(J0-1)
              YE = WL(J0)
              CALL LINEAR(TIMEB,YB,TIMEE,YE,T1,Y0)
              TIME_NEW(I) = T1
              WL_NEW(I) = Y0
            ELSEIF((GAP .GT. CRITERIA1) .AND.(GAP .LE. CRITERIA2)) THEN 
              IF(METHOD .EQ. 0) THEN
                CALL SPLINE(NFIT,X,Y,T1,YNEW)
              ELSE IF(METHOD .EQ. 1) THEN
                CALL SVD(NFIT,MA,X,Y,T1,YNEW)
              ENDIF
              TIME_NEW(I) = T1
              WL_NEW(I) = YNEW
            ELSE IF(GAP .GT. CRITERIA2) THEN
              TIME_NEW(I) = T1
              WL_NEW(I) = -999.0
            ENDIF
          ENDIF           
        ENDIF   
40      CONTINUE
        T1 = HOUR_BEGIN+I*DELTA
        IF(T1 .LT. TIME(J0)) J0 = MAX0(1,J0-1)
      ENDDO 

      RETURN
      END

