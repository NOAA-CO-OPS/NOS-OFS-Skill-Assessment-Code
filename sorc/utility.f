      SUBROUTINE VELDIR(AAVN,AAVE,AANGLE,AVEL)
      REAL*8 RADDEG
      DATA RADDEG /57.295779513802D0/

C      IF(AAVN .EQ. 0.0) THEN
C        IF(AAVE .LT. 0.0) THEN
C          AANGLE = 270.0
C        ELSE IF(AAVE .EQ. 0.0) THEN
C          AANGLE =   0.0
C        ELSE
C          AANGLE =  90.0
C        END IF
C      ELSE
C        AANGLE = ATAN(AAVE/AAVN)
C        AANGLE = AANGLE * RADDEG
C        IF(AAVN .LE. 0.0) THEN
C          AANGLE = AANGLE + 180.0
C        ELSE
C          IF(AAVE .LE. 0.0) THEN
C            AANGLE = AANGLE + 360.0
C          END IF
C        END IF
C      END IF
C      AVEL = SQRT(AAVE**2 + AAVN**2)
C      AANGLE = MOD(AANGLE + 360.0,360.0)

C  NOTE: The AANGLE CAN BE DETERMINED AS FOLLOWING:
      AVEL = SQRT(AAVE**2 + AAVN**2)
      IF(AVEL .EQ. 0.0) THEN
        AANGLE = 0.0
      ELSE
        AANGLE = ATAN2(AAVN, AAVE) * RADDEG
        AANGLE = 90.0 - AANGLE
      END IF
      IF(AANGLE .LT. 0.0) AANGLE = 360.0 + AANGLE

      RETURN
      END


C  Returns the last non blank character position in a string
      SUBROUTINE NCRGHT(LINE,NC)
      CHARACTER*(*) LINE

      ILIM = LEN (LINE)
      DO 100 I = 1,ILIM
        IF(LINE(I:I) .NE. ' ') THEN
          NC = I
        END IF
100   CONTINUE

      RETURN
      END


      SUBROUTINE GREGORIAN(JDAY,YR,MONTH,DAY,HOUR)
      REAL*8 JDAY,YR,MONTH,DAY,HOUR
      INTEGER NDP(13)
      INTEGER NDM(12)

      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/

      KD   = INT(JDAY)
      HOUR = (JDAY - KD) * 24.0
      
C  KD = 1 CORRESPONDS TO JANUARY 1, 0000
C  GIVEN THE (GREGORIAN) DAY#, KD, AS CALCULATED ABOVE IN THIS ROUTINE,
C  ENTRY DMY RETURNS THE (GREGORIAN) DAY, MONTH, YEAR AND CENTURY.
      IF(KD .LE. 0) THEN
        WRITE(11,5040) KD
        RETURN
      END IF
      KKD = KD

C  CALCULATE ICC AND SUBTRACT THE NUMBER OF DAYS REPRESENTED BY ICC
C  JFH IS THE NUMBER OF 400 YEAR INTERVALS UP TO KKD
C  JCC IS THE NUMBER OF ADDITIONAL CENTURIES UP TO KKD
      JFH = KKD/146097
      KKD = KKD - JFH*146097
      IF(KKD .LT. 36525) THEN
        JCC = 0
      ELSE
        KKD = KKD - 36525
        JCC = 1 + KKD/36524
        KKD = KKD - (JCC-1) * 36524
      END IF

      ICC = 4*JFH + JCC
      IF(KKD .EQ. 0) THEN
        ICC = ICC-1
        IYY = 99
        IMM = 12
        IDD = 31
        YR  = ICC*100 + IYY
	MONTH = IMM
	DAY = IDD
	RETURN
      ENDIF

C  CALCULATE IYY. JFY IS THE NUMBER OF FOUR YEAR INTERVALS IN THE
C  CURRENT CENTURY. THE FIRST FOUR YEAR INTERVAL IS SHORT (1460
C  DAYS RATHER THAN 1461)IF THE CURRENT CENTURY IS NOT DIVISIBLE 
C  BY 4, AND IN THIS CASE JCC.NE.0 AS CALCULATED ABOVE.
C  CALCULATE JFY:
      JFY = 0
      IF((JCC .NE. 0) .AND. (KKD .GE. 1460)) THEN
        JFY = 1
        KKD = KKD - 1460
      END IF
      KK  = KKD/1461
      JFY = JFY + KK
      KKD = KKD - KK*1461

C  CALCULATE JYY, THE REMAINING YEARS OF THE CURRENT CENTURY UP
C  TO THE CURRENT DAY:
C  THE NEXT YEAR IS NOT A LEAP YEAR IF JFY=0 AND JCC.NE.0.
      JYY = 0
      IF(JFY .EQ. 0 .AND. JCC .NE. 0) GOTO 20
      IF(KKD .LT. 366) GOTO 30
      JYY = 1
      KKD = KKD - 366
20    JYYY = KKD / 365
      JYY = JYY + JYYY
      KKD = KKD - JYYY * 365
30    IYY = 4 * JFY + JYY
      IF(KKD .EQ. 0) THEN
        IYY = IYY - 1
        IMM = 12
        IDD = 31
        YR  = ICC * 100 + IYY
        MONTH = IMM
        DAY = IDD
        RETURN
      END IF

C  SET L=1 IF WE HAVE A LEAP YEAR.
      L = 0
      IF(IYY-(IYY/4)*4 .NE. 0) GOTO 40
      IF(IYY .EQ. 0 .AND. (ICC-(ICC/4)*4) .NE. 0) GOTO 40
      L = 1

C  CALCULATE IMM AND IDD
40    IF(KKD .GT. 31) GOTO 50
      IMM = 1
      IDD = KKD
      YR = ICC * 100 + IYY
      month = IMM
      day = IDD
      RETURN

50    IF(KKD .GT. 59) GOTO 60
      IMM = 2
      IDD = KKD - 31
      YR  =ICC * 100 + IYY
      MONTH = IMM
      DAY = IDD
      RETURN

60    IF(KKD .GT. 60) GOTO 70
      IF(L .EQ. 0) GOTO 70
      IMM = 2
      IDD = 29
      YR  = ICC * 100 + IYY
      MONTH = IMM
      DAY = IDD
      RETURN

70    IF(L .EQ. 1) KKD = KKD - 1
      DO 80 I = 4, 13
        IF(KKD .GT. NDP(I)) GOTO 80
        IMM = I - 1
        IDD = KKD - NDP(I-1)
        YR  = ICC * 100 + IYY
        MONTH = IMM
        DAY = IDD
        RETURN
80    CONTINUE
90    WRITE(11,5050)
5050  FORMAT(' ERROR IN GREGORIAN')
5040  FORMAT(' KD = ',I7,'  INVALID INPUT. DMY STOP.')

      RETURN
      END


C  GIVEN DAY,MONTH,YEAR AND CENTURY(EACH 2 DIGITS), GDAY RETURNS
C  THE DAY#, KD BASED ON THE GREGORIAN CALENDAR.
C  THE GREGORIAN CALENDAR, CURRENTLY 'UNIVERSALLY' IN USE WAS
C  INITIATED IN EUROPE IN THE SIXTEENTH CENTURY. NOTE THAT GDAY
C  IS VALID ONLY FOR GREGORIAN CALENDAR DATES.
C
C  KD=1 CORRESPONDS TO JANUARY 1, 0000 
C	
C  NOTE THAT THE GREGORIAN REFORM OF THE JULIAN CALENDAR 
C  OMITTED 10 DAYS IN 1582 IN ORDER TO RESTORE THE DATE
C  OF THE VERNAL EQUINOX TO MARCH 21 (THE DAY AFTER
C  OCT 4, 1582 BECAME OCT 15, 1582), AND REVISED THE LEAP 
C  YEAR RULE SO THAT CENTURIAL YEARS NOT DIVISIBLE BY 400
C  WERE NOT LEAP YEARS.
C
C  THIS ROUTINE WAS WRITTEN BY EUGENE NEUFELD AT IOS, IN JUNE 1990.
      FUNCTION JULIAN(YR,MONTHB,DAYB,HOURB)

      REAL*8 YR, MONTHB, DAYB, HOURB
      REAL*8 YB, MB, JULIAN
      INTEGER NDP(13)
      INTEGER NDM(12)

      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/

      IDD = INT(DAYB)
      IMM = INT(MONTHB)
      ICC = INT(YR/100)
      IYY = INT(YR - ICC * 100)

      LP = 6
      IF(ICC .LT. 0) THEN
        WRITE(11,5000) ICC
        RETURN
      ENDIF

      IF(IYY .LT. 0 .OR. IYY .GT. 99) THEN
        WRITE(11,5010) IYY
        RETURN
      ENDIF

      IF(IMM .LE. 0 .OR. IMM .GT. 12) THEN
        WRITE(11,5020) IMM
        RETURN
      ENDIF

      IF(IDD .LE. 0) THEN
        WRITE(11,5030) IDD
        RETURN
      ENDIF

      IF(IMM .NE. 2 .AND. IDD .GT. NDM(IMM)) THEN
        WRITE(11,5030) IDD
        STOP
      ENDIF

      IF(IMM .EQ. 2 .AND. IDD .GT. 29) THEN
        WRITE(11,5030) IDD
        STOP
      ENDIF

      IF(IMM .EQ. 2 .AND. IDD .GT. 28 .AND. ((IYY/4)*4-IYY .NE. 0
     .   .OR. (IYY .EQ. 0 .AND. (ICC/4)*4-ICC .NE. 0))) THEN
        WRITE(11,5030) IDD
        STOP
      ENDIF

5000  FORMAT(' INPUT ERROR. ICC = ',I7)
5010  FORMAT(' INPUT ERROR. IYY = ',I7)
5020  FORMAT(' INPUT ERROR. IMM = ',I7)
5030  FORMAT(' INPUT ERROR. IDD = ',I7)

C  CALCULATE DAY# OF LAST DAY OF LAST CENTURY:
      KD = ICC*36524 + (ICC+3)/4

C  CALCULATE DAY# OF LAST DAY OF LAST YEAR:
      KD = KD + IYY*365 + (IYY+3)/4

C  ADJUST FOR CENTURY RULE:
C  (VIZ. NO LEAP-YEARS ON CENTURYS EXCEPT WHEN THE 2-DIGIT
C  CENTURY IS DIVISIBLE BY 4.)
      IF(IYY.GT.0.AND.(ICC-(ICC/4)*4).NE.0) KD = KD - 1

C  KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF LAST YEAR.
C  CALCULATE DAY# OF LAST DAY OF LAST MONTH:
      KD = KD + NDP(IMM)

C  ADJUST FOR LEAP YEARS:
      IF(IMM.GT.2.AND.((IYY/4)*4-IYY).EQ.0.AND.((IYY.NE.0).OR.
     .   (((ICC/4)*4-ICC).EQ.0))) KD = KD + 1

C  KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF THE LAST
C  MONTH.  CALCULATE THE CURRENT DAY#:
      KD = KD + IDD
      JULIAN = dble(KD*1.0) + dble(HOURB/24.0)

      RETURN
      END


      SUBROUTINE LINEAR(X1,Y1,X2,Y2,X,Y)

      SLOPE = (Y2-Y1)/(X2-X1)
      Y     = Y1 + SLOPE * (X-X1)

      RETURN
      END


      SUBROUTINE LINEAR_DP(X1,Y1,X2,Y2,X,Y)

      REAL*8 X1,X2,X,SLOPE,YTMP
      SLOPE = dble((Y2-Y1))/(X2-X1)
      YTMP  = dble(Y1) + SLOPE * (X-X1)
      Y     = Real(YTMP)

      RETURN
      END


C  NOTE: THIS SUBROUTINE IS WORKING ONLY FOR X1 INCREASES
C        IF X1 DECREASES, IT NEEDS TO REVERSE THE ARRAY.
      SUBROUTINE LINEARARRAY(N1,X,Y,N2,X1,Y1)
      DIMENSION X(N1),Y(N1),X1(N2),Y1(N2)
      DIMENSION TPMX(N2),TPMY(N2)

      IF(X1(N2) .GE. X1(1)) THEN
        DO I = 1, N2
          TPMX(I) = X1(I)
          TPMY(I) = Y1(I)
        END DO
      ELSE
        DO I = 1, N2
          K = N2 -I + 1
          TPMX(I) = X1(K)
          TPMY(I) = Y1(K)
        END DO
      END IF

      DO I = 1, N1	 
        IF(X(I) .LE. TPMX(1)) THEN
          Y(I) = TPMY(1)
        ELSE IF(X(I) .GE. TPMX(N2)) THEN
	  Y(I) = TPMY(N2)
	ELSE
	  DO N = 1, N2-1
            IF((X(I) .GT. TPMX(N)) .AND. 
     1	       (X(I) .LE. TPMX(N+1))) THEN
              N0 = N
              GOTO 10
	    ENDIF
	  ENDDO
10	  CONTINUE  	     
          SLOPE = (TPMY(N0+1)-TPMY(N0))/(TPMX(N0+1)-TPMX(N0) )
          Y(I) = TPMY(N0) + SLOPE*(X(I)-TPMX(N0))
	ENDIF
      ENDDO	   

      RETURN
      END


      SUBROUTINE SPLINE(NUM,X,Y,XNEW,YNEW)
      DIMENSION X(NUM),Y(NUM),Y2(NUM),U(NUM)

      YP1 = (Y(2)-Y(1))/(X(2)-X(1) )
      YPN = (Y(NUM)-Y(NUM-1))/(X(NUM)-X(NUM-1) ) 
C  IBC = 1 FOR NATURAL CUBIC-SPLINE BOUNDARY CONDITION, WHICH HAS
C          ZERO SECOND DERIVATIVE ON ONE OR BOTH END POINTS
C  IBC = 2 SET EITHER OF Y2(1) OR Y2(N) TO VALUES CALCULATED FROM
C          EQUATION (3.3.5) TO MAKE THE FIRST DERIVATIVE OF THE
C          INTERPOLATING FUNCTION HAs A SPECIFIED VALUE ON EITHER
C          OR BOTH END POINTS.
      IBC = 1
      IF(IBC .EQ. 1) THEN
        Y2(1) = 0.0
        U(1)  = 0.0
      ELSE
        Y2(1) = -0.5
        U(1)  = (3.0/(X(2)-X(1))) * ((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

      DO I = 2, NUM - 1
        SIG = (X(I)-X(I-1)) / (X(I+1)-X(I-1))
        P   = SIG*Y2(I-1) + 2.0
        Y2(I) = (SIG-1.0) / P
        U(I) = (6.0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/ 
     1         (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO

      IF (IBC .EQ. 1) THEN
        QN = 0.0
        UN = 0.0
      ELSE
        QN = 0.5
        UN = (3.0/(X(NUM)-X(NUM-1)))*(YPN-(Y(NUM)-Y(NUM-1))/
     1       (X(NUM)-X(NUM-1)))
      END IF

      Y2(NUM) = (UN-QN*U(NUM-1))/(QN*Y2(NUM-1)+1.0)
      DO K = NUM-1,1,-1
        Y2(K) = Y2(K)*Y2(K+1)+U(K)
      END DO

      XA  = XNEW
      KLO = 1
      KHI = NUM
1     IF(KHI-KLO .GT. 1) THEN
        K = (KHI+KLO)/2
        IF(X(K) .GT. XA) THEN
          KHI = K
        ELSE
          KLO = K
        END IF
        GOTO 1
      END IF

      H = X(KHI) - X(KLO)
      IF(ABS(H) .LT. 1.0E-5) THEN
        PRINT *, 'Bad XA input in SPLINE'
        STOP
      END IF

      A = (X(KHI)-XA)/H
      B = (XA-X(KLO))/H
      YNEW = A*Y(KLO)+B*Y(KHI)+((A**3-A)*Y2(KLO)+
     1       (B**3-B)*Y2(KHI))*(H**2)/6.0 

      RETURN
      END


C  THIS PROGRAM SEARCHES MAXIMUM CONTINUOUS TIME DURATION IN A 
C  TIME SERIES
      SUBROUTINE CONTINUOUS(Z,IMX,GAP,IS,IFRST,ILAST)
      DIMENSION Z(IMX),M(IMX),IFRST(IMX),ILAST(IMX)

      M(1) = 1
      IPR = 0
      DO I = 2, IMX
        GP = Z(I) - Z(I-1)
        M(I) = 1
        IF(GP .GE. GAP) M(I) = 0     
      ENDDO

      IS = 0
      DO I = 1, IMX - 1
        IF(I .EQ. 1 .AND. M(I) .EQ. 1) THEN
          IS = IS + 1
          IFRST(IS) = I
        ENDIF

        IF(I .GE. 1 .AND. M(I) .EQ. 0 .AND. M(I+1) .EQ. 1) THEN
          IS = IS + 1
          IFRST(IS) = I + 1
        ENDIF
      ENDDO

      IF(IPR .EQ. 1) THEN
        WRITE(*,*) ' IS = ', IS
        IF(IS .GT. 0) THEN
          DO I = 1, IS
            WRITE(*,*) ' I = ', I,' IFRST = ', IFRST(I)
          ENDDO
        ENDIF
      ENDIF

C  LOOK FOR END OF A SERIES OF '1'S.  ILAST IS I OF LAST '1'
      IE = 0
      DO I = 1, IMX - 1
        IF(M(I) .EQ. 1 .AND. M(I+1) .EQ. 0) THEN
          IE = IE + 1
          ILAST(IE) = I
        ENDIF
      ENDDO

      IF(M(IMX) .EQ. 1) THEN
        IE = IE + 1
        ILAST(IE) = I
      ENDIF

      IF(IPR .EQ. 1) THEN
        WRITE(*,*) ' IE= ', IE
        IF(IE .GT. 0) THEN
          DO I = 1, IE
            WRITE(*,*) ' I= ', I,' ILAST= ', ILAST(I)
          ENDDO
        ENDIF
      ENDIF

C        FIND LONGEST DURATION
      IF(IS .NE. IE) THEN
        WRITE(*,"(/,3X,'**FUNCTION TMDO. STARTS NOT EQUAL TO ENDS**')")
        WRITE(*,"(5X,'IS= ',I4,' IE= ',I4)") IS, IE
        IF(IS .GT. 0) THEN
          DO I = 1, IS
            WRITE(*,*) ' I= ', I,' IFRST= ', IFRST(I)
          ENDDO
        ENDIF
        IF(IE .GT. 0) THEN
          DO I = 1, IE
            WRITE(*,*) ' I= ', I,' ILAST= ', ILAST(I)
          ENDDO
        ENDIF

        DO I = 1, IMX
          WRITE(*,"(' I=',I4,' M=',I2)") I, M(I)
        ENDDO
        STOP
      ELSE IF(IS .GT. 0) THEN
        MDO = 0
        DO I = 1, IS
          NDIF = ILAST(I) - IFRST(I)
          IF(NDIF .GT. MDO) THEN
            MDO = NDIF
            IMAX = I
          ENDIF  
        ENDDO
      ENDIF

C      ISTART = IFRST(IMAX)
C      IEND   = ILAST(IMAX)

      RETURN
      END


C  FUNCTION:
C  Calculate the arc length for two points on the spherical plane
C  input:
C     XSG,XSL,XLATL,XLONGL: Latitude and Longitude of two points
C  output:
C     D:  arc length of two points in spherical plane in meters
      SUBROUTINE DIST(XSG,XSL,XLATL,XLONGL,D)

      REAL LA0,LA1,LO0,LO1,LB,LL
      DATA CONV/57.29578/

      HAV(X) = (SIN(0.5*X))**2
      AHAV(X) = 2.0*ASIN(SQRT(X))
      LA0 = XLATL / CONV
      LO0 = XLONGL / CONV
      LA1 = XSG / CONV
      LO1 = XSL / CONV
      LL = LA0 + LA1
      LB = LA0 - LA1
      R = AHAV(HAV(LB)+COS(LA0)*COS(LA1)*HAV(LO0-LO1))
      IF(R .EQ. 0.0) THEN
        D = 0.0
      ELSE
        D = 3437.7468*R   ! in nautical miles = 1.852 km
        D = 6371.0E03*R   ! in meters
      END IF

      RETURN
      END


      SUBROUTINE SIGMA2Z_ROMS(CW_1,CW_2,H_CUT,C_CUT,CFF_M,CFF_C,
     1   HC,THETAS,SIGMA,H,ELE,KB,ZSIGMA)
      DIMENSION SIGMA(KB),ZSIGMA(KB)

      THETAB = CFF_M*TANH(C_CUT*(H-H_CUT))+CFF_C
      DO K = 1, KB
        PTHETA = SINH(THETAS*SIGMA(K))/SINH(THETAS)
        RTHETA = TANH(THETAS*(SIGMA(K)+0.5))/(2.0*TANH(THETAS*0.5))-0.5
        CSIGMA = CW_1*((1.0-THETAB)*PTHETA+THETAB*RTHETA) + 
     1           CW_2*ABS((1.0-(1.0+SIGMA(K))**(2.0*THETAS))**THETAB)
        Z0 = HC*SIGMA(K) + (H-HC)*CSIGMA
        ZSIGMA(K) = Z0 + (1.0+Z0/H)*ELE
        IF(ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
      ENDDO

      RETURN
      END


      SUBROUTINE SIGMA2Z_POM(SIGMA,H,ELE,KB,ZSIGMA)
      DIMENSION SIGMA(KB),ZSIGMA(KB)

      DO K = 1, KB
        ZSIGMA(K) = SIGMA(K)*(H+ELE) + ELE
        IF(ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
      ENDDO

      RETURN
      END


      SUBROUTINE SIGMA2Z_ROMS_FIX(SIGMA,H,ELE,KB,ZSIGMA,
     1    HC,THETAS,THETAB,TC)
      DIMENSION SIGMA(KB),ZSIGMA(KB)

      DO K = 1, KB
        PTHETA = SINH(THETAS*SIGMA(K))/SINH(THETAS)
        RTHETA = TANH(THETAS*(SIGMA(K)+0.5))/(2.0*TANH(THETAS*0.5))-0.5
        CSIGMA = (1.0-THETAB)*PTHETA+THETAB*RTHETA
        Z0 = HC*SIGMA(K)+(H-HC)*CSIGMA
        ZSIGMA(K) = Z0+(1.0+Z0/H)*ELE
        IF(ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
      ENDDO

      RETURN
      END


      SUBROUTINE SIGMA2Z_SELFE(SIGMA,H,ELE,KB,ZSIGMA,
     1    H_S,H_C,THETAS,THETAB)
      DIMENSION SIGMA(KB),ZSIGMA(KB),CS(KB),DCS(KB)

      NVRT = KB
      KZ = 1
      NSIG = NVRT-KZ+1 !# of S levels (including "bottom" & f.s.)
      S_CON1 = SINH(THETAS)

!     Compute C(s) and C'(s)
      DO K = 1, NSIG
        CS(K) = (1-THETAB)*SINH(THETAS*SIGMA(K))/SINH(THETAS)+
     1     THETAB*(TANH(THETAS*(SIGMA(K)+0.5))-TANH(THETAS*0.5))/
     2     2/TANH(THETAS*0.5)
        DCS(K) = (1-THETAB)*THETAS*COSH(THETAS*SIGMA(K))/
     1     SINH(THETAS)+THETAB*THETAS/2/TANH(THETAS*0.5)/
     2     COSH(THETAS*(SIGMA(K)+0.5))**2
      ENDDO 

      DO K = KZ, NVRT
        KIN = K-KZ+1
        HMOD2 = AMIN1(H,H_S)
        IF(HMOD2 .LE. H_C) THEN
          ZSIGMA(K) = SIGMA(KIN)*(HMOD2+ELE)+ELE
        ELSE
          ZSIGMA(K) = ELE*(1+SIGMA(KIN))+H_C*SIGMA(KIN)+
     1       (HMOD2-H_C)*CS(KIN)
        ENDIF
        IF(ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
      ENDDO

      RETURN
      END


      SUBROUTINE ONE2TWOD(A,B,IM,JM)
      DIMENSION A(IM*JM),B(IM,JM)

      DO J = 1, JM
        DO I = 1, IM
          N = (J-1)*IM + I
          B(I,J) = A(N)
        ENDDO
      ENDDO    	       

      RETURN
      END


      SUBROUTINE TWO2ONED(A,B,IM,JM)
      DIMENSION A(IM*JM),B(IM,JM)

      DO J = 1, JM
        DO I = 1, IM
          N = (J-1)*IM + I
          A(N) = B(I,J)
        ENDDO
      ENDDO    	       

      RETURN
      END

