C  Modified by Lianyuan Zheng on 03/01/2017

      PROGRAM REFWL
      PARAMETER (NMAX = 500000)
      CHARACTER*80 FILEINPUT,FILEOUTPUT,BUFFER
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      REAL*8 TDAY,TDAY0,TMP1(NMAX)
      REAL TMP2(NMAX),ONED1(NMAX)

      CALL GETARG(1,BUFFER)
      READ(BUFFER,*) IYRS, IMMS, IDDS, IHHS, MNS
      TIME = IHHS + MNS/60.0
      CALL GETARG(2,BUFFER)
      READ(BUFFER,*) IYRE, IMME, IDDE, IHHE, MNE
      TIMEL = IHHE + MNE/60.0
      CALL GETARG(3,FILEINPUT)
      CALL GETARG(4,FILEOUTPUT)
      YEARB = IYRS
      MONTHB = 1.0
      DAYB = 1.0
      HOURB = 0.0
      JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      TDAY0 = -99999.0   

      OPEN(9,FILE=TRIM(FILEOUTPUT))
      OPEN(10,FILE=TRIM(FILEINPUT))
      N = 0
10    READ(10,*,END=19) YEAR, RMONTH, DAY, HOUR, RMIN, OBS
      IF(OBS .LE. -99.9) GOTO 10
      YEARB  = YEAR
      MONTHB = RMONTH
      DAYB   = DAY
      HOURB = HOUR + RMIN/60.0
      JDAY = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      TDAY = JDAY - JBASE_DATE
      IYEAR = INT(YEAR)
      IMON  = INT(RMONTH)
      IDAY  = INT(DAY)
      IHR   = INT(HOUR)
      IMN   = INT(RMIN)
      IF(TDAY .GT. TDAY0) THEN
        N = N+1
	TMP1(N) = TDAY
	TMP2(N) = OBS
        TDAY0 = TDAY
      ENDIF
      GOTO 10

19    NTOL = N
      IF(NTOL .GE. 10) THEN
        PRINT*,' The number of the observation data NTOL = ',NTOL
C  Calculate Mean Value
        AVG = 0.0
        NTMP = 0
        DO N = 1, NTOL
          NTMP = NTMP + 1
          AVG  = AVG + TMP2(N)
          ONED1(NTMP) = TMP2(N)
        ENDDO
        IF(NTMP .GT. 0) AVG = AVG/NTMP

C  Calculate Standard Deviation Value
        SD = 0.0
        DO N = 1, NTMP
          SD = SD + (ONED1(N)-AVG)**2
        ENDDO
        SD = SQRT(SD/(NTMP-1))
        BOUND_L = AVG-5.0*SD
        BOUND_U = AVG+5.0*SD
        WRITE(*,*) ' Lower and Upper limits and Average =',
     1     BOUND_L,BOUND_U,AVG

        DO N = 1, NTOL
          IF((TMP2(N).GE.BOUND_L).AND.(TMP2(N).LE.BOUND_U)) THEN
            JDAY = TMP1(N) + JBASE_DATE
            CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
            IYEAR = INT(YEARB)
            ICM   = INT(MONTHB+0.001)
            ICD   = INT(DAYB+0.001)
            IHR   = INT(HOURB+0.001)
            IMN   = INT((HOURB-IHR)*60+0.1)
            TDAY = TMP1(N) + 1
            WRITE(9,100) TDAY, IYEAR, ICM, ICD, IHR, IMN, TMP2(N)
          ENDIF
        ENDDO
      ELSE
        PRINT*,' The number of the observation data is < 10!!!'
      END IF

100   FORMAT(F13.8,I5,4I3,F10.4)

999   STOP
      END

