      PROGRAM REFORMAT_USGS
      CHARACTER*80 FILEINPUT,FILEOUTPUT,BUFFER
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB

      CALL GETARG(1,BUFFER)
      READ(BUFFER,*) IYRS, IMMS, IDDS, IHHS, MNS
      TIME = IHHS + MNS/60.0
      CALL GETARG(2,BUFFER)
      READ(BUFFER,*) IYRE, IMME, IDDE, IHHE, MNE
      TIMEL = IHHE + MNE/60.0
      CALL GETARG(3,BUFFER)
      READ(BUFFER,*) KINDAT
      CALL GETARG(4,FILEINPUT)
      CALL GETARG(5,FILEOUTPUT)

      YEARB  = IYRS
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)

      YEARB  = IYRS * 1.0
      MONTHB = IMMS * 1.0
      DAYB   = IDDS * 1.0
      HOURB  = IHHS * 1.0
      JDAY0 = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JBASE_DATE + 1.0

      YEARB  = IYRE * 1.0
      MONTHB = IMME * 1.0
      DAYB   = IDDE * 1.0
      HOURB  = IHHE * 1.0
      JDAY1 = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JBASE_DATE + 1.0
   
      OPEN(9,FILE=TRIM(FILEOUTPUT))
      OPEN(10,FILE=TRIM(FILEINPUT))

10    CONTINUE
      IF(KINDAT .EQ. 1) THEN
        READ(10,*,END=999) YEAR,RMONTH,DAY,HOUR,RMIN,SPD,DIR
      ELSE
        READ(10,*,END=999) YEAR,RMONTH,DAY,HOUR,RMIN,OBS
      ENDIF
      YEARB  = YEAR
      MONTHB = RMONTH
      DAYB   = DAY
      HOURB  = HOUR + RMIN/60.0
      JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JBASE_DATE + 1.0

      IF(JDAY .GE. JDAY0 .AND. JDAY .LE. JDAY1) THEN
        TDAY   = JDAY
        IYEAR = INT(YEAR)
        IMON  = INT(RMONTH)
        IDAY  = INT(DAY)
        IHR   = INT(HOUR)
        IMN   = INT(RMIN)
        IF(KINDAT .EQ. 1) THEN
          SPD = SPD * 0.51444        !!  Convert knots into m/s
          U   = SPD * SIN(DIR*3.1415926/180.)
          V   = SPD * COS(DIR*3.1415926/180.)
          WRITE(9,100) TDAY,IYEAR,IMON,IDAY,IHR,IMN,SPD,DIR,U,V
        ELSE
          IF(KINDAT .EQ. 2) OBS = OBS*0.3048 !! Convert feet to meters
  	  WRITE(9,101) TDAY,IYEAR,IMON,IDAY,IHR,IMN,OBS    
        ENDIF
        GOTO 10
      ELSE
        GOTO 10
      END IF

100   FORMAT(F10.5,I5,4I3,4F10.4)
101   FORMAT(F10.5,I5,4I3, F10.4)

999   CLOSE(9)
      CLOSE(10)

      STOP
      END

