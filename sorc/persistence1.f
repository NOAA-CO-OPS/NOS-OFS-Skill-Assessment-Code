CC    Program persistence1.f
CC    This program is used to make persistence water level forecast
CC    forecast=tide prediction + offset(observation -tide prediction, at last point)

      PARAMETER(NUM=200000)
      CHARACTER*200 FILEO,FILET,FILENAME
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      REAL*8 JBEGIN,JEND
      REAL*4  TIMEO(NUM), TIMET(NUM), OBS(NUM,2), TIDES(NUM,2)
      REAL*4, ALLOCATABLE ::  TIME(:,:,:)
      REAL*4, ALLOCATABLE ::  FORE(:,:,:,:)
      LOGICAL FEXIST

      READ(5,'(A200)') FILEO
      READ(5,'(A200)') FILET
      READ(5,'(A200)') FILENAME
      READ(5,*) KINDAT

      READ(5,*) IYRS,MMS,IDDS,IHHS,IMINS 
      YEARB  = IYRS*1.0
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JDAY0  = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      YEARB  = IYRS*1.0
      MONTHB = MMS*1.0
      DAYB   = IDDS*1.0
      HOURB  = IHHS*1.0 + IMINS/60.0
      JBEGIN = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0

      READ(5,*) IYRE,MME,IDDE,IHHE,IMINE
      YEARB  = IYRE*1.0
      MONTHB = MME*1.0
      DAYB   = IDDE*1.0
      HOURB  = IHHE*1.0 + IMINE/60.0
      JEND   = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
      NDAYS  = INT(JEND-JBEGIN+0.1) + 1
      IJS    = INT(JBEGIN)
      IJE    = INT(JEND)
      WRITE(*,*) 'BEGIN TIME=',IJS,' END TIME=',IJE,' NDAYS=',NDAYS

      READ(5,*) DELT
      DELT = DELT/60.0
      READ(5,*) NCYCLE
      READ(5,*) NFDURATION
      NCUT = INT(NFDURATION/DELT+0.1) + 1
      WRITE(*,*) 'DELT = ', DELT,'NCYCLE = ', NCYCLE,'NCUT = ', NCUT

      ALLOCATE(TIME(NCUT,NCYCLE,NDAYS))
      ALLOCATE(FORE(NCUT,NCYCLE,NDAYS,2))
      OPEN(20, FILE = FILENAME)
      INQUIRE(FILE=TRIM(FILEO), EXIST=FEXIST)
      IF(.NOT. FEXIST) THEN
        WRITE(*,*) 'filein=',TRIM(FILEO), ' does not exist '
        WRITE(*,*) 'No observation is available'
        WRITE(*,*) 'use tidal prediction as persistence forecast'
        NUMO = (JEND-JBEGIN)*24.0/DELT + 1
	DO N = 1, NUMO
	  TIMEO(N) = JBEGIN + (N-1)*DELT/24.0 - JDAY0
          OBS(N,1) = -999.9
          OBS(N,2) = -999.9
	ENDDO  
	GOTO 45
      ENDIF

      OPEN(10, FILE = FILEO)
      DO N = 1, NUM
        IF(KINDAT .EQ. 1) THEN
          READ(10,*,END=40) TIMEO(N),IYR,ICM,ICD,IHR,IMN,SP,DR
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMEO(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JDAY0
          OBS(N,1) = SP*SIN(DR*3.1415926/180.)
          OBS(N,2) = SP*COS(DR*3.1415926/180.)
        ELSE IF(KINDAT .GE. 2) THEN
          READ(10,*,END=40) TIMEO(n),IYR,ICM,ICD,IHR,IMN,obs(N,1)
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMEO(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JDAY0
        ENDIF
      ENDDO
40    CLOSE(10)

      NUMO = N - 1
      WRITE(*,*)'Number of observation is ',NUMO
      WRITE(*,*)'Fist and last times are ',TIMEO(1),TIMEO(NUMO)
45    CONTINUE
      INQUIRE(FILE=TRIM(filet),EXIST=FEXIST)
      IF(.NOT. FEXIST) THEN
           WRITE(*,*) 'filein=',trim(filet), ' does not exist '
           WRITE(*,*) 'Tidal prediction has to be existed in order to'
           WRITE(*,*) 'make persistence forecast ...'
           STOP
      ENDIF

      OPEN(10,FILE = FILET)
      DO N = 1, NUM
        IF(KINDAT .eq. 1) THEN
          READ(10,*,END=50) TIMET(N),IYR,ICM,ICD,IHR,IMN,sp,dr
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMET(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JDAY0
          TIDES(N,1) = SP*SIN(DR*3.1415926/180.0)
          TIDES(N,2) = SP*COS(DR*3.1415926/180.0)
        ELSEIF(KINDAT .GE. 2) THEN
          READ(10,*,END=50) TIMET(N),IYR,ICM,ICD,IHR,IMN,TIDES(N,1)
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMET(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB) - JDAY0
        ENDIF
      ENDDO

50    CLOSE(10)
      NUMT = N - 1
      WRITE(*,*) 'The number of tide prediction is ',NUMT
      WRITE(*,*) 'Fist and last time ',timet(1),timet(NUMT)
      DO I = 1, NDAYS
        DO J = 1, NCYCLE
          DO N = 1, NCUT
            TIME(N,J,I) = IJS+(I-1)+(J-1)*1.0/NCYCLE+(N-1)*DELT/24.0
            IF(N .EQ. 1) THEN
              DO IT = 1, NUMT
                IF(ABS(TIME(N,J,I)-TIMET(IT)) .LT. 0.001 ) THEN
                  IT0 = IT
                  GOTO 433
                ELSE
                  IF(IT .GE. NUMT) THEN
		    WRITE(*,*) 'TIME=',TIME(N,J,I),I,J,N
                    WRITE(6,*) 'No prediction time match with forecast'
                    STOP
                  ENDIF 
                ENDIF
              ENDDO 

433           CONTINUE
              DIFF1 = 0.0
              DIFF2 = 0.0
              DO IO=1,NUMO
                IF(ABS(TIMET(IT0)-TIMEO(IO)) .LT. 0.001) THEN
                  WRITE(*,"('N,OI,OT=',5I6)") I,J,N,IO,IT
                  WRITE(*,*) TIME(N,J,I),TIMEO(IO),TIMET(IT)
                  IF(OBS(IO,1) .GT. -900) THEN
                    IF(KINDAT .EQ. 1) THEN
                      DIFF1 = OBS(IO,1)-TIDES(IT0,1)                     
                      DIFF2 = OBS(IO,2)-TIDES(IT0,2)                     
                    ELSE IF(KINDAT .GE. 2)T HEN
                      DIFF1 = OBS(IO,1)-TIDES(IT0,1)                     
                    ENDIF
                  ELSE
                    DIFF1 = 0.0
                    DIFF2=0.0
                  ENDIF 
                  GOTO 110
                ELSE
                  IF(IO .GE. NUMO) THEN
                    WRITE(*,*) 'No observation data at', TIMET(IT0)
                    DIFF1 = 0.0
                    DIFF2 = 0.0
                  ENDIF 
                ENDIF
              ENDDO
            ENDIF

110         CONTINUE
            IF(KINDAT .EQ. 1) THEN
              FORE(N,J,I,1) = TIDES(IT0,1) + DIFF1
              FORE(N,J,I,2) = TIDES(IT0,2) + DIFF2
            ELSE IF (KINDAT .GE. 2)THEN
              FORE(N,J,I,1) = TIDES(IT0,1) + DIFF1
            ENDIF
            IF(N .EQ. 1) WRITE(*,399) TIME(N,J,I),FORE(N,J,I,1)
            IT0 = IT0 + 1
          ENDDO
        ENDDO
      ENDDO
399   FORMAT(F10.5,50F10.4)

      IYR = IYRS
      DO I = 1, NDAYS
        DO J = 1, NCYCLE
          DO N = 1, NCUT
            JDAY = TIME(N,J,I) + JDAY0
            CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
            IYEAR = INT(YEARB+0.001)
            ICM   = INT(MONTHB+0.001)
            ICD   = INT(DAYB+0.001)
            IHR   = INT(HOURB+0.001)
            IMN   = INT((HOURB-IHR)*60.0+0.1)
            DAYJ  = JDAY - JDAY0 + 1.
            IF(KINDAT .eq. 1) THEN
              u = fore(N,J,I,1)
              v = fore(N,J,I,2)
              CALL VELDIR (v,u,AANGLE,AVEL)
              IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
              sp = AVEL
              dr = AANGLE
              WRITE(20,100)dayj,IYEAR,ICM,ICD,IHR,IMN,sp,dr,u,v
            ELSE IF(KINDAT .GE. 2) THEN
              WRITE(20,100)dayj,IYEAR,ICM,ICD,IHR,IMN,fore(N,J,I,1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
100   FORMAT(F10.5,I5,4I3,10F10.4)
      CLOSE(20)

      STOP
      END

