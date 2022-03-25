CC    Program Persistence.f
CC    This program is used to make persistence water level forecast
CC    forecast=tide prediction + offset(observation -tide prediction, at last point)
      PARAMETER(NUM = 200000)
      CHARACTER*200 FILEO,FILET,FILENAME,FILEF
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      REAL*8 JBEGIN,JEND
      REAL*4 TIMEO(NUM)
      REAL*4 TIMET(NUM)
      REAL*4 OBS(NUM,2)
      REAL*4 TIDES(NUM,2)
      REAL*4, ALLOCATABLE ::  TIME(:,:,:)
      REAL*4, ALLOCATABLE ::  FORE(:,:,:,:)
      LOGICAL FEXIST

      READ(5,'(A200)') FILEO
      READ(5,'(A200)') FILET
      READ(5,'(A200)') FILEF
      READ(5,'(A200)') FILENAME
      READ(5,*) KINDAT
      READ(5,*) IYRS, MMS, IDDS, IHHS, IMINS 
      READ(5,*) IYRE, MME, IDDE, IHHE, IMINE
      READ(5,*) DELT
      DELT = DELT/60.0
      READ(5,*) NCYCLE
      READ(5,*) NFDURATION

      NCUT = INT(NFDURATION/DELT+0.1)+1
      WRITE(*,*) 'DELT = ',DELT,' NCYCLE = ',NCYCLE,' NCUT = ',NCUT
      YEARB  = IYRS*1.0
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JDAY0 = JULIAN(YEARB,MONTHB,DAYB,HOURB)

      YEARB  = IYRS*1.0
      MONTHB = MMS*1.0
      DAYB   = IDDS*1.0
      HOURB  = IHHS*1.0 + IMINS/60.0
      JBEGIN = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0

      YEARB  = IYRE*1.0
      MONTHB = MME*1.0
      DAYB   = IDDE*1.0
      HOURB  = IHHE*1.0 + IMINE/60.0
      JEND   = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
      NDAYS  = INT(JEND-JBEGIN+0.1)+1
      IJS = INT(JBEGIN)
      IJE = INT(JEND)
      WRITE(*,*) 'BEGIN TIME = ',IJS,' END TIME = ',IJE,
     1   ' NDAYS = ',NDAYS
      ALLOCATE(TIME(NCUT,NCYCLE,NDAYS))
      ALLOCATE(FORE(NCUT,NCYCLE,NDAYS,2))

      OPEN(20, FILE = TRIM(FILENAME))
      INQUIRE(FILE=TRIM(FILEO),EXIST=FEXIST)
      IF(.NOT. FEXIST) THEN
        WRITE(*,*) 'filein=',TRIM(FILEO), ' does not exist '
        WRITE(*,*) 'no observation is available'
        WRITE(*,*) 'use tidal prediction as persistence forecast'
        NUMO = INT((JEND-JBEGIN)*24.0/DELT+0.1) + 1
	DO N = 1, NUMO
          TIMEO(N) = JBEGIN+(N-1)*DELT/24.0-JDAY0
          OBS(N,1) = -999.9
          OBS(N,2) = -999.9
        ENDDO  
        GOTO 45
      ENDIF

      OPEN(10, FILE = TRIM(FILEO))
      DO N = 1, NUM
        IF(KINDAT .EQ. 1) THEN
          READ(10,*,END = 40) TIMEO(N),IYR,ICM,ICD,IHR,IMN,SP,DR
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMEO(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
          OBS(N,1) = SP*SIN(DR*3.1415926/180.0)
          OBS(N,2) = SP*COS(DR*3.1415926/180.0)
        ELSE IF(KINDAT .GE. 2)THEN
          READ(10,*,END = 40) TIMEO(n),IYR,ICM,ICD,IHR,IMN,OBS(N,1)
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMEO(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
        ENDIF
      ENDDO
40    CLOSE(10)
      NUMO = N - 1
      WRITE(*,*) 'NUMBER OF OBSERVATION IS ',NUMO
      WRITE(*,*) 'FIST AND LAST TIME ',TIMEO(1), TIMEO(NUMO)

45    CONTINUE
      INQUIRE(FILE = TRIM(FILET), EXIST=FEXIST)
      IF(.NOT. FEXIST) THEN
        WRITE(*,*) 'filein=',trim(filet), ' does not exist '
        WRITE(*,*) 'tidal prediction has to be existed in order to'
        WRITE(*,*) 'make persistence forecast ...'
        STOP
      ENDIF

      OPEN(10, FILE = TRIM(FILET))
      WRITE(*,*) FILET
      DO N = 1, NUM
        IF(KINDAT .EQ. 1) THEN
          READ(10,*,END=50) TIMET(N),IYR,ICM,ICD,IHR,IMN,SP,DR
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMET(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
          TIDES(N,1) = SP*SIN(DR*3.1415926/180.)
          TIDES(N,2) = SP*COS(DR*3.1415926/180.)
        ELSE IF(KINDAT .GE. 2) THEN
          READ(10,*,END=50) TIMET(N),IYR,ICM,ICD,IHR,IMN,TIDES(N,1)
          YEARB  = IYR*1.0
          MONTHB = ICM*1.0
          DAYB   = ICD*1.0
          HOURB  = IHR*1.0 + IMN/60.0
          TIMET(N) = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
        ENDIF
      ENDDO
50    CLOSE(10)

      NUMT = N - 1
      WRITE(*,*)'Number of tide prediction is ',NUMT
      WRITE(*,*)'Fist and last times are ',TIMET(1),TIMET(NUMT)
      OPEN(30,FILE = TRIM(FILEF))
      DO I = 1, NDAYS
        DO J = 1, NCYCLE
          DO N= 1, NCUT
            READ(30,*,END=140) TIMETMP,IYR,ICM,ICD,IHR,IMN
            YEARB  = IYR*1.0
            MONTHB = ICM*1.0
            DAYB   = ICD*1.0
            HOURB  = IHR*1.0 + IMN/60.0
            TIME(N,J,I) = JULIAN(YEARB,MONTHB,DAYB,HOURB)-JDAY0
            IF(N .EQ. 1) THEN
              DO IT = 1, NUMT
                IF(ABS(TIME(N,J,I)-TIMET(IT)) .LT. 0.001 ) THEN
                  IT0 = IT
                  GOTO 433
                ELSE
                  IF(IT .GE. NUMT) THEN
		    WRITE(*,*) 'Time = ',TIME(N,J,I),I,J,N
                    WRITE(*,*) 'No prediction time match with forecast'
                    STOP
                  ENDIF 
                ENDIF
              ENDDO 
433           CONTINUE
              DIFF1 = 0.0
              DIFF2 = 0.0
              DO IO = 1, NUMO
                IF(ABS(TIMET(IT0)-TIMEO(IO)) .LT. 0.001) THEN
                  WRITE(*,"('N,OI,OT=',5I6)") I,J,N,IO,IT
                  WRITE(*,*) TIME(N,J,I),TIMEO(IO),TIMET(IT)
                  IF(OBS(IO,1) .GT. -900) THEN
                    IF(KINDAT .EQ. 1) THEN
                      DIFF1 = OBS(IO,1) - TIDES(IT0,1)                     
                      DIFF2 = OBS(IO,2) - TIDES(IT0,2)                     
                    ELSE IF(KINDAT .GE. 2) THEN
                      DIFF1 = OBS(IO,1) - TIDES(IT0,1)                     
                    ENDIF
                  ELSE
                    DIFF1 = 0.0
                    DIFF2 = 0.0
                  ENDIF 
                  GOTO 110
                ELSE
                  IF(IO .GE. NUMO) THEN
                    WRITE(*,*) 'NO OBSERVATION DATA AT',TIMET(IT0)
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
            ELSE IF(KINDAT .GE. 2) THEN
              FORE(N,J,I,1)=TIDES(IT0,1)+DIFF1
            ENDIF
            IT0 = IT0 + 1
          ENDDO
        ENDDO
      ENDDO 

140   CONTINUE
      IMAXF = I - 1
      DO I = 1, IMAXF
        DO J = 1, NCYCLE
          DO N = 1, NCUT
            IYR = IYRS
            JDAY = TIME(N,J,I)+JDAY0
            CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
            IYEAR = INT(YEARB)
            ICM   = INT(MONTHB)
            ICD   = INT(DAYB)
            IHR   = INT(HOURB)
            IMN   = INT((HOURB-IHR)*60+0.1)
            DAYJ  = TIME(N,J,I) + 1.0
            IF(KINDAT .EQ. 1) THEN
              U = FORE(N,J,I,1)
              V = FORE(N,J,I,2)
              CALL VELDIR (V,U,AANGLE,AVEL)
              IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
              SP = AVEL
              DR = AANGLE
              WRITE(20,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,SP,DR,U,V
            ELSE IF(KINDAT .GE. 2) THEN
              WRITE(20,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,FORE(N,J,I,1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
100   FORMAT(F10.5,I5,4I3,10F10.4)
      CLOSE(20)

      STOP
      END

