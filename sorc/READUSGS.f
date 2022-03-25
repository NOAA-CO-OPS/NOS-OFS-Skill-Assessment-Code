C   Documentation for FORTRAN main program READUSGS.f
C
C-------------------------------------------------------------------------------------------------------------------
C
C Fortran Program Name: READUSGS.f
C
C Directory Location:   /COMF/oqcs/sorc
C
C Technical Contact(s):   Name:  Aijun Zhang             Org:  NOS/CO_OPS
C                         Phone: 301-713-2890 X 127      E-Mail: Aijun.zhang@noaa.gov
C
C Abstract:
C               This program is used to replace READUSGS.pl to process data set from USGS web site.
C               The time is converted from local time zone to UTC/GMT
C               Converstion from feet to MKS done here
C               
C Usage:
C      ./READUSGS.x FIN "varlist" FOUT
C      FIN:   Input file name generated from wget
C      varlist: variable list to process. e.g. "DISCHARGE  TEMP COND"
C      FOUT   :  Output file name              
C Compile:
C   lf95 READUSGS.f -o READUSGS.x
C
C variables:  
CCCC  USGS parameter ID  ****************************************
C    01   00065     Gage height, feet
C    09   00055     Stream velocity, feet per second(UVM Path 3)
C    15   00300     Dissolved oxygen, water, unfiltered, milligrams per liter
C    16   00010     Temperature, water, degrees Celsius
C    17   00095     Specific conductance, water, unfiltered, microsiemens per centimeter at 25 degrees Celsius
C    18   00060     Discharge, cubic feet per second
C    19   63680     Turbidity, water, unfiltered, monochrome near infra-red LED light, 780-900 nm, 
C                   detection angle 90 +/ -2.5 degrees, formazin nephelometric units (FNU)
C    20   99905     99905(Chlorophyll u/l)
C    21   99900     99900(Blue Green Alage c/ml)
C***************************************************************
C TEMP   	00010   - TEMPERATURE, WATER (DEG. C)
C COND  	00095   - SPECIFIC CONDUCTANCE (MICROSIEMENS/CM AT 25 DEG. C)
C DISCHARGE   	00060   - DISCHARGE, CUBIC m PER SECOND
C GAGE   	00065   - GAGE HEIGHT, m
C Returns  -9999.00000 if the station doesn't have the data type
C Usage:
C  
C Language:  Fortran
C Revisions:
C         Date          Author         Description
C        11/20/2009      Aijun Zhang    to replace READUSGS.pl   
C  Modified by Lianyuan Zheng on 03/01/2017

      PROGRAM READUSGS
      PARAMETER (NMAX = 90000)
      CHARACTER*120 FIN,FOUT,CORMSLOG,CTMP(200),VARID_USGS(20)
      CHARACTER*120 CVARIN(20),TZONE*3,VARID_IN(20),CTAB*1
      CHARACTER*400 BUFFER,BUFFER1
      DIMENSION VAL(100),IORDER(20),VAL1(100)
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      LOGICAL FEXIST
      INTEGER DAYS_PER_MONTH(12)

      DATA (DAYS_PER_MONTH(i),I = 1,12) /
     &   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/ 

      CTAB = CHAR(9)   !!  value of tab separate is 9
      CALL GETARG(1,BUFFER)
      FIN = TRIM(adjustL(BUFFER))
      CALL GETARG(2,BUFFER)
      BUFFER = TRIM(adjustL(BUFFER))
      LEN0 = LEN_TRIM(BUFFER)
      N = 0
      DO WHILE (LEN0 > 0)
        N = N+1
        LL = INDEX(TRIM(BUFFER),' ')
	IF(LL .EQ. 0 ) THEN
	  READ(BUFFER(1:LEN0),*) CVARIN(N)
	  EXIT
        ELSE
	  READ(BUFFER(1:LL-1),*) CVARIN(N)
	  BUFFER = BUFFER(LL:LEN0)
          BUFFER = TRIM(adjustL(BUFFER))
          LEN0 = LEN_TRIM(BUFFER)
	ENDIF
      ENDDO
      NVAR = N

      DO N = 1, NVAR
        IF(CVARIN(N) .EQ. "DISCHARGE" ) VARID_IN(N) = "00060"
        IF(CVARIN(N) .EQ. "GAGE" ) VARID_IN(N) = "00065"
        IF(CVARIN(N) .EQ. "COND" ) VARID_IN(N) = "00095"
        IF(CVARIN(N) .EQ. "TEMP" ) VARID_IN(N) = "00010"
        IF(CVARIN(N) .EQ. "SALT" ) VARID_IN(N) = "90860"
        IF(CVARIN(N) .EQ. "DO" )   VARID_IN(N) = "00300"
        IF(CVARIN(N) .EQ. "TURBIDITY" ) VARID_IN(N) = "63680"
        IF(CVARIN(N) .EQ. "SEDIMENT" )  VARID_IN(N) = "99409"
      ENDDO	
 
      CALL GETARG(3, BUFFER)
      FOUT = TRIM(adjustL(BUFFER))

      INQUIRE(FILE = FIN, EXIST = FEXIST)
      IF(.NOT. FEXIST) STOP
      OPEN(10,FILE = TRIM(FIN),STATUS = 'OLD')
      OPEN(20,FILE = TRIM(FOUT))

      READ(10,'(a400)') BUFFER
      LL=INDEX(TRIM(BUFFER),' parameter ')
      DO WHILE (LL .EQ. 0)
        READ(10,'(A400)', END=99) BUFFER
        LL=INDEX(TRIM(BUFFER),' parameter ')
      ENDDO

      I = 0
      READ(10,'(A400)', END = 99) BUFFER
      LEN0 = LEN_TRIM(BUFFER)
      DO WHILE(LEN0 .GT. 1)
	I = I+1
	READ(BUFFER,*) CTMP(1), CTMP(2), CTMP(3), CTMP(4)
	VARID_USGS(I) = TRIM(adjustL(CTMP(3)))
        LTMP = INDEX(TRIM(BUFFER),'Temperature')
	IF (LTMP .GT. 0) THEN
          DO N = 1, NVAR
            IF(CVARIN(N) .EQ. "TEMP") VARID_IN(N) = VARID_USGS(I)
          ENDDO	
	ENDIF

        LTMP = INDEX(TRIM(BUFFER),'conductance')
	IF (LTMP .GT. 0) THEN
          DO N = 1, NVAR
            IF(CVARIN(N) .EQ. "COND" ) VARID_IN(N) = VARID_USGS(I)
          ENDDO	
	ENDIF

        LTMP = INDEX(TRIM(BUFFER),'Salinity')
	IF (LTMP .GT. 0) THEN
          DO N = 1, NVAR
            IF(CVARIN(N) .EQ. "SALT" ) VARID_IN(N) = VARID_USGS(I)
          ENDDO	
	ENDIF
	 
        READ(10,'(A400)') BUFFER
        LEN0 = LEN_TRIM(BUFFER)
      ENDDO

      NVAR_USGS = I
      NCOL = NVAR_USGS + 5

      READ(10,'(A400)', END = 99) BUFFER
      LL = INDEX(TRIM(BUFFER),'agency')
      DO WHILE (LL .EQ. 0)
        READ(10,'(A400)', END = 99) BUFFER
        LL = INDEX(TRIM(BUFFER),'agency')
      ENDDO

      READ(BUFFER,*) (CTMP(I),I = 1, NVAR_USGS*2 + 4)
      DO N = 1, NVAR
	IORDER(N) = 99
	DO I = 1, NVAR_USGS
	  IF(TRIM(VARID_IN(N)) .EQ. TRIM(VARID_USGS(I))) THEN
	    IORDER(N) = I
	    EXIT
	  ENDIF 	
	ENDDO     
      ENDDO

30    READ(10,'(A200)',END = 99,ERR = 99) BUFFER
      IF(BUFFER(1:4) .EQ. 'USGS') THEN
        LEN0 = LEN_TRIM(BUFFER)
        N = 0
        DO I = 1, LEN0
          IF(BUFFER(I:I) .EQ. CTAB) THEN
            N = N + 1
            BUFFER1(N:N+4) = ' TAB '  !! convert <tab> into "blank space
            N = N + 4
          ELSE
            N = N + 1
            BUFFER1(N:N) = BUFFER(I:I)
          ENDIF             	      
        ENDDO  

        BUFFER = BUFFER1
        LEN0 = LEN_TRIM(BUFFER)
        LL = INDEX(TRIM(BUFFER),'TAB P ')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+5) = "      "
          LL = INDEX(TRIM(BUFFER),'TAB P ')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'TAB P:e ')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+7) = "        "
          LL = INDEX(TRIM(BUFFER),'TAB P:e ')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'TAB P:< ')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+7) = "        "
          LL = INDEX(TRIM(BUFFER),'TAB P:< ')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'***')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+2) = "999"
          LL = INDEX(TRIM(BUFFER),'***')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'TAB A ')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+5) = "      "
          LL = INDEX(TRIM(BUFFER),'TAB A ')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'TAB  TAB ')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+7) = "-999.99  "
          LL = INDEX(TRIM(BUFFER),'TAB  TAB ')
        ENDDO
        LL = INDEX(TRIM(BUFFER),'TAB')
        DO WHILE (LL .NE. 0)
          BUFFER(LL:LL+2) = "   "
          LL = INDEX(TRIM(BUFFER),'TAB')
        ENDDO

        READ(BUFFER,*) (CTMP(I),I = 1,NCOL)
	LTMP = LEN_TRIM(CTMP(NCOL))
        IF(LTMP .LE. 0) THEN
	  PRINT*,'data less colomn: ',TRIM(BUFFER)
          STOP
	ENDIF 

        LL = INDEX(TRIM(BUFFER),':')
	READ(CTMP(3),100) IYR, IMM, IDD
	READ(CTMP(4),'(I2,1x,I2)') IHH, IMN
	READ(CTMP(5),'(a3)') TZONE
        YEARB  = IYR
        MONTHB = IMM
        DAYB   = IDD
        HOURB  = IHH + IMN/60.0

Cc convert local time to UTC/GMT time	  
        IF(TRIM(TZONE) .EQ. 'EST') HOURB = HOURB + 5.0
        IF(TRIM(TZONE) .EQ. 'EDT') HOURB = HOURB + 4.0
        IF(TRIM(TZONE) .EQ. 'CST') HOURB = HOURB + 6.0  !Central Standard Time
        IF(TRIM(TZONE) .EQ. 'CDT') HOURB = HOURB + 5.0  !Central Daylight Time
        IF(TRIM(TZONE) .EQ. 'MST') HOURB = HOURB + 7.0  
        IF(TRIM(TZONE) .EQ. 'MDT') HOURB = HOURB + 6.0  
        IF(TRIM(TZONE) .EQ. 'PST') HOURB = HOURB + 8.0  
        IF(TRIM(TZONE) .EQ. 'PDT') HOURB = HOURB + 7.0  
        IF(TRIM(TZONE) .EQ. 'HST') HOURB = HOURB + 10.0 
        IF(TRIM(TZONE) .EQ. 'HADT')HOURB = HOURB + 9.0 !! not work since TZONE is defined as 3 characters
        JDAY = JULIAN(YEARB,MONTHB,DAYB,HOURB)
        CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
        IYR  = INT(YEARB)
        IMM  = INT(MONTHB+0.001)
        IDD  = INT(DAYB+0.001)
        IHH  = INT(HOURB+0.001)
        IMN  = INT((HOURB-IHH)*60+0.1)
        ISEC = 0

C AJ the following are used to get rid of 60 for seconds and minute, 24 for hours 	
	IF(ISEC .EQ. 60) THEN
	  ISEC = 0
	  IMN  = IMN + 1
 	ENDIF
	IF(IMN .EQ. 60) THEN
	  IMN = 0
	  IHH = IHH + 1
	ENDIF      
        IF(IHH .EQ. 24) THEN
          IHH = 0
	  IDD = IDD + 1
	  IF(MOD(IYR,4) .EQ. 0) DAYS_PER_MONTH(2) = 29   !!    Leap Year
	  IF(IDD .GT. DAYS_PER_MONTH(IMM) ) THEN
	    IDD = IDD - DAYS_PER_MONTH(IMM)
	    IMM = IMM + 1
	    IF(IMM .GT. 12) THEN
	      IMM = IMM - 12
	      IYR = IYR + 1
	    ENDIF   
          ENDIF
        ENDIF

	DO I = 1, NVAR_USGS
	  READ(CTMP(I+5),*) VAL(I)
	ENDDO

	DO L = 1,NVAR
	  VAL1(L)=-999.9
	  IF(IORDER(L) .LE.NVAR_USGS) THEN
	    VAL1(L) = VAL(IORDER(L))
            IF(CVARIN(L) .EQ. "DISCHARGE") VAL1(L) = VAL1(L)*0.3048**3
            IF(CVARIN(L) .EQ. "GAGE") VAL1(L) = VAL1(L)*0.3048
            IF(CVARIN(L) .EQ. "COND") VAL1(L) = VAL1(L)*0.001
            IF(CVARIN(L) .EQ. "TEMP") VAL1(L) = VAL1(L)*1.0
            IF(CVARIN(L) .EQ. "SALT") VAL1(L) = VAL1(L)*1.0
	  ELSE
	    VAL1(L) = -999.9
	  ENDIF   
	ENDDO
	IF(VAL1(1) .LE. -0.999.AND.CVARIN(1) .EQ. 'COND') GOTO 30
	IF(VAL1(1) .LE. -999.0) GOTO 30
	WRITE(20,200) IYR,IMM,IDD,IHH,IMN,0,(VAL1(L),L = 1,NVAR) 
      ENDIF	
      GOTO 30
      
100   FORMAT(I4,1x,I2,1x,I2)
200   FORMAT(I4.4,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,20F15.3)

99    CLOSE(10)
      CLOSE(20)

      STOP
      END   


      SUBROUTINE GREGORIAN(jday,yr,month,day,hour)
      REAL*8 JDAY,YR,MONTH,DAY,HOUR,DAYOWEEK,WEEK

      A = AINT(JDAY+0.5)
      B = A+1537
      C = AINT((B-122.1)/365.25)
      D = AINT(365.25*C)
      E = AINT((B-D)/30.6001)
      DAY = B-D-AINT(30.6001*E)+MOD(JDAY+0.5D00,1.0D00)
      MONTH = E-1-12*AINT(E/14)
      YR = C-4715-AINT((7+MONTH)/10)
      HOUR = MOD(JDAY+0.5D00,1.0D00)*24.0
      DAY = AINT(DAY)

      IF(MONTH .EQ. 2.0 .AND. DAY .EQ. 31.0) THEN
        YR = YR-1.0
        DAY = 29.0
      ENDIF

      RETURN
      END


      FUNCTION JULIAN(YR,MONTHB,DAYB,HOURB)
      REAL*8 YR,MONTHB,DAYB,HOURB
      REAL*8 YB,MB,JULIAN

      YB = YR
      IF(YR .LT. 100.0 .AND. YR .GT. 50.0) THEN
        YB = YR + 1900.0
      ELSEIF(YR .LT. 100.0 .AND. YR .LE. 50.0) THEN
        YB = YR + 2000.0
      ENDIF

      IF(MONTHB .LE. 2.0) THEN
        YB = YB-1.0
        MB = MONTHB + 12.0
      ELSE
        YB = YB
        MB = MONTHB
      ENDIF

      JULIAN = AINT(365.25*YB) + AINT(30.6001*(MB+1))
     1    + DAYB + HOURB/24.0 + 1720981.50

       RETURN
       END

