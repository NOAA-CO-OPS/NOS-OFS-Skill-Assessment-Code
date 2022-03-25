C lf95 --staticlink READCDEC.f -o READCDEC.x
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
C    19   63680     Turbidity, water, unfiltered, monochrome near infra-red LED light, 780-900 nm, detection angle 90 +/ -2.5 degrees, formazin nephelometric units (FNU)
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

      parameter (NMAX=90000)
      character*120 FIN,FOUT,CORMSLOG,CTMP(200),VARID_USGS(20)
      character*120 BUFFER*400,CVARIN(20),TZONE*3,VARID_IN(20),CTAB*1
      dimension VAL(100),IORDER(20),VAL1(100)
      real*8 jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      LOGICAL FEXIST
      INTEGER DAYS_PER_MONTH(12)

      DATA (DAYS_PER_MONTH(I),I = 1,12) /
     &  31,28,31,30,31,30,31,31,30,31,30,31/ 

      CTAB = CHAR(9)   !!  value of tab separate is 9
      CALL GETARG(1,BUFFER)
      FIN = trim(adjustL(BUFFER))
      print *,'fin= ',TRIM(FIN)
      CALL GETARG(2,BUFFER)
      BUFFER = trim(adjustL(BUFFER))
      LEN = LEN_TRIM(BUFFER)
      N=0
      DO WHILE (LEN > 0)
        N = N+1
        LL = INDEX(TRIM(BUFFER),' ')
	IF(LL .EQ. 0 ) THEN
	  READ(BUFFER(1:LEN),*) CVARIN(N)
	  exit
        ELSE
	  READ(BUFFER(1:LL-1),*) CVARIN(N)
	  BUFFER=BUFFER(LL:LEN)
          BUFFER=trim(adjustL(BUFFER))
          LEN=LEN_TRIM(BUFFER)
	ENDIF
      ENDDO
      NVAR = N	        

      DO N=1,NVAR
        IF(CVARIN(N) .EQ. "DISCHARGE" ) VARID_IN(N) = "00060"
        IF(CVARIN(N) .EQ. "GAGE" ) VARID_IN(N) = "00065"
        IF(CVARIN(N) .EQ. "COND" ) VARID_IN(N) = "00095"
        IF(CVARIN(N) .EQ. "TEMP" ) VARID_IN(N) = "00010"
        IF(CVARIN(N) .EQ. "DO" )   VARID_IN(N) = "00300"
        IF(CVARIN(N) .EQ. "TURBIDITY" ) VARID_IN(N) = "63680"
      ENDDO	
      CALL GETARG(3,BUFFER)
      FOUT=trim(adjustL(BUFFER))

C20    CONTINUE
      INQUIRE(FILE=FIN,EXIST=FEXIST)
      IF(.NOT.FEXIST) STOP
      OPEN(10,FILE=trim(FIN),STATUS='OLD')
      OPEN(20,FILE=trim(FOUT))
      read(10,'(a400)') BUFFER
      LL = INDEX(TRIM(BUFFER),'<tr><td nowrap>')

      DO WHILE (LL .EQ. 0)
        read(10,'(a400)',end=99) BUFFER
        LL = INDEX(TRIM(BUFFER),'<tr><td nowrap>')
      ENDDO

      I = 0
      TZONE = 'PST'

30    read(10,'(a120)',end=99,err=99) BUFFER  
      LL=INDEX(TRIM(BUFFER),'<tr><td nowrap>')
      if (LL.eq.0) goto 30
      LEN=LEN_TRIM(BUFFER)
      LL=INDEX(TRIM(BUFFER),'--')

      if (LL.GT.0) goto 30
      LL=INDEX(TRIM(BUFFER),':')

      if (LL.eq.0) stop
      READ(BUFFER(LL-13:LL+2),101) IMM,IDD,IYR,IHH,IMN
      BUFFER = BUFFER(LL+3:LEN)
      BUFFER = trim(adjustL(BUFFER))
      LL = INDEX(TRIM(BUFFER),"  ")
      READ(BUFFER(LL:LL+9),*) vvv

      YEARB  = IYR
      MONTHB = IMM
      DAYB   = IDD
      HOURB  = IHH + IMN/60.0
         
      IF(TRIM(TZONE) .EQ. 'PST') HOURB = HOURB + 8.0  
      JDAY = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
      IYR = INT(YEARB)
      IMM = INT(MONTHB+0.001)
      IDD = INT(DAYB+0.001)
      IHH = INT(HOURB+0.001)
      IMN = INT((HOURB-IHH)*60+0.1)
      ISEC = 0

C AJ the following  are used to get rid of 60 for seconds and minute, 24 for hours 	
      IF(ISEC .EQ. 60) THEN
        ISEC=0
        IMN=IMN+1
      ENDIF
      IF(IMN .EQ. 60) THEN
        IMN = 0
        IHH = IHH + 1
      ENDIF      
      IF(IHH .EQ. 24) THEN
        IHH = 0
        IDD = IDD + 1
        IF(MOD(IYR,4) .EQ. 0) DAYS_PER_MONTH(2) = 29   !!    Leap Year
        IF(IDD .GT. DAYS_PER_MONTH(IMM)) THEN
          IDD = IDD - DAYS_PER_MONTH(IMM)
          IMM = IMM + 1
	  IF(IMM .GT. 12) THEN
            IMM = IMM - 12
            IYR = IYR + 1
          ENDIF   
        ENDIF
      ENDIF

      IF(CVARIN(1) .EQ. "COND" ) vvv=vvv*0.001
      IF(CVARIN(1) .EQ. "TEMP" ) vvv=vvv*1.0
      WRITE(20,200) IYR,IMM,IDD,IHH,IMN,0,vvv 
      GOTO 30
101   FORMAT(I2,1x,I2,1x,I4,1x,I2,1x,I2)
200   FORMAT(I4,1x,I2,1x,I2,1x,I2,1x,I2,1x,I2,1x,20F15.3)
      
C860   CONTINUE
C99    CONTINUE

99    STOP
      END   


      SUBROUTINE GREGORIAN(jday,yr,month,day,hour)
      real*8 jday,yr,month,day,hour,dayoweek,week

      a = aint(jday+.5)
      b = a+1537
      c = aint((b-122.1)/365.25)
      d = aint(365.25*c)
      e = aint((b-d)/30.6001)
      day = b-d-aint(30.6001*e)+mod(jday+0.5d00,1.0d00)
      month = e-1-12*aint(e/14)
      yr = c-4715-aint((7+month)/10)
      hour = mod(jday+0.5d00,1.0d00)*24.
      day = aint(day)

!      dayoweek = mod(aint(jday+0.5),7)
!      week = aint((jday-2444244.5)/7)
!  2000 leap year crashed.
      if(month.eq.2. .and. day .eq. 31.) then
       yr = yr - 1.
       day = 29.
      endif

      RETURN
      END


      FUNCTION JULIAN(yr,monthb,dayb,hourb)
      real*8 yr,monthb,dayb,hourb
      real*8 Yb,mb,julian

      Yb = yr
      if(yr.lt.100. .and. yr.gt.50.) then
        Yb = yr + 1900.
      elseif(yr.lt.100. .and. yr.le.50.) then
        Yb = yr + 2000.
      endif

      if(monthb.le. 2.) then
        Yb = yb - 1
        mb = monthb + 12.
      else
        Yb = Yb
        mb = monthb
      endif

      JULIAN = aint(365.25*Yb) + aint(30.6001*(mb+1))
     1    + dayb + hourb/24. + 1720981.50

       RETURN
       END

