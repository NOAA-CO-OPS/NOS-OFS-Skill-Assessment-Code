C   Documentation for FORTRAN program fillnan.f
C
C--------------------------------------------------------------------------------------
C
C Fortran Program Name: fillnan.f
C
C Directory Location:   /COMF/oqcs/sorc  
C
C Technical Contact(s):   Name:  Aijun Zhang             Org:  NOS/COOPS
C                         Phone: 301-713-2890x127      E-Mail: Aijun.Zhang@noaa.gov
C Abstract:	
C 	  Reads in a file with lines of dates,data and fills missing
C 	  time delta's with NAN
C
C	  This FORTRAN program is used to replace fillnan.c, which doesn't work well.
C			   
C Usage:  
C	  cat data.dat | fillnan.x "$t1" "$t2" dt "Nan Nan"  > filled_data.dat
C	  cat data.dat | fillnan.x "2005 03 10 12 00" "2005 03 12 12 00" 0.10 "Nan Nan"  > filled_data.dat
C	Ex.,
C         cat data.dat | ./fillnan.x "2005 03 10 13 30" "2005 03 10 15 00" 0.10 "Nan Nan"  > filled_data.dat
C         cat data.dat
C                  2005 03 10 14 06
C                  2005 03 10 14 12
C                  2005 03 10 14 18
C                  2005 03 10 14 24
C                  2005 03 10 14 30
C         cat filled_data.dat
C                  2005 03 10 13 30 Nan Nan 
C                  2005 03 10 13 36 Nan Nan 
C                  2005 03 10 13 42 Nan Nan 
C                  2005 03 10 13 48 Nan Nan 
C                  2005 03 10 13 54 Nan Nan 
C                  2005 03 10 14 00 Nan Nan 
C                  2005 03 10 14 06  
C                  2005 03 10 14 12  
C                  2005 03 10 14 18  
C                  2005 03 10 14 24  
C                  2005 03 10 14 30  
C                  2005 03 10 14 36 Nan Nan 
C                  2005 03 10 14 42 Nan Nan 
C                  2005 03 10 14 48 Nan Nan 
C                  2005 03 10 14 54 Nan Nan 
C
C Input Parameters:  
C
C Language:	FORTRAN	  
C
C Compiling/Linking Syntax: 
C   lf95 --staticlink fillnan.f -L/opt/comf/oqcs/binlinux -loqcs -o fillnan.x
C
C Target Computer:  Runs on COMF computers, such as ofsdev  
C
C Estimated Execution Time: 
C
C Suboutines/Functions Called:
C                  Name         Directory Location        Description
C
C Input Files:
C  Unit No.        Name         Directory Location        Description
C
C Output Files:
C  Unit No.        Name         Directory Location        Description
C
C Libraries Used:     
C
C Error Conditions:
C
C Creation Date:   2013-09-27
C Modified by Lianyuan Zheng on 03/01/2017
C
C Remarks: 
C
C Subroutines/Functions Called:
C      Name       Directory Location          Description
C     julian      /COMF/oqcs/sorc/library     Convert Gregorian dates to Julian days.
C     gregorian   /COMF/oqcs/sorc/library     Convert Julian days to Gregorian.
C---------------------------------------------------------------------------------------

      PROGRAM FILLNAN
      IMPLICIT NONE
      
      INTEGER NARRAY
      PARAMETER(NARRAY=900000)
      CHARACTER*120 FILEOBS,FILEOUT,ARGV,BUFFER,SFILLVALUES,SOBS(NARRAY)
      INTEGER I,N,N0,NUMTIMES,IYEAR,NOBS,IARGC
      INTEGER IMONTH,IDAY,IHOUR,IMIN,K,FH
      INTEGER NCOL,LEN0,LL
      REAL*8 IST_YR, IST_MON, IST_DAY, IST_HR, IST_MIN
      REAL*8 LST_YR, LST_MON, LST_DAY, LST_HR, LST_MIN
      REAL*8 JDAYFIRST,JDAYLAST
      REAL*8 DT,DTJ,WLSTART,FHOUR
      REAL*8 OBS(NARRAY),JDAYOBS(NARRAY),OBS1(NARRAY),OBS2(NARRAY)
      REAL*8 OBSOUT(NARRAY),JDAYOUT(NARRAY)
      REAL*8 JULIAN, DY ,YEAR,DAY,MONTH,HOUR,MIN,SEC,YDAY   
      REAL*8 JDAYFIRST_OBS,FILLVALUES
      REAL*8 JDAYLAST_OBS
      REAL TMP1,TMP2,TMP3

      N = IARGC()      
      CALL GETARG(1,ARGV)
      ARGV = TRIM(ADJUSTL(ARGV))
      READ(ARGV,*) IST_YR, IST_MON, IST_DAY, IST_HR, IST_MIN
      JDAYFIRST = JULIAN(IST_YR,IST_MON,IST_DAY,IST_HR+IST_MIN/60.0)

      CALL GETARG(2,ARGV)
      ARGV = TRIM(ADJUSTL(ARGV))
      READ(ARGV,*) LST_YR, LST_MON, LST_DAY, LST_HR, LST_MIN
      JDAYLAST = JULIAN(LST_YR,LST_MON,LST_DAY,LST_HR+LST_MIN/60.0)

      CALL GETARG(3,ARGV)
      ARGV = TRIM(ADJUSTL(ARGV))
      READ(ARGV,*) DT

      CALL GETARG(4,ARGV)
      SFILLVALUES = TRIM(ADJUSTL(ARGV))
      READ(SFILLVALUES,*) FILLVALUES

      NCOL = 1
      LL   = 1
      BUFFER = TRIM(ADJUSTL(ARGV))
      DO WHILE (LL .GT. 0)
        LEN0 = LEN_TRIM(BUFFER)
        LL  = INDEX(TRIM(BUFFER),' ')      
        IF(LL .GT. 1 .AND. LL .LE. LEN0) THEN
          BUFFER = BUFFER(LL:LEN0)
          BUFFER = TRIM(ADJUSTL(BUFFER))
          READ(BUFFER,*) TMP1
          IF(TMP1 .GT. -99.9) GOTO 12
          NCOL = NCOL+1
        ENDIF
      ENDDO
12    CONTINUE

C  Build the output time variable
      DTJ = DT/1440.0     ! input dt is in minutes in SA tool
      NUMTIMES = (JDAYLAST-JDAYFIRST)/DTJ +1
      DO I = 1,NUMTIMES
        JDAYOUT(I) = JDAYFIRST + DTJ*(I-1)
      ENDDO

      I = 1
      DO WHILE (I .LT. NARRAY) 
        IF(NCOL .EQ. 1) THEN
          READ(*,*,END=1000) YEAR,MONTH,DAY,HOUR,MIN,OBS(I)
        ELSEIF(NCOL .EQ. 2) THEN
          READ(*,*,END=1000) YEAR,MONTH,DAY,HOUR,MIN,OBS(I),OBS1(I)
        ELSEIF(NCOL .EQ. 3) THEN
          READ(*,*,END=1000) YEAR,MONTH,DAY,HOUR,MIN,OBS(I),OBS1(I),
     1      OBS2(I)
        ENDIF

        JDAYOBS(I) = JULIAN(YEAR,MONTH,DAY,HOUR+MIN/60.0D00)
        I = I + 1
      ENDDO
1000  NOBS = I - 1
       
C  Output the lines of data
      DO I = 1, NUMTIMES
        BUFFER = SFILLVALUES
        TMP1 = FILLVALUES
        TMP2 = FILLVALUES
        TMP3 = FILLVALUES
        N0 = 1
        DO N = N0, NOBS
	  IF(ABS(JDAYOUT(I)-JDAYOBS(N)) .LE. 0.001) THEN
            N0 = N
            IF(NCOL .EQ. 1) THEN
              TMP1 = OBS(N)
            ELSEIF(NCOL .EQ. 2) THEN
              TMP1 = OBS(N)
              TMP2 = OBS1(N)
            ELSEIF(NCOL .EQ. 3) THEN
              TMP1 = OBS(N)
              TMP2 = OBS1(N)
              TMP3 = OBS2(N)
            ENDIF
	    GOTO 80
	  ENDIF  
        ENDDO

80      CONTINUE	
	CALL GREGORIAN(JDAYOUT(I),YEAR,MONTH,DAY,HOUR)
        IYEAR  = YEAR + 0.1
        IMONTH = MONTH + 0.1
        IDAY   = DAY + 0.1
        HOUR   = HOUR + 0.01/3600  
        IHOUR  = HOUR
        IMIN   = (HOUR-IHOUR)*60.0 + 1.0/60.0
        IF(NCOL .EQ. 1) THEN
          WRITE(*,77) IYEAR,IMONTH,IDAY,IHOUR,IMIN,TMP1
        ELSEIF(NCOL .EQ. 2) THEN
          WRITE(*,77) IYEAR,IMONTH,IDAY,IHOUR,IMIN,TMP1,TMP2
        ELSEIF(NCOL .EQ. 3) THEN
          WRITE(*,77) IYEAR,IMONTH,IDAY,IHOUR,IMIN,TMP1,TMP2,TMP3
        ENDIF
      ENDDO
77    FORMAT(I4,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,2X,10F12.4)

      STOP      
      END
 
