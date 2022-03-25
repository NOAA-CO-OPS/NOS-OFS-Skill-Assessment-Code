C      compile : lf95 lsqha.f -o lsqha.x
C      run by command: lsqha.x $KINDAT  $NCON $DELT $LONGITUDE $FILEIN
CC     KINDAT:  =1 for vector data (current speed and direction)
C               =2 dor scalar (water level)
C      NCON     = number of constituents to be analyzed by H.A., maximum=37
C      DELT:    equally-spaced data time interval in hours, =0.1 for 6 minutes time interval
C      Longitude:  longitude of the station in decimal degrees (positive for west longitude)
C      FLOODDIR    flood direction in degrees
C      FILEIN:    input data file name                
C                
      PARAMETER(NT = 200000)
      CHARACTER*10 LABEL,IDENS,TKIND,SKIND,LABLE
      character*80 filestation,outputf0,FNAME,AQ,BUFFER*100 
      DIMENSION TIME(NT),VIN(NT),DIN(NT),SIYR(NT),SMON(NT)
     1          ,SDD(NT),SHH(NT),SMIN(NT),ifrst(NT),ilast(NT)
     2          ,XTMP(NT),YTMP(NT),U(NT),V(NT)
      real*8 jday,jdayb,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      REAL*8 SPEED,SHOUR,DELTT,DELTAX
      COMMON/PARAM1/NCONST,NYYYY,NDELTT,CUTOFF,TZERO,MINTRM,MAXT
     1RM,NCOL,NROW,TEBAR
      COMMON/PARAM2/NOBS,MORE,JYY,JMM,JDD,JHH,JMIN,THBAR,XLONG,TIMEZ
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(NT,2),RATIO
      COMMON/PARAM4/IYY(20),IMM(20),IDD(20),IHH(20),IMIN(20)
      COMMON/PARAM5/YEPOCH(180),YNODE(180),SPEED(180),P(180)
      COMMON/PARAM6/IDENS(16),LABEL(180),TKIND,SKIND
      COMMON/BGROUP/DELTT,SHOUR(20),DELTAX
      COMMON/MEANL1/MNYEAR,JYEAR,JDAZ,DAY,LDAYY,LLDAYY,MONTH
      COMMON/MEANL2/TIMEX,JDAAYI,JDAYI,JDAYY,LDAA,LXYERE
      COMMON/MEANL3/IALL,LPYER,LYDAY,LJDAA,LNDAY,LGYER,LZDAY
      COMMON/MEANL4/NJOBX,MMONTH,NXDAY,JDAZZ,JDAAZ,FST,INDATA
      COMMON/PARAMC/RAMP(180),RKAPP(180),LABLE(180),RKAP(180)
      WRITE(*,*)'RUN LSQHA.f'
      TM=0.0
      CALL GETARG(1,BUFFER)
      READ(BUFFER,*)KINDAT
      CALL GETARG(2,BUFFER)
      READ(BUFFER,*)NCON 
      CALL GETARG(3,BUFFER)
      READ(BUFFER,*)DELT 
      DELT=DELT/60.0  !! 11/30/2006  convert from minutes into hours
      CALL GETARG(4,BUFFER)
      READ(BUFFER,*)GONL  
!      CALL GETARG(5,BUFFER)
!      READ(BUFFER,*)FLOODDIR  
      CALL GETARG(5,FNAME)
      if (GONL .lt. 0.0)GONL=-GONL
!      call ncrght(FNAME,ncut)
      AQ='free format ascii file'
11    continue
100     format(a7,1x,a4,2x,a40,2x,f8.3)
110     format(I2,1x,F5.2,1x,a1,1x,I3,1x,F5.2,1x,a1) 
      CLOSE(2)
      OPEN(2,file=trim(FNAME))
      N=0
      IF(KINDAT .EQ. 1)THEN
         I=1
5        READ(2,*,END=10)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I),DIN(I) 
         IF(VIN(I) .GT. 10.0)GOTO 5      
         IF (I .EQ. 1)THEN
          yearb=SIYR(1)
          monthb=1.
          dayb=1.
          hourb=0.
          jbase_date=JULIAN(yearb,monthb,dayb,hourb)
         ENDIF
         yearb=SIYR(I)
         monthb=SMON(I)
         dayb=SDD(I)
         hourb=SHH(I)+SMIN(I)/60.
         jdayb=JULIAN(yearb,monthb,dayb,hourb)
         TIME(I)=(jdayb-jbase_date)*24.
         I=I+1
         GOTO 5
10      NMAX=I-1
!! fill up the small gaps in the original time series, i.e.  criteria1 < gap < criteria2 
         tstart=TIME(1)/24.
	 tfinish=TIME(NMAX)/24.
	 IMAXB=INT((tfinish-tstart)*24/delt)
	 method=1     ! using Singular Value Decomposition (SVD)
	 criteria1=2.0  ! in hours
	 criteria2=6.0  ! in hours
	 print *,'nmax= ',nmax,'Imaxb= ',Imaxb
         DO I=1,NMAX
            sp=VIN(I)
            dr=DIN(I)
            u(I)=sp*sin(dr*3.1415926/180.)
            v(I)=sp*cos(dr*3.1415926/180.)
         ENDDO 
         CALL equal_interval(tstart,tfinish,delt,delt,method,
     1        criteria1,criteria2,TIME,u,XTMP,YTMP,NMAX,IMAXB)
         DO I=1,IMAXB
            U(I)=YTMP(I)
	 ENDDO   
         CALL equal_interval(tstart,tfinish,delt,delt,method,
     1        criteria1,criteria2,TIME,v,XTMP,YTMP,NMAX,IMAXB)
         DO I=1,IMAXB
            V(I)=YTMP(I)
	    TIME(I)=XTMP(I)
	 ENDDO   

        DO I=1,IMAXB
           CALL VELDIR(V(I),U(I),DIN(I),VIN(I))
        ENDDO
        II=0
	NMAX=IMAXB
        DO I=1,NMAX
          IF ((U(i) .GT. -900.) .and. (V(i) .GT. -900.) )then
            II=II+1
            TIME(II)=TIME(I)
            jday=time(i)/24.+jbase_date
            call GREGORIAN(jday,yearb,monthb,dayb,hourb)
            SIYR(II)=INT(yearb)
            SMON(II)=int(monthb+0.001)
             SDD(II)=INT(dayb+0.001)
             SHH(II)=INT(hourb+0.001)
            SMIN(II)=INT((hourb-INT(hourb+0.001))*60+0.1)
             VIN(II)=VIN(I)
             DIN(II)=DIN(I)
	     U(II)=U(I)
	     V(II)=V(I)
          ENDIF
        ENDDO
        NMAX=II
!        DO I=1,NMAX
!	  tday=TIME(I)/24.
!	  IYY1=INT(SIYR(I))
!	  IMM1=INT(SMON(I))
!	  IDD1=INT(SDD(I))
!	  IHH1=INT(SHH(I))
!	  IMN1=INT(SMIN(I))
!	  WRITE(46,344)tday,IYY1,IMM1,IDD1,IHH1,IMN1
!     1                   ,VIN(I),DIN(I),U(I),V(I)
!        ENDDO
344   FORMAT(F10.5,I5,4I3,4F10.4)	
      ELSE IF(KINDAT .EQ. 2)THEN
         I=1
15       READ(2,*,END=20)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I)
         IF(VIN(I) .GT. 10.0)GOTO 15      
         IF (I .EQ. 1)THEN
          yearb=SIYR(1)
          monthb=1.
          dayb=1.
          hourb=0.
          jbase_date=JULIAN(yearb,monthb,dayb,hourb)
         ENDIF
         yearb=SIYR(I)
         monthb=SMON(I)
         dayb=SDD(I)
         hourb=SHH(I)+SMIN(I)/60.
         jdayb=JULIAN(yearb,monthb,dayb,hourb)
         TIME(I)=(jdayb-jbase_date)*24.
         I=I+1
         GOTO 15
20      NMAX=I-1
!! fill up the small gaps in the original time series, i.e.  criteria1 < gap < criteria2 
         tstart=TIME(1)/24.
	 tfinish=TIME(NMAX)/24.
	 IMAXB=INT((tfinish-tstart)*24/delt)
	 method=1     ! using Singular Value Decomposition (SVD)
	 criteria1=2.0  ! in hours
	 criteria2=6.0  ! in hours
         CALL equal_interval(tstart,tfinish,delt,delt,method,
     1        criteria1,criteria2,TIME,VIN,XTMP,YTMP,NMAX,IMAXB)
         DO I=1,IMAXB
            VIN(I)=YTMP(I)
	    TIME(I)=XTMP(I)
	 ENDDO   
        II=0
	NMAX=IMAXB
        DO I=1,NMAX
          IF (VIN(i) .GT. -900.  )then
            II=II+1
            TIME(II)=TIME(I)
            jday=time(i)/24.+jbase_date
            call GREGORIAN(jday,yearb,monthb,dayb,hourb)
            SIYR(II)=INT(yearb)
            SMON(II)=int(monthb+0.001)
             SDD(II)=INT(dayb+0.001)
             SHH(II)=INT(hourb+0.001)
            SMIN(II)=INT((hourb-INT(hourb+0.001))*60+0.1)
             VIN(II)=VIN(I)
          ENDIF
        ENDDO
        NMAX=II
      ENDIF
      gap=1.5*DELT
      call continuous(time,nmax,gap,Nsegments,ifrst,ilast)
      NSMAX=0
      do NG=1,Nsegments
         Istart=ifrst(NG)
         IEND=ilast(NG)
         NDIF=IEND-Istart
         IF (NDIF .GT. NSMAX)THEN
           NSMAX=NDIF
           INDEX=NG
         ENDIF
	 Day_duration=time(IEND)/24.-time(ISTART)/24.
	 print *,'NG= ',NG,ISTART,IEND,NDIF,time(ISTART)/24.
     1	 ,time(IEND)/24.,Day_duration
      ENDDO
      NMAX=ilast(INDEX)-ifrst(INDEX)+1
      IBEGIN=ifrst(INDEX)
      IYRS=SIYR(IBEGIN)
      MONS=SMON(IBEGIN)
      IDDS=SDD(IBEGIN)
      IHHS=SHH(IBEGIN)
      MINS=SMIN(IBEGIN)
      IF(KINDAT .EQ. 1)THEN
        DO I=1,NMAX
          RDATA(I,1)=VIN(IBEGIN+I-1)
          RDATA(I,2)=DIN(IBEGIN+I-1)
          XTMP(I)=VIN(I)*sin(DIN(I)*3.14159265/180.)
          YTMP(I)=VIN(I)*cos(DIN(I)*3.14159265/180.)
        ENDDO 
C    calculate principle current direction
        CALL prcmp(NMAX,XTMP,YTMP,pcd,rxy,RATIO)
!        pcd=FLOODDIR
      ELSE IF(KINDAT .GE. 2)THEN
        DO I=1,NMAX
         RDATA(I,1)=VIN(IBEGIN+I-1)
         RDATA(I,2)=0.0
        ENDDO
      ENDIF
C    prepare control file for LSQHA.f
      CUTOFF=0.000001
      NBLK=1
      NJOBX=INT(NMAX*delt/24.)
      AZI=pcd
      N=NMAX
      NSPH=int(1./delt+0.0001)
      CVAR=0.0
      umean=0.0
      vmean=0.0
      INDATA=1
      VFAC=1.0
      ITYPE=2
C*****************************
      NJ=1
      IAND=0
      IREF=0
      IIT=0
      ISKIP9=0
      IEL=KINDAT-1
C****************************
      WRITE(*,*)'first data point is',IYRS,MONS,IDDS,IHHS,MINS
      IYRE=SIYR(IBEGIN+NMAX-1)
      MONE=SMON(IBEGIN+NMAX-1)
      IDDE=SDD(IBEGIN+NMAX-1)
      IHHE=SHH(IBEGIN+NMAX-1)
      MINE=SMIN(IBEGIN+NMAX-1)
      WRITE(*,*)'last  data point is',IYRE,MONE,IDDE,IHHE,MINE
      write(*,*)'total data length is',NJOBX,' days'
        OPEN(20,file='ha_analysis.ctl')
C      IF (NJOBX .GT. 60)THEN
        WRITE(20,'(a80)')FNAME
        WRITE(20,500)CUTOFF,NBLK,NJOBX
        WRITE(20,510)AZI,N,NSPH,CVAR,UMEAN,VMEAN
        WRITE(20,520)IYRS,MONS,IDDS,IHHS,MINS,TM,GONL
        WRITE(20,530)INDATA,KINDAT,VFAC,NCON,ITYPE
        WRITE(20,'(a40)')AQ
        CLOSE(20)
        OPEN(25,file='ha_analysis.ctl')
        CALL LSQHA3
        CLOSE(25)
!      ELSE IF (NJOBX .LE. 60)THEN
!        WRITE(20,'(a80)')FNAME
!        WRITE(20,600)NJ,IAND,NSPH,IREF,IIT
!        WRITE(20,610)AZI,TM,GONL,CVAR,IEL
!        WRITE(20,'(a40)')AQ
!        WRITE(20,620)IYRS,MONS,IDDS,IHHS,MINS,ISKIP9,N
!        OPEN(25,file='ha_analysis.ctl')
!        CALL HARM29D
!        CLOSE(25)
!      ENDIF
500   FORMAT(f10.6,2I6)
510   FORMAT(f8.3,I8,I4,1x,3f8.3)
520   FORMAT(I5,4I4,f5.1,f9.3)
530   FORMAT(2I4,f8.3,2I4)
600   FORMAT(5I6)
610   FORMAT(f8.3,1x,3f8.3,I5)
620   FORMAT(I5,4I4,I4,I8)
      STOP
      END
C*****************************************************************************
      SUBROUTINE CMSV(A,NOBS,AMEAN,UMEAN,LMEAN,ASD,USD,LSD)
      REAL A(NOBS)
      REAL AMEAN,UMEAN,LMEAN,ASD,USD,LSD
      EXTERNAL TDIS
C
C  COMPUTE THE MEAN VALUE
C
      RNOBS=FLOAT(NOBS)
      SUM = 0.0
      DO 40 I = 1, NOBS
      SUM = SUM + A(I)
   40 CONTINUE
      AMEAN = SUM/RNOBS
C
C  DE-MEAN TIME SERIES
C
      DO 50 I = 1, NOBS
      A(I) = A(I) - AMEAN
   50 CONTINUE
C
C  CALCULATION OF THE STANDARD DEVIATION
C
      SUM = 0.0
      DO 60 I = 1, NOBS
      SUM = SUM + A(I)*A(I)
   60 CONTINUE
      ASD = SQRT(SUM/FLOAT(NOBS - 1))
      AVAR=ASD*ASD
      N=NOBS-1
      CALL CONINT(N,CUPPER,CLOWER)
      USD=SQRT(CUPPER*AVAR)
      LSD=SQRT(CLOWER*AVAR)
      RN=FLOAT(N)

      UMEAN=AMEAN+ASD*TDIS(RN,1)/SQRT(RNOBS)
      LMEAN=AMEAN-ASD*TDIS(RN,1)/SQRT(RNOBS)
C
C RESTORE INPUT TIME SERIES
C
      DO 70 I = 1, NOBS
      A(I) = A(I) + AMEAN
   70 CONTINUE
      RETURN
      END
        FUNCTION TDIS(EF,IAL)
C
C --- PERCENTILE FOR T-DISTRIBUUTION WITH EF DEGREES OF FREEDOM
C ---    IAL (1,2) FOR 95% OR 80% PROBABILITY LEVEL. FROM
C ---    ABRAMOWITZ P949
C
        X=1.96
        IF(IAL.EQ.2)X=1.282
        X3=X**3
        X5=X**5
        X7=X**7
        G1=(X3+X)/4
        G2=(5*X5+16*X3+3*X)/96
        G3=(3*X7+19*X5+17*X3-15*X)/384
        TDIS=X+(G1+(G2+G3/EF)/EF)/EF
        RETURN
        END
      SUBROUTINE CONINT(NDOF,CUPPER,CLOWER)
C
C  THIS SUBROUTINE CALCULATES A 95% CONFIDENCE INTERVAL FOR SPECTRAL 
C  ESTIMATES GIVEN THE NUMBER OF DEGREES OF FREEDOM. THE MULTIPLICATIVE 
C  CONFIDENCE INTERVAL FACTORS ARE CALCULATED FROM VALUES OF THE CHI 
C  SQUARE DISTRIBUTION.
C
C  DEGREES OF FREEDOM MUST BE LESS THAN OR EQUAL TO 240.
C
C      VARIABLE LIST
C
C  NAVG          THE NUMBER OF AVERAGING OVER SPECTRAL WINDOW
C  NDOF          THE NUMBER OF DEGREES OF FREEDOM.
C  CUPPER        THE MULTIPLICATIVE FACTOR USED IN DETERMINING THE UPPER
C                VALUE OF THE CONFIDENCE INTERVAL.
C  CLOWER        THE MULTIPLICATIVE FACTOR USED IN DETERMINING THE LOWER
C                VALUE OF THE CONFIDENCE INTERVAL.
C  RNDOF         SAME AS NDOF BUT REAL.
C  CP025         ARRAY CONTAINING VALUES OF THE CHI SQUARE DISTRIBUTION
C                USED IN CALCULATING THE LOWER VALUE OF THE CONFIDENCE
C                INTERVAL.
C  CP975         ARRAY CONTAINING VALUES OF THE CHI SQUARE DISTRIBUTION
C                USED IN CALCULATING THE UPPER VALUE OF THE CONFIDENCE
C                INTERVAL.
C
      REAL CP025(33),CP975(33)
      DATA CP975/0.00098,0.0506,0.216,0.484,0.831,1.24,1.69,2.18,2.70,
     *           3.25,3.82,4.40,5.01,5.63,6.26,6.91,7.56,8.23,8.91,9.59,
     *          10.28,10.98,11.69,12.40,13.12,13.84,14.57,15.31,16.05,
     *          16.79,24.43,40.48,91.58/
      DATA CP025/5.02, 7.38, 9.35,11.14,12.83, 14.45,16.01,17.53,19.02,
     *          20.48,21.92,23.34,24.47,26.12, 27.49,28.85,30.19,31.53,
     *          32.85,34.17,35.48,36.78,38.08, 39.36,40.65,41.92,43.19,
     *          44.46,45.72,46.98,59.34,83.30,152.21/ 
C
C --- COMPUTE THE NUMBER OF DENSITY SPECTRA AVERAGED OVER THE SPECTRAL 
C     WIDTH
C
      NAVG=NDOF
C
C  IF LESS THAN 61 DEGREES OF FREEDOM, SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS USING VALUES OF THE CHI SQUARE DIS-
C  TRIBUTION.
C
      RNAVG=FLOAT(NAVG)
      IF(NAVG.GT.30)GO TO 100
      CUPPER=RNAVG/CP975(NAVG)
      CLOWER=RNAVG/CP025(NAVG)
      GO TO 300
C
C  IF LESS THAN 80 DEGREES OF FREEDOM(BUT MORE THAN 60),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
  100 IF(NAVG.GT.40) GO TO 200
      J=30
      CUPPER=RNAVG/(CP975(J)+AMOD(RNAVG,10.)/10.*(CP975(J+1)-CP975(J)))
      CLOWER=RNAVG/(CP025(J)+AMOD(RNAVG,10.)/10.*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF LESS THAN 120 DEGREES OF FREEDOM(BUT MORE THAN 80),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
  200 IF(NAVG.GT.60) GO TO 210
      J=31
      CUPPER=RNAVG/(CP975(J)+AMOD(RNAVG,20.)/20.*(CP975(J+1)-CP975(J)))
      CLOWER=RNAVG/(CP025(J)+AMOD(RNAVG,20.)/20.*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF LESS THAN 240 DEGREES OF FREEDOM(BUT MORE THAN 120),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
  210 IF(NAVG.GT.120) GO TO 220
      J=32
      CUPPER=RNAVG/(CP975(J)+AMOD(RNAVG,60.)/60.*(CP975(J+1)-CP975(J)))
      CLOWER=RNAVG/(CP025(J)+AMOD(RNAVG,60.)/60.*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF GREATER THAN 240 DEGREES OF FREEDOM, USE THE FORMULA LISTED ON
C  P.388, BENDAT & PIERSOL, 1971 
C
  220 Zalpha=1.96
      TERM=2./9./RNAVG
      TERMSQ=SQRT(TERM)
      CUPPER=RNAVG/(RNAVG*(1.0-TERM-ZALPHA*TERMSQ)**3)
      CLOWER=RNAVG/(RNAVG*(1.0-TERM+ZALPHA*TERMSQ)**3)
  300 CONTINUE
      RETURN
      END
      subroutine atan3(u,v,ampl,alpha)
c
c --- This subroutine calculates the angle counted clockwise form 
c     true north.
c
      radpdeg = acos(-1.0)/180.
      degprad = 1./radpdeg
      alpha = atan2(v,u) * degprad
      if(alpha.le.0.) then
        alpha = 90.- alpha
      else
        if(alpha.lt.90) then
          alpha = 90.- alpha
        else
          alpha = 450.- alpha
        end if
      end if
      ampl = sqrt(u*u+v*v)
      return
      end
                                                                     

C*******************************************************************************
!      PROGRAM LSQHA3
      SUBROUTINE LSQHA3
C
C  REVISED AND SIMPLIFIED BY CHRIS ZERVAS   DECEMBER 1995
C
C  REVISED AND DOCUMENTED 1983
C  USES FORTRAN77
C  PROGRAM MODIFIED TO HANDLE TIDES AND TIDAL CURRENTS: JULY, 1985
C  PROGRAM MODIFIED JULY, 1985 BY E. E. LONG
C  COMPUTES (U) AND (V) VECTOR TIDAL CURRENTS ON ONE PASS
C  ORIGINAL LSQHA2 PROGRAM SAVED INTACT
C  -------------------------------------------------------------------
C
C    UNIVERSITY OF FLORIDA  17 FEBRUARY 1984
C
C    THIS IS THE LATEST VERSION OF LSQHA2
C     IT CONTAINS A NEW VERSION OF SUBROUTINE -SCREEN-
C     THIS VERSION IS COPIED TO J FANCHER/NOAA ON THIS DATE
C____________________________________________________________________
C
C
C  PURPOSE:
C     PROGRAM LSQHA2, A REVISION OF LSQHA, EXPANDS A GIVEN TABULAR FUNC-
C  TION IN THE FORM  Y(M) = SUM OF A(N)COS(F(N)M-M1), WHERE F(N) IS AN
C  ARBITRARY SET OF FREQUENCIES.  THE SOLUTION IS A BEST FIT IN A LEAST
C  SQUARES SENSE.  LSQHA2 IS SPECIFICALLY DESIGNED FOR THE ANALYSIS OF
C  HOURLY TIDE DATA IN TERMS OF THE 37 TRIGONOMETRIC FUNCTIONS TRADI-
C  TIONALLY USED BY THE N.O.S..  THE TIME INTERVAL MAY ALSO BE CHANGED
C  BY THE USER, OR THE STANDARD LIST OF FREQUENCIES MAY BE REPLACED WITH
C  OTHER FREQUENCIES BY THE USER.  NO PROGRAM CHANGES ARE REQUIRED WHEN
C  USING FEWER THAN 37 CONSTITUENTS, BUT WHEN FEWER CONSTITUENTS ARE
C  USED, SPACE MAY BE SAVED BY REWRITING THE DIMENSION STATEMENTS IN THE
C  MAIN PROGRAM AND IN SUBROUTINE SCREEN.  REVISION OF THESE DIMENSION
C  STATEMENTS, AND ENTERING A NEW LIST OF FREQUENCIES AND IDENTIFIERS
C  ARE THE ONLY CHANGES REQUIRED FOR INCREASING THE NUMBER OF
C  FREQUENCIES.
C     SUBROUTINE --ENT--  IS USED TO CONSTRUCT THE DATA FIELDS BY USING
C  THE FOLLOWING VARIABLES:
C
C     NOBS    = NUMBER OF OBSERVATIONS
C     MORE    = ZERO, THE FINAL DATA IS BEING READ FOR ANY RUN
C             = 1, ADDITIONAL DATA ARE TO BE READ LATER
C     JYY     = YEAR, 4 DIGITS                X
C     JMM     = MONTH                          X  FIRST
C     JDD     = DAY                             > DATA
C     JHH     = HOUR                           X  SAMPLE
C     JMIN    = MINUTE (IN TENTHS OF HOURS)   X
C     TH(I)   = TIDAL HEIGHTS, I=1,NOBS
C     THBAR   = APPROXIMATE MEAN WATER LEVEL. MUST BE A
C               CONSTANT FOR ANY MACHINE RUN
C     XLONG   = LONGITUDE OF STATION (+ FOR WEST, - FOR EAST),
C               TO HUNDREDTHS OF A DEGREE.
C     TIMEZ   = TIME ZONE OF STATION.  REMEMBER, + IS ASSUMED WHEN
C               BLANK, BUT, - MUST BE USED WHEN REQUIRED.
C
C     ANY OTHER PROCEDURE FOR ESTABLISHING THIS DATA FIELD IS
C  ACCEPTABLE.    ALL OTHER SUBROUTINES ARE HIGHLY STRUCTURED
C  AND SHOULD NOT BE CHANGED WITHOUT VERY CAREFUL CONSIDERATION.
C     PROGRAM LSQHA2 IS USED FOR ALL INPUT OTHER THAN THE DATA TO BE
C  ANALYZED, AND FOR THE INITIALIZATION OF FIELDS USED BY TWO OR MORE
C  SUBROUTINES.
C     PROGRAM LSQHA2 SHOULD ALSO BE USED FOR ANY REQUIRED LINKAGE TO
C  OTHER PROGRAMS, WITH ONE EXCEPTION.  ALL FINAL COMPUTATIONAL RESULTS
C  ARE ACCESSIBLE THROUGH PROGRAM LSQHA2 BY MEANS OF COMMON STATEMENTS
C  OR CALL STATEMENTS.  THE EXCEPTION IS DESCRIBED IN CONNECTION WITH
C  SUBROUTINE -SCREEN-.
C     TWO TYPES OF CONTROL STATEMENTS ARE USED, VARIABLE AND BINARY.
C  1) VARIABLE CONTROL PARAMETERS ARE LISTED IN INPUT ORDER
C     AS FOLLOWS:
C     NCONST   = NUMBER OF CONSTITUENTS TO BE CONSIDERED
C     NYYYY    = YEAR NUMBER OF THE DATA BEING ANALYZED.  THIS PARAMETER
C                CONTROLS THE ADJUSTMENT FOR THE EQUILIBRIUM ARGUMENTS
C                AND NODE FACTORS.  (0000 MEANS NO ADJUSTMENT.)
C     NDELTT   = NUMBER OF SAMPLES PER HOUR
C     CUTOFF   = CALCULATIONS ARE TERMINATED WHEN THE NEXT PREDICTOR
C                SELECTED ACCOUNTS FOR A FRACTION OF THE TOTAL VARIANCE
C                LESS THAN -CUTOFF-.
C     MINTRM   = MINIMUM NUMBER OF   :   (BOTH ARE USED IN
C                TERMS WANTED.       :   -SCREEN- FOR PRINTING
C     MAXTRM   = MAXIMUM NUMBER OF   :   PORTIONS OF THE
C                TERMS WANTED.       :   OUTPUT EQUATION.)
C
C  2) THE BINARY CONTROL FACTORS ARE LISTED BELOW.  IN EACH CASE, A ZERO
C  INDICATES THAT THE OPTIONAL ACTION IS NOT WANTED.  A ONE INDICATES
C  THAT IT IS WANTED.
C
C        . . . PROGRAM OPTIONS . . .
C
C     ICNTL(1)  = NOT USED
C     ICNTL(2)  = NOT USED
C     ICNTL(3)  = NOT USED
C     ICNTL(4)  = PRINT TABLE OF MEANS AND STANDARD DEVIATIONS FROM
C                    SUBROUTINE -CSTAT2-
C     ICNTL(5)  = USE MATRIX OF CORRELATION COEFFICIENTS INSTEAD OF
C                    COVARIANCE TO IMPROVE STABILITY
C     ICNTL(6)  = NOT USED
C     ICNTL(7)  = DISPLAY INTERMEDIATE RESULTS FROM SUBROUTINES -SCREEN-
C                    AND -CSTAT2-
C     ICNTL(8)  = NOT USED
C     ICNTL(9)  = CONTINUE THROUGH LIST OF ALL PREDICTORS INSTEAD OF
C                    TERMINATING CALCULATIONS AFTER FINDING FIRST
C                    PREDICTOR WHOSE VARIANCE IS BELOW CUTOFF
C     ICNTL(10) = CALCULATE PREDICTORS IN ORDER GIVEN INSTEAD OF REARRANGING
C                    MATRIX TO CHOOSE PREDICTORS BASED ON THE LARGEST
C                    REDUCTION OF VARIANCE
C     AT THIS TIME -ICNTL- IS DIMENSIONED TO HAVE UP TO 10 OPTIONS.
C        C     OTHER VARIABLE NAMES USED IN -LSQHA2-                                     C
C     FREQ(I)   = CHANGE IN THETA IN ONE TIME INCREMENT.
C  + LABEL(I)   = ARRAY CONTAINING STANDARD SYMBOLS FOR TIDAL CONSTIT-
C                 UENTS OR OTHER DATA TYPE IDENTIFIER USED FOR OUTPUT.
C     NDEXY     = UPPER BOUND OF DATA MATRIX + 1; THIS COLUMN CONTAINS
C                 (TH(I) - THBAR) FOR STANDARD HARMONIC ANALYSIS.  WHEN
C                 ADDITIONAL TYPES OF INPUT DATA ARE USED THIS COLUMN
C                 CONTAINS ...................
C     NDEXM1    = UPPER BOUND OF DATA MATRIX
C     NDEXP1    = (NDEXY + 1); THIS COLUMN CONTAINS THE CONSTITUENT
C                 INDEX NUMBER.
C     NDEXP2    = (NDEXY + 2); THIS COLUMN CONTAINS A "2" INDICATING
C                 AMPLITUDE PHASE PAIRS OF HARMONIC CONSTITUENTS, OR
C                 "1" INDICATING NON PERIODIC PREDICTOR AMPLITUDES
C                 WITH NO ASSOCIATED PHASE.
C     NDR       = THE NUMBER OF DATA SETS READ.  MORE = 1, INDICATES
C                 THAT MORE DATA ARE TO BE READ.  -NDR- WILL
C                 INCREMENT UNTIL ALL DATA SETS ARE READ, I.E.,
C                 UNTIL -MORE- IS MADE EQUAL TO ZERO.  IF MORE = 0
C                 AT THE START, NDR = 1.
C     NRECRD    = NUMBER OF OBSERVATIONS OR TIDAL HEIGHTS USED IN THE
C                 ANALYSIS; EQUIVALENT TO -NOBS- IF THERE ARE NO
C                 BREAKS IN THE DATA SERIES.
C     NV        = NUMBER OF SINE AND COSINE VECTORS FOR PERIODIC
C                 CONSTITIENTS; (2 * NCONST).
C     SETS      = TOTAL NUMBER OF OBSERVATIONS (DATA POINTS) IN SERIES
C                 BEING ANALYZED.  ONLY IF DATA SERIES IS ONE UNBROKEN
C                 SET IS -SETS- EQUAL TO -NRECRD-.  IF SERIES TO BE
C                 ANALYZED HAS BREAKS, -SETS- IS LESS THAN -NRECRD-.
C  +  SPEED(I)  = CONSTITUENT FREQUENCY IN DEGREES-PER-SOLAR-HOUR.
C     YEPOCH(I) = EQUILIBRIUM ARGUMENTS TO MODIFY -EPOCH- TO A STANDARD
C                 YEAR.
C     YNODE(I)  = NODE FACTORS TO MODIFY -HSUBN- TO A STANDARD YEAR.
C
C     (ADD THEM HERE)
C
C      .  .  .  DIMENSION AND DECLARATION STATEMENTS .  .  .
C
C **********************************************************************
C  M1 THRU M6 MUST BE REPLACED BY NUMERICAL VALUES IN PROGRAM LSQHA2 AND
C  IN THE DIMENSION STATEMENT IN SUBROUTINE -SCREEN- BEFORE CALCULATIONS
C  ARE ATTEMPTED.
C     M1 = NCONST + 1               M3 = 0                  M5 = 2 * M2
C     M2 = 2 * NCONST + 3           M4 = 12 (USUALLY)       M6 = NCONST
C **********************************************************************
      INTEGER M1,M2,M3,M4,M5,M6
c     PARAMETER ( M1=146, M2=293, M3=0000, M4=20, M5=586, M6=145 )
      PARAMETER ( M1=176, M2=353, M3=0000, M4=20, M5=706, M6=175 )
c     PARAMETER ( M1=38, M2=77, M3=2560, M4=12, M5=154, M6=37)
      PARAMETER ( NRECRD= 200000 )
      REAL*8 FREQ(M1),FHOUR(M4),TIM(NRECRD),SPEED,DELTT,SHOUR,DELTAX
      REAL   PP(M2,M2),QQ(M2,3)
      REAL   PE(M2,M2),QE(M2,3)
      REAL   TH(NRECRD),SIG(M2),R(M2),
     *     TSTORE(M1),SX(M2),V(M2)
      REAL   TE(NRECRD),SIE(M2),RE(M2),
     *     TSTORN(M1),       SY(M2),VE(M2)
      REAL AEPOCH(M1),AHSUBN(M1),COLA(M1),EPOCH(M1),HSUBN(M1),
     *     KAPPA(M1),RT(M2),R1(M2),RU(M2),R2(M2)
C      REAL TIME(15)
      CHARACTER*10 LABEL,IDENS,TKIND,SKIND,LABLE
      CHARACTER*4 BANAME(8)*4
      CHARACTER*80 FNAME
      DIMENSION JR(M2),JRE(M2)
      DIMENSION VIN(NRECRD),DIN(NRECRD)
      DIMENSION KAYE(M2),MCNSTE(M2),MORDER(M2),NQPE(M2),INOBS(M4)
      DIMENSION ICNTL(10),KAY(M2),MCNST(M2),NORDER(M2),NQP(M2)
      EQUIVALENCE(EPOCH(1),TIM(1)),(HSUBN(1),TIM(200))
      EQUIVALENCE(KAPPA(1),TIM(400)),(COLA(1),TIM(600))
      EQUIVALENCE(AEPOCH(1),TIM(800)),(AHSUBN(1),TIM(1000))
      EQUIVALENCE(R(1),TIM(1200)),(RE(1),TIM(1600))
      EQUIVALENCE(RT(1),TIM(2000)),(RU(1),TIM(2400))
      EQUIVALENCE(R1(1),TIM(2800)),(R2(1),TIM(3200))
      EQUIVALENCE(SIG(1),TIM(3600)),(SIE(1),TIM(4000))
      EQUIVALENCE(TH(1),VIN(1)),(TE(1),DIN(1))
      EQUIVALENCE(QQ(1,1),TIM(4400)),(QE(1,1),TIM(5600))
      EQUIVALENCE(TSTORN(1),TIM(6800)),(TSTORE(1),TIM(7000))
      EQUIVALENCE(JRE(1),JR(1)),(NQPE(1),NQP(1)),(MORDER(1),NORDER(1))
      EQUIVALENCE(MCNSTE(1),MCNST(1)),(KAYE(1),KAY(1))
C
C
C      . . . COMMON STATEMENTS . . .
C
      COMMON/PARAM1/NCONST,NYYYY,NDELTT,CUTOFF,TZERO,MINTRM,MAXT
     1RM,NCOL,NROW,TEBAR
      COMMON/PARAM2/NOBS,MORE,JYY,JMM,JDD,JHH,JMIN,THBAR,XLONG,TIMEZ
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(NRECRD,2)
     1 ,RATIO                                             
      COMMON/PARAM4/IYY(20),IMM(20),IDD(20),IHH(20),IMIN(20)
      COMMON/PARAM5/YEPOCH(180),YNODE(180),SPEED(180),P(180)
      COMMON/PARAM6/IDENS(16),LABEL(180),TKIND,SKIND
      COMMON/BGROUP/DELTT,SHOUR(20),DELTAX
      COMMON/MEANL1/MNYEAR,JYEAR,JDAZ,DAY,LDAYY,LLDAYY,MONTH
      COMMON/MEANL2/TIMEX,JDAAYI,JDAYI,JDAYY,LDAA,LXYERE
      COMMON/MEANL3/IALL,LPYER,LYDAY,LJDAA,LNDAY,LGYER,LZDAY
      COMMON/MEANL4/NJOBX,MMONTH,NXDAY,JDAZZ,JDAAZ,FST,INDATA
      COMMON/PARAMC/RAMP(180),RKAPP(180),LABLE(180),RKAP(180)
      NAMELIST /NML0/FNAME,CUTOFF,NBLK
      NAMELIST /NICNTL/ICNTL
C
C     .  .  .  DATA STATEMENTS  .  .  .
C
      DATA (BANAME(I),I=1,8)/'V-12', 'U-12', 'V-18', 'U-18',
     * 'V-24', 'U-24', 'V-30', 'U-30'/
C
C  .  .  .  PROGRAM BEGINS HERE  .  .  .
C
C * * * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * *
C
C  SET UP CALL TO TSKTME TO MEASURE TIMES DURING EXECUTION -
C  STRICTLY FOR DEBUGGING PURPOSES.
C
C      CALL TSKTME
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***********************************************************************
C
C     . . . STEP 1 . . .
C
C     READ AND ECHO PRINT PROGRAM CONTROL FIELDS
C
C***********************************************************************
C
c     READ(5,NML=NML0)
c     WRITE(6,NML=NML0)
      READ(25,'(a80)')FNAME
      READ(25,*)CUTOFF,NBLK,NJOBX
c     READ(5,NML=NICNTL)
c     WRITE(6,NML=NICNTL)
!  330 OPEN(UNIT=9,FILE=FNAME,FORM='FORMATTED',
!     * ACCESS='SEQUENTIAL')
      DO 6005 I=1,16
      IDENS(I)='          '
 6005 CONTINUE
       IDENS(1) = 'Harmonic A'
       IDENS(2) = 'nalysis of'
!       IDENS(3) = ' Data in  '
       IDENS(3) = FNAME(1:10)
       IDENS(4) = FNAME(11:20)
       IDENS(5) = FNAME(21:30)
       IDENS(6) = FNAME(31:40)
!       IDENS(7) = FNAME(31:40)
      write(idens(7),171) RATIO
  171 format(3H R=,F6.3,X)
      DO 333 I = 1,10
  333 ICNTL(I) = 0
      TZERO = CUTOFF
      NROW =   M6
      NCOL = NROW
      NDR = 0
      avgazi = 0.
      MORE = 1
      DO 250 IBLK = 1,NBLK
      IF(IBLK.EQ.NBLK) MORE = 0
  444 CALL ENT(NRECRD,NDR,TH,TE,TIM,FREQ,VIN,DIN)
      avgazi = avgazi + azi
      IYY(NDR) = JYY
      IMM(NDR) = JMM
      IDD(NDR) = JDD
      IHH(NDR) = JHH
      IMIN(NDR) = JMIN
      IF(NDR.GT.1) GO TO 240
C********************************************************************
C
      MAXTRM = 0
      MINTRM = 0
c     MAXTRM = 10
c     MINTRM = 1
C
C********************************************************************
C______________________________________________________________________
C
C     TEMPORARY ALTERATION
C
      PRINT 12, NCONST,NYYYY,NDELTT,CUTOFF,TZERO,MINTRM,
     * MAXTRM,NCOL,NROW
   12 FORMAT ( ' NCONST=',I4,' NYYYY=',I5,
     * ' NDELTT=',I3,' CUTOFF=',E10.2,' TZERO=',E10.2,
     * ' MINTRM=',I3,' MAXTRM=',I3,' NCOL=',I4,' NROW=',I4,/)
C______________________________________________________________________
      DELTT = 1.0/FLOAT(NDELTT)
      DELTAX = DFLOAT(NDELTT)
      PRINT 16, NCONST,DELTT,CUTOFF,NYYYY
   16 FORMAT(' Harmonic Analysis considering',I4,' constituents.',/,
     * ' The time interval between observations is',F5.2,' hours.',
     * '  Predictors which account for less than',F10.7,' of the',
     * ' variance',/,' of the input data series are dropped.',
     * '  Harmonic constants are adjusted for the year',I5,/)
      PRINT 22, (ICNTL(I),I=1,10)
   22 FORMAT (' Binary control fields:',10I2 /)
C
C***********************************************************************
C
C
C     ESTABLISH CONTROL CONSTANTS  .  .  .
C
 183  NV = 2 * NCONST
      NDEXY = NV + 1
      I1 = 1
      ITEST = 0
      IF( NDEXY .LE. M2 ) GO TO 199
      PRINT 195, NDEXY, M1
  195 FORMAT (1X,' THE NUMBER OF VARIABLES',I4,' EXCEEDS THE ALLOTTED ',
     * 'ORDER OF THE CORRELATION MATRIX',I4,'.  THE TOTAL NUMBER OF ',
     * 'VARIABLES HAS BEEN',/,' REDUCED TO FIT THE ALLOTTED SPACE.',/)
      NDEXY = M2
      IF(NV .LT. NDEXY) GO TO 199
      PRINT 197, M1, NCONST
 197  FORMAT(1X,'THE ALLOTTED SIZE OF THE MATRIX',I4,' IS TOO SMALL FOR 
     *THE ',I3,' CONSTITUENTS REQUESTED... CALCULATIONS TERMINATED.'
     *,/)
      GO TO 777
C
C     ZERO ALL SUMS, INITIALIZE COUNTERS AND INDICES
C
  199 SETS = 0.0
      DO 130 I=1,NDEXY
      SY(I) = 0.0
      SX(I)=0.0
      NORDER(I) = 0
      MORDER(I) = 0
      MCNST(I) = 0
      DO 125 J=1,NDEXY
      PE(I,J) = 0.0
  125 PP(I,J) = 0.0
  130 CONTINUE
      NDEXP1 = NDEXY + 1
      NDEXP2 = NDEXY + 2
      NDEXM1 = NDEXY - 1
C
C     LABEL ROWS AND COLUMNS IN MATRIX
      DO 135 L=1,NDEXY
      PE(L,NDEXP1) = L
  135 PP(L,NDEXP1) = L
      DO 235 L=1,NV
      PE(L,NDEXP2) = 2
 235  PP(L,NDEXP2) = 2
      DO 236 L=NV+1,NDEXM1
      PE(L,NDEXP2) = 1
 236  PP(L,NDEXP2) = 1
C
C***********************************************************************
C
  240 INOBS(NDR) = NOBS
C      TIME(4) = TTIME(0)
C  ACCUMULATE TOTAL NUMBER OF OBSERVATIONS.
      SETS = SETS + FLOAT(NOBS)
C
C
C      TIME(5) = TTIME(0)
      IF(KAZI.EQ.2) GO TO 245
  242 CALL MODT1(PP,FREQ,FHOUR,SX,TH,THBAR,V,I1,
     * M2,NDEXY,NOBS,NV,TIM,NDR)
      IF(INDATA.EQ.2) CALL TRGSARAN(PP,FREQ,SX,NOBS,M2,NV,TIM)
  245 IF(KINDAT.EQ.2) GO TO 250
      IF(KAZI.EQ.1) GO TO 250
      CALL MODT1(PE,FREQ,FHOUR,SY,TE,TEBAR,VE,I1,
     * M2,NDEXY,NOBS,NV,TIM,NDR)
      IF(INDATA.EQ.2) CALL TRGSARAN(PE,FREQ,SY,NOBS,M2,NV,TIM)
  250 CONTINUE
C      TIME(6) = TTIME(0)
C***********************************************************************
C
C     . . . STEP 8 . . .
C
C***********************************************************************
      IF(INDATA.EQ.2) GO TO 260
      IF(KAZI.EQ.2) GO TO 255
      CALL TRGSA ( PP,FREQ,FHOUR,SX,NDR,INOBS,M2,     NV)
  255 IF(KINDAT.EQ.2) GO TO 260
      IF(KAZI.EQ.1) GO TO 260
      CALL TRGSA (PE,FREQ,FHOUR,SY,NDR,INOBS,M2,     NV)
C      TIME(7) = TTIME(0)
  260 IF(KAZI.EQ.2) GO TO 265
      CALL CSTAT2 (PP,OMEAN,OSD,OVAR,R,SIG,SX,SETS,THBAR,TZERO,
     * ICNTL,ITEST,KAY,M2,M4,       NDEXY,       NV,NCOL)
C      TIME(8) = TTIME(0)
      IF(ITEST.EQ.1) GO TO 290
  265 IF(KINDAT.EQ.2) GO TO 270
      IF(KAZI.EQ.1) GO TO 270
      CALL CSTAT2 ( PE,PMEAN,PSD,PVAR,RE,SIE,SY,SETS,TEBAR,TZERO,
     * ICNTL,ITEST2,KAYE,M2,M4,       NDEXY,       NV,NCOL )
      IF(ITEST2.EQ.1) GO TO 290
C***********************************************************************
C
C     . . . STEP 9 . . .
C
C***********************************************************************
  270 NECOM = KAZI
      IF(NECOM.EQ.3) KAZI = 1
      IF(KAZI.EQ.2) GO TO 275
      avgazi = avgazi/ndr
      if(kindat.eq.1)then
      write(idens(15),271) nint(avgazi)
  271 format(6Halong ,I3,X)
      idens(16)='degrees   '
      end if
      CALL HEADER(NECOM,KAZI,AZI,KINDAT,TKIND,SKIND)
      CALL SCREEN (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * CUTOFF,EPOCH,HSUBN,ICNTL,                    INOBS,    JR,
     * K,KAPPA,MAXTRM,MCNST,MINTRM,   M2,     NCOL,NCONST,NDEXM1,
     * NDEXP1,NDEXP2,NDEXY,NDR,NORDER,NQP,     NV,NYYYY,OMEAN,OSD,
     * OVAR,  QQ,RT,R1,SETS,SIG,SX,THBAR,TIMEZ,TSTORE,TZERO,
     * XLONG )
      open(10,file='cons.out',status='unknown')
      call cards(rkapp,ramp,omean,ndeltt,1.,idens,4,avgazi)
C     TIME(11) = TTIME(0)
  275 IF(KINDAT.EQ.2) GO TO 290
      IF(NECOM.EQ.3) KAZI = 2
      IF(KAZI.EQ.1) GO TO 290
      if(kindat.eq.1)then
      azi90 = avgazi + 90.
      if(azi90.ge.360.) azi90 = azi90 - 360.
      write(idens(15),271) nint(azi90)
      idens(16)='degrees   '
      end if
      CALL HEADER(NECOM,KAZI,AZI,KINDAT,TKIND,SKIND)
      DO 277 III = 1,NDEXY
      MCNSTE(III)=0
  277 MORDER(III) = 0
      CALL SCREEN (PE,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * CUTOFF,EPOCH,HSUBN,ICNTL,                    INOBS,    JRE,
     * KK,KAPPA,MAXTRM,MCNSTE,MINTRM,   M2,      NCOL,NCONST,NDEXM1,
     *NDEXP1,NDEXP2,NDEXY,NDR,MORDER,NQPE,     NV,NYYYY,PMEAN,PSD,
     * PVAR,   QE,RU,R2,SETS,SIE,SY,TEBAR,TIMEZ,TSTORN,TZERO,
     * XLONG)
      call cards(rkapp,ramp,pmean,ndeltt,1.,idens,4,azi90)
C************************************************************************
C
C     . . . STEP 10 . . .
C
C*****************************************************************************
C     CALL ONE80
C
C     ANY LOCAL PROGRAMS TO SAVE RESULTS IN COMPUTER COMPATIBLE FORMAT.
C     OR TO USE THE HARMONIC CONSTANTS FOR ADDITIONAL CALCULATIONS
C     SHOULD BE CALLED HERE.
C
C
C     ALL FILES ESTABLISHED BY -LSQHA2- SHOULD BE CLOSED.
C
c 280 PRINT 285
  285 FORMAT (5X,'END OF OUTPUT FOR THIS JOB . . . . .THANK YOU')
 290  CONTINUE
C      TIME(12) = TTIME(0)
C      PRINT 770, (I,TIME(I),I=1,12)
C  770 FORMAT (1H0,'TIME(',I2,')=',F8.4)
  777 CONTINUE
  773 CLOSE(9)
      STOP
      END
      SUBROUTINE ASTRO(XYER,DAYB,DAYM,GRBS,GRMS,JOBX)                   
      IMPLICIT REAL*4(A-H,O-Z)
      DIMENSION CXX(30),OEX(5)                                          
      COMMON/BOXA/S,XL,PM,PL,SL,PS,PLM,SKYN,VI,V,XI,VPP                 
      COMMON/BOXB/VP,P,AUL,AUM,CRA,CQA                                  
      COMMON/VEE/TML,CON,U,Q,UI                                         
      COMMON/BOXS/AW,AI,AE,AE1,ASP                                      
    6 FORMAT( //   5X,  7HYEAR = , I5  // )                             
 2110 FORMAT(1H0, 38HERROR.....PARAMETER (I) EXCEEDS LIMITS )           
 2111 FORMAT(1H0,15X,1HR,10X,1HQ,15X,9HU OF M(2))                       
 2112 FORMAT( 8X, F10.3, 3X, F10.3, 9X, F10.4)                          
 2040 FORMAT(1H ,10F10.3//// 1X,10F10.3//// 1X,10F10.3)                 
 2141 FORMAT(10X /// 5X, 6HDAYB =, F6.0, 5X, 6HDAYM =, F6.0 )           
 2142 FORMAT(10X /// 5X, 6HGRBS =, F7.2, 5X, 6HGRMS =, F7.2)            
      PINV = 57.29578
      NYEAR = IFIX(XYER)                                                
      CALL ORBIT(XCEN,XSX,XPX,XHX,XP1X,XNX,OEX,T,XYER,5)                
      XW = 23.4522944 - .0130125*T - .00000164*T**2 + .000000503*T**3
      XI = 5.14537628
      AW = XW*0.0174533
      AI = XI*0.0174533
      AE = 0.0548997
      AE1 = 0.01675104 - 0.0000418*T - .000000126*T**2
      ASP = .46022931
      DO 30 NOE = 1,30                                                  
   30 CXX(NOE) = 0.0                                                    
      IF(DAYB.GT.0.0) DAYB = DAYB - 1.0
      IF(DAYM.GT.0.0) DAYM = DAYM - 1.0
      DOBY = 0.0
      AMIT = 0.0
      AMI = XYER - XCEN
      CPLEX = XCEN/400.0 + 0.0001
      DICF = CPLEX - AINT(CPLEX)
      IF(AMI.EQ.0.0) GO TO 32
      XCET = XCEN + 1.0
      CDIF = XYER - XCET
      DOBY = CDIF/4.0 + 0.0001
      AMIT = AINT(DOBY)
      IF(DICF.LT.0.001) AMIT = AMIT + 1.0
   32 FARM = 0.25*AMI
      FARX = FARM - AINT(FARM)                                              
      CXX(1) = XSX + 129.384820*AMI + 13.1763968*(DAYB + AMIT) + 0.54901
     16532*GRBS
      CXX(2) = XPX + 40.6624658*AMI + 0.111404016*(DAYB + AMIT) + 0.0046
     141834*GRBS
      CXX(3) = XHX - 0.238724988*AMI + 0.985647329*(DAYB + AMIT) + 0.041
     1068639*GRBS
      CXX(4) = XP1X + 0.01717836*AMI + 0.000047064*(DAYB + AMIT) + 0.000
     1001961*GRBS
      CXX(5) = XPX + 40.6624658*AMI + 0.111404016*(DAYM + AMIT) + 0.0046
     141834*GRMS
      CXX(6) = XNX - 19.3281858*AMI - 0.052953934*(DAYM + AMIT) - 0.0022
     106414*GRMS
   40 CXX(7) = XPX + 40.6624658*AMI + 0.111404016*(DAYM + AMIT) + 0.0046
     141834*GRBS
      CXX(8) = XNX - 19.328185764*AMI - 0.0529539336*(DAYB + AMIT) - 0.0   
     1022064*GRBS
      CALL TWOPI(CXX, 8)                                                   
   41 DO 100 II = 1,8                                                      
  100 CXX(II) = FLOAT(IFIX(CXX(II)*100.0 + 0.5))*0.01                      
      ANG = CXX(8)                                                         
      CALL TABLE6(VIB,VB,XIB,VPB,VPPB,XX,XX,XX,XX,XX,ANG,ANB,ATB)
      CXX(26) = VIB                                                        
      CXX(27) = VB                                                         
      CXX(28) = XIB                                                        
      CXX(29) = VPB                                                        
      CXX(30) = VPPB                                                       
      CXSB = CXX(1)
      CXPB = CXX(2)
      CXHB = CXX(3)
      CXP1B= CXX(4)
      ANG = CXX(6)                                                         
      CALL TABLE6(VI,V,XI,VP,VPP,CIG,CVX,CEX,PVC,PVCP,ANG,AN,AT)           
      CXX(9 ) = VI                                                         
      CXX(10) = V                                                          
      CXX(11) = XI                                                         
      CXX(12) = VP                                                         
      CXX(13) = VPP                                                        
  230 DO 333 II = 9,13                                                     
  333 CXX(II) = FLOAT(IFIX(CXX(II)*100.0 + 0.5))*0.01                      
      PGX = CXX(5) - CXX(11)                                               
      PGX = FLOAT(IFIX(PGX*100.0 + 0.5))*0.01                              
      CALL TWOPI(PGX, 1)                                                  
      XPG = PGX*0.0174533                                                 
      CXX(14) = PGX                                                       
      RAXE = SIN(2.0*XPG)                                                 
      RAXN = (COS(0.5*AT)**2/(6.0*SIN(0.5*AT)**2)) - COS(2.0*XPG)
      RXX = 0.0
      IF(RAXE.EQ.0.0.OR.RAXN.EQ.0.0) GO TO 232                            
      RAX = RAXE/RAXN                                                     
      IF(RAX.GT.3450.0) GO TO 232                                         
        RXX   = ATAN(RAX )*PINV
      CXX(22) = RXX                                                       
  232 CRA = SQRT(1.0 - 12.0*(SIN(0.5*AT)**2/COS(0.5*AT)**2)*COS(2.0*XPG)
     1 + 36.0*(SIN(0.5*AT)**4/COS(0.5*AT)**4))
      UM2 = 2.0*(CXX(11) - CXX(10))                                       
      CXX(21) = UM2                                                       
      CXX(24) = CRA                                                       
      UL2 = UM2 - RXX
  404 UL2 = UL2 + 180.0
  405 CXX(15) = UL2
      ZES = (5.0*COS(AW) - 1.0)*SIN(XPG)                                  
      ZEC = (7.0*COS(AW) + 1.0)*COS(XPG)                                  
      CALL FITAN(ZES,ZEC,QXX,SPXX,2)                                      
      CXX(23) = QXX                                                       
      CRAV = 0.5*UM2 + QXX + 090.0
      CXX(16) = CRAV                                                      
      CQA = SQRT(0.25 + 1.5*((COS(AW)/COS(0.5*AW)**2)*COS(2.0*XPG)) + 2.
     125*(COS(AW)**2/COS(0.5*AW)**4))
      CXX(25) = CQA                                                       
      DO 444 III = 14,23                                                  
  444 CXX(III) = FLOAT(IFIX(CXX(III)*100.0 + 0.5 ))*0.01                  
      IF(JOBX.GT.15.AND.JOBX.GT.0 ) GO TO 88                              
      PGXX =  CXX(7) - CXX(11)                                            
      CALL TWOPI(PGXX,1)                                                  
      XXPG = PGXX*0.0174533                                               
      CXX(17) = PGXX                                                      
      BATX = CXX(8)*0.0174533                                             
      EYEX = COS(AI)*COS(AW) - SIN(AI)*SIN(AW)*COS(BATX)                  
      CXX(20) = ACOS(EYEX)*PINV
      UM2X = 2.0*(CXX(11) - CXX(10))                                        
      EYIT = CXX(20)*0.0174533                                              
      ZEXS = (5.0*COS(AW) - 1.0)*SIN(XXPG)                                  
      ZEXC = (7.0*COS(AW) + 1.0)*COS(XXPG)                                  
      CALL FITAN(ZEXS,ZEXC,QXXX,SPXXX,2)                                    
      CXX(26) = QXXX                                                        
      CRAVX = 0.5*UM2X + QXXX + 090.0
      CXX(19) = CRAVX                                                       
      CQA = SQRT(0.25 + 1.5*((COS(AW)/COS(0.5*AW)**2)*COS(2.0*XXPG)) +2.    
     125*( COS(AW)**2/COS(0.5*AW)**4))                                      
      CXX(25) = CQA                                                         
   88 CONTINUE                                                              
      PM   = CXX(1)                                                         
      PL   = CXX(2)                                                         
      SL   = CXX(3)                                                         
      PS   = CXX(4)                                                         
      PLM  = CXX(5)                                                         
      SKYN = CXX(6)                                                         
      VI   = CXX(9)                                                         
      V    = CXX(10)                                                        
      XI   = CXX(11)                                                        
      VP   = CXX(12)                                                        
      VPP  = CXX(13)                                                        
      P    = CXX(14)                                                        
      AUL  = CXX(15)                                                        
      AUM  = CXX(16)                                                        
      CRA = CXX(24)                                                         
      CQA = CXX(25)                                                         
      U = V*0.0174533                                                       
      Q = P*0.0174533                                                       
      UI = VI*0.0174533                                                     
      RETURN                                                                
      END                                                                   
C
      SUBROUTINE AZIM(DX,VR,A,VN,VE,NEC,J,COMPV)                            
      DIMENSION DX(J), VR(J),VN(J),VE(J)                                    
      DO 100 I = 1,J                                                        
      DIR = (DX(I) + COMPV - A)*0.0174533
      VEL = VR(I)
      VE(I) = 0.0
      VN(I) = 0.0
      GO TO (10,20,10,10,10),NEC
   10 VN(I) = VEL*COS(DIR)
      GO TO (100,100,20,20,100),NEC
   20 VE(I) = VEL*SIN(DIR)
  100 CONTINUE
      RETURN                                                                
      END                                                                   
      SUBROUTINE CARDS(W,D,AVE,NSPH,TFAC,IDENS,K,AZI)
      CHARACTER*16 IDENS(10)
      DIMENSION W(180),D(180),IW(37),ID(37)
      write(10,1)(IDENS(I),I=1,10)
      IAVE = nint((AVE) * 1000.)
      IAZI = AZI
      NOS=0
      IF(K.EQ.4) NOS=1
      IF(K.EQ.3.AND.IAZI.EQ.0) NOS=1
      IF(K.EQ.3.AND.IAZI.NE.0) NOS=2
C     write(10,2)IAVE,NSPH,NOS,TFAC
      write(10,2)IAVE
C
      DO 5 I=1,37
      IW(I) = nint(W(I)*10.)
      ID(I) = nint((D(I))*1000.)
    5 CONTINUE
      do 11 n=1,5
      nn=7*(n-1)
      write(10,6)N,(ID(nn+j),IW(nn+j),j=1,7)
   11 continue
      write(10,6)6,ID(36),IW(36),ID(37),IW(37)
C
      IEBB=IAZI+180
      if(iebb.gt.360)iebb=iebb-360
      IF(NOS.NE.1) write(10,7)IAZI,IEBB
      RETURN
    1 FORMAT(5A16)
    2 FORMAT(I6,I2,I2,2H 0,F10.8)
    6 FORMAT(7X,I1,7(I5,I4))
    7 FORMAT(2I5)
      END
      SUBROUTINE CPUNCH( AMP,CKAPP,NUM,IDENS,K,AZI,KINDAT)
      CHARACTER*6 ITAG
      CHARACTER*10 IDENS
      DIMENSION IDENS(16)
      DIMENSION AMP(180),CKAPP(180),IAMP(180),KAPP(180)
    1 FORMAT(1H1)
    2 FORMAT(    2I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4, A6)
    3 FORMAT(5X, 2I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4,I5,I4, A6)
    4 FORMAT(8A10/8A10)
      IAZ = IFIX(AZI)
      DO 100 I = 1,NUM
      IAMP(I) = IFIX(AMP(I)*1000.0 + 0.5)
  100 KAPP(I) = IFIX(CKAPP(I)*10.0 + 0.5)
      IF(KINDAT.EQ.2) GO TO 102
      IF(K.EQ.0.AND.IAZ.EQ.0) ITAG = ' Major'
      IF(K.EQ.1.AND.IAZ.EQ.0) ITAG = ' North'
      IF(K.EQ.2.AND.IAZ.EQ.0) ITAG = ' East'
      IF(K.EQ.1.AND.IAZ.GT.0) ITAG = ' Major'
      IF(K.EQ.2.AND.IAZ.GT.0) ITAG = ' Minor'
      IF(K.EQ.1.AND.IAZ.GT.0) ICARD = 1111
      IF(K.EQ.0.AND.IAZ.EQ.0) ICARD = 1111
      IF(K.EQ.1.AND.IAZ.EQ.0) ICARD = 2222
      IF(K.EQ.2.AND.IAZ.GT.0) ICARD = 2222
      IF(K.EQ.2.AND.IAZ.EQ.0) ICARD = 3333
      GO TO 105
  102 ITAG = ' Tides'
      ICARD = 4444
  105 NOX = 0
      ISTOP = NUM/7 + 1
c     PRINT 1
C     PUNCH 4, (IDENS(LL),LL = 1,16)
      DO 200 J = 1,NUM,7
      NOX = NOX + 1
      MM = J/7 + 1
      IEND = J + 6
      IF(IEND.GT.NUM) IEND = NUM
C     PUNCH 2, ICARD, MM, (IAMP(IO),KAPP(IO),IO = J,IEND), ITAG
      PRINT 2, ICARD, MM, (IAMP(IO),KAPP(IO),IO = J,IEND), ITAG
      IF(NOX.EQ.ISTOP) GO TO 888
  200 CONTINUE
  888 RETURN
      END
      SUBROUTINE CSTAT2 (PP,OMEAN,OSD,OVAR,R,SIG,SX,SETS,THBAR,TZERO,
     * ICNTL,ITEST,KAY,M2,M4,       NDEXY,       NV,NCOL)
C
C_______________________________________________________________________
C
C     1 DECEMBER 1983
C
C     LABEL:  CSTAT2T
C
C     PURPOSE:  THIS A MODIFICATION OF CSTAT2P TO TEST THE STABILITY
C               IMPROVEMENT OBTAINED BY WORKING WITH THE CORRELATION
C               MATRIX RATHER THAN THE COVARIANCE MATRIX.
C_______________________________________________________________________
C
C  WRITTEN MAY 1983
C
C  PURPOSE:
C     SUBROUTINE CSTAT2 COMPUTES THE FOLLOWING STATISTICAL VALUES:
C       MEAN
C       STANDARD DEVIATION
C       CORRELATION COEFFICIENTS BETWEEN SPECIFIED VARIABLES
C
C
C
C  THE FOLLOWING ITEMS WILL BE COMPUTED FOR OUTPUT:
C     OVAR   = OBSERVED VARIANCE
C     OMEAN  = OBSERVED MEAN
C     OSD    = STANDARD DEVIATION OF THE RECORD
C     KAY(I) = AN ARRAY FOR INDEXING OUTPUT.
C     R(I)   = CORRELATION COEFFICIENT
C     SIG(I) = STANDARD DEVIATION OF THE CROSS PRODUCT MATRIX
C     SX(L)  = THE SUMS OF ALL VARIABLES, REPLACED BY MEANS
C     PP(L,L) = INITALLY THE SUMS OF CROSS PRODUCTS, REPLACED BY
C               CROSS PRODUCTS WITH THE MEANS REMOVED
C
C________________________________________________________________________
C
C  THIS SUBROUTINE USES THE RELATION
C
C       SUM(X-XBAR)(Y-YBAR) = SUM(X*Y) - SETS * XBAR * YBAR
C
C  WHERE -SETS- IS THE NUMBER OF PAIRED OBSERVATIONS, -XBAR- AND -YBAR-
C  ARE THE MEAN VALUES OF X AND Y.  TO REMOVE THE MEANS FROM ALL CROSS
C  PRODUCTS, MEANS AND STANDARD DEVIATIONS ARE COMPUTED AND MAY BE
C  PRINTED IF DESIRED.  THE CORRELATION MATRIX MAY BE COMPUTED AND
C  PRINTED ALSO.
C
C_______________________________________________________________________
C
C      . . . DIMENSION STATEMENT  .  .  .
C
      REAL    PP(M2,M2)
      REAL SIG(M2),R(M2),SX(M2)
      INTEGER ICNTL(10),KAY(10)
C___________________________________________________________________________
C
c     PRINT 10,(PP(N,N),N=NV+1,NDEXY)
   10 FORMAT(' FROM CSTAT, PP(N,N) =',/,8F14.6,/,8F14.6)
c     PRINT 20, SETS,THBAR,M2,M4,NDEXY
   20 FORMAT (' SETS=',F8.0,'THBAR=',F8.4,' M2=',I4,' M4=',
     * I4,' NDEXY=',I4,/)
C
C******************************************************************************
C
C     . . . STEP 1 . . .
C
C     COMPUTE MEANS
C
C*****************************************************************************
C
      DO 30 I = 1,NDEXY
      KAY(I) = I
   30 SX(I) = SX(I)/SETS
      DO 50 I = 1,NDEXY
      DO 40 J = I,NDEXY
      PP(I,J) = PP(I,J) - SETS * SX(I) * SX(J)
 40   CONTINUE
 50   CONTINUE
C
C     CORRECT CROSS PRODUCTS FOR MEANS, COMPLETE THE CROSS PRODUCT MATRIX,
C     AND COMPUTE THE STANDARD DEVIATION.
C
C     FORM THE CORRELATION MATRIX, IF DESIRED, AND COMPLETE THE LOWER
C     TRIANGULAR MATRIX.
C
C**************************************************************************
C
C     . . . STEP 2 . . .
C
C     FORM STANDARD DEVIATIONS OF ALL VARIABLES
C
C*****************************************************************************
C
   65 DO 60 I = 1,NDEXY
      IF (PP(I,I) .LE. TZERO) THEN
      PRINT 70, I,PP(I,I)
   70 FORMAT (1H0,'IN SUBROUTINE CSTAT2 WITH I =',I4,', PP(I,I) = '
     * , F14.6 / )
      PP(I,I) = TZERO
c     ITEST = 1
      ITEST = 0
      END IF
      SIG(I) =  SQRT(PP(I,I)/SETS)
 60   CONTINUE
C
      OVAR  = PP(NDEXY,NDEXY)/SETS
      OMEAN = SX(NDEXY) + THBAR
      OSD   = SIG(NDEXY)
C
C****************************************************************************
C
C     . . . STEP 3
C
C     COMPUTE CORRELATION MATRIX
C     COMPLETE THE LOWER TRIANGULAR MATRIX
C
C***************************************************************************
C
      DO 100 I = 1,NDEXY
      DO 80 J = I,NDEXY
C
C
      IF (ICNTL(5) .EQ. 1) PP(I,J)=PP(I,J)/(SIG(I)*SIG(J)*SETS)
   80 CONTINUE
C
  100 CONTINUE
C
      DO 120 I = 2,NDEXY
      DO 110 J = 1,I-1
      PP(I,J) = PP(J,I)
  110 CONTINUE
  120 CONTINUE
C
C
      IF (ICNTL(4) .EQ. 0) GO TO 150
C
C**********************************************************************
C
C     . . . STEP 4 . . .
C
C     PRINT MEAN AND STANDARD DEVIATIONS IF THESE ARE WANTED.
C
C**********************************************************************
C
      PRINT 130
      PRINT 140, (KAY(I),SX(I),SIG(I),I=1,NDEXY)
 130  FORMAT(1H1,2X,'VARIABLE NO.',5X,'MEAN VALUE',5X,'STANDARD ',
     *'DEVIATIONS',15X,'VARIABLE NO.',5X,'MEAN VALUE',5X,'STANDARD ',
     *'DEVIATIONS',/)
 140  FORMAT(5X,I4,11X,F8.5,13X,F9.5,24X,I4,11X,F8.5,13X,F9.5)
C
C****************************************************************************
C
C     . . . STEP 5 . . .
C
C     COMPUTE AND PRINT TABLE OF CORRELATION COEFFICIENTS,
C     IF THESE ARE WANTED.
C
C*******************************************************************************
C
C
  150 IF (ICNTL(5).EQ.0 .AND. ICNTL(7) .EQ. 1) THEN
      NCOL = 15
      NP1 = NDEXY
 160  NP2 = NP1 - NCOL+1
      PRINT 170, ( KAY(N),N=NP1,NP2,-1 )
  170 FORMAT (1H1,'CORRELATION COEFFICIENTS BETWEEN INDICATED VARIABLES'
     *,/,(1H ,7X,15(I5,3X)),//)
      IF(NP2 .LT. 1) NP2 = 1
      N2 = NP1 + 1
      DO 200 I = 1,NP1
      N = N2 - I
      JONE = I
      IF(I .GT. NCOL) JONE = NCOL
      IF(JONE .GT. NCOL) JONE=NCOL
      DO 180 J = 1,JONE
      M = N2 - J
      R(J)=PP(N,M)/(SIG(N)*SIG(M)*SETS)
 180  CONTINUE
      PRINT 190, N,(R(L),L=1,JONE)
 190  FORMAT(1H ,I4,3X,15F8.4)
 200  CONTINUE
      IF(NP2 .EQ. 1) GO TO 210
      NP1 = NP1 - NCOL
      GO TO 160
C
      END IF
C
C     PRINT 205, ((PP(I,J) J=41,44) I=J,44)
C 205 FORMAT (1H , ' FROM CSTAT',/,4F12.6,/)
C
  210 RETURN
      END
      SUBROUTINE DATEX(BJOBX,FSTT,TM,JJDAII,DAYB,DAYM,GRBS,GRMS,CP)
 8800 FORMAT(/ 5X, 8H DAYB = ,F5.0,10X, 8H DAYM = ,F5.0 / 5X,  8H GHBS =
     1 ,F7.2,10X, 8H GRMS = ,F7.2 / 5X,17H Original T.M. = ,F7.2,10X,18H
     2 Corrected T.M. = ,F7.2 )
      DAYLEN = 24.00
      DAYB = FLOAT(JJDAII)
      HSER = BJOBX/2.0
      TEST = HSER - AINT(HSER)
      DAYH = AINT(HSER)
  231 CP = TM
  235 GRBS = CP/15.0 + FSTT
      IF(GRBS.LT.24.00) GO TO 240
      GRBS = GRBS - 24.00
      DAYB = DAYB + 1.0
  240 GRMS = CP/15.0 + TEST*DAYLEN
      DAYM = DAYB + DAYH
      IF(GRMS.LT.24.00) GO TO 242
      GRMS = GRMS - 24.00
      DAYM = DAYM + 1.0
  242 PRINT 8800,DAYB,DAYM,GRBS,GRMS,TM,CP
      RETURN
      END
      SUBROUTINE ENT(NPT,ISET,VN,VE,T,SPEEDF,VIN,DIN)
C
C     PROGRAM PREPARED BY E. E. LONG, SEPT., 1979 *** MODIFIED JULY, 1985
C
C ----------------------------------------------------------------------------
C - REVISION: JAN 8 1991
C
C                   REVISED   BY CHRIS ZERVAS        NOVEMBER 1995
C
      REAL*8 SPEED,SPEEDF,VFAC,SPEEDX
      REAL*8 T,TIMX,TXX,STTX,XNSPH,PIRAD
      real*8 jdayb,jbase_date,JULIAN,yearb,monthb,dayb,hourb
C
      CHARACTER*30 IDENST
      CHARACTER*10 LABEL,IDENS,LABEL2
      CHARACTER*10 SKIND,TKIND,LABL
      CHARACTER*40 AQ
      DIMENSION LABEL2(37)
      DIMENSION VN(NPT),VE(NPT),T(NPT)
      DIMENSION VIN(NPT),DIN(NPT)
      DIMENSION SUMN(50),SUME(50),INUM(50)
      DIMENSION SPEEDF(180),ZEPOCH(180)
      DIMENSION NODAYS(12)
c     DIMENSION OEX(5),CXX(30),ZEPOCH(180)
      COMMON/PARAM1/NCONST,NYYYY,NDELTT,CUTOFF,TZERO,MINTRM,MAXT
     1RM,NCOL,NROW,TEBAR
      COMMON/PARAM2/NOBS,MORE,JYY,JMM,JDD,JHH,JMIN,THBAR,XLONG,TIMEZ
!      COMMON/PARAM3/KINDAT,KAZI,AZI
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(200000,2)
     1      ,RATIO                                             
      COMMON/PARAM5/YEPOCH(180),YNODE(180),SPEED(180),P(180)
      COMMON/PARAM6/IDENS(16),LABEL(180),TKIND,SKIND
      COMMON/VEE/TML,CON,U,Q,UI                                               
      COMMON/LOCAT/TM,GONL
      COMMON/MEANL1/MNYEAR,JYEAR,JDAZ,DAY,LDAYY,LLDAYY,MONTH
      COMMON/MEANL2/TIMEX,JDAAYI,JDAYI,JDAYY,LDAA,LXYERE
      COMMON/MEANL3/IALL,LPYER,LYDAY,LJDAA,LNDAY,LGYER,LZDAY
      COMMON/MEANL4/NJOBX,MMONTH,NXDAY,JDAZZ,JDAAZ,FST,INDATA

      NAMELIST /NML1/AZI,JOBX,N,NSPH,CVAR,UMEAN,VMEAN
      NAMELIST /NML2/XYER,MONTH,DAY,STT,STTM,TM,GONL
      NAMELIST /NML3/INDATA,KINDAT,VFAC,AQ,NCON,ITYPE
      NAMELIST /TSPEED/SPEED
      NAMELIST /TLABEL/LABEL

      DATA(NODAYS(J),J = 1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
! zaj 12/21/2006 rearrange the order of tidal consttituents for using less constituents for ITYPE=2 
!	DATA (LABEL2(I),I=1,37)/'M(2)','S(2)','N(2)','K(2)','O(1)',
!     1     'K(1)','P(1)','S(1)','T(2)','R(2)','2N(2)','Q(1)',
!     2     'M(1)','J(1)','S(4)','M(6)','S(6)','M(8)','2Q(1)',
!     3     'L(2)','MK(3)','M(4)','MS(4)','M(3)','2SM(2)','2MK(3)',
!     4     'MN(4)','OO(1)','RHO(1)','MU(2)','NU(2)',
!     5     'LAMBDA(2)','MF','MSF','MM','SSA','SA'/

	DATA (LABEL2(I),I=1,37)/'2Q(1)','O(1)','K(1)','OO(1)','2N(2)',
     1     'M(2)','S(2)','2SM(2)','2MK(3)','MK(3)','M(4)','MS(4)',
     2     'S(4)','M(6)','S(6)','M(8)','Q(1)','M(1)','J(1)','N(2)',
     3     'L(2)','M(3)','MN(4)','RHO(1)','P(1)','MU(2)','NU(2)',
     4     'LAMBDA(2)','K(2)','S(1)','T(2)','R(2)','MF','MSF','MM',
     5     'SSA','SA'/

C
C
    1 FORMAT(1H1)
    8 FORMAT(5F4.0)
   10 FORMAT(32H HARMONIC CONSTANTS DERIVED FROM,I5,43H-DAY  (N + E) HAR
     1MONIC ANALYSES  NOAA - NOS   )
   11 FORMAT(2X,I4,2X,A10,1X,F9.4,4X,F9.3,4X,F9.3,4X,F9.3,F13.5,I8)
   12 FORMAT(// 8X,8HCONSTIT., 7X, 5H( H ), 6X, 7H(KAPPA), 6X, 6H(KP-K),
     1 7X, 8HK(PRIME), 5X, 5HSPEED / )                                       
   13 FORMAT(/// 14X,19HHARMONIC CONSTANTS,5X,13H(H) AND KAPPA  )            
   14 FORMAT(  5X,' HARMONIC CONSTANTS DERIVED FROM A', I5,     
     1 '-DAY HARMONIC ANALYSIS  NOAA - NOS')                                           
   15 FORMAT(/11X,44HELIMINATION OF CONSTITUENT EFFECT CONSIDERED  )         
   23 FORMAT(/11X,48HELIMINATION OF CONSTITUENT EFFECT NOT CONSIDERED  )     
   16 FORMAT(1H1,4X,2A10)                                                    
   17 FORMAT(   5X,8A10/ 5X, 8A10)                                           
   18 FORMAT(/// 6X,11HCONSTITUENT,9X, 3H(R),9X, 4HZETA,7X, 5HSPEED / )      
   24 FORMAT(/// 8X, 8HCONSTIT.,9X, 4HR(P),7X, 7HZETA(P),6X, 5HSPEED / )     
   19 FORMAT(4X,I4,2X,A10,1X,F9.4,4X,F9.3,4X,F9.5,4X,F9.5,F13.5,I8)
   20 FORMAT(8A10/8A10)                                                      
   21 FORMAT( 3I5, 2I10, 2I6, 18X, A10, A5 )                                 
   22 FORMAT( A10, 5I10, 5X, A10, A5 )                                       
   76 FORMAT( F5.0, 7I5,       2F5.0, 2F6.3 )
  221 FORMAT(I5)                                                             
  222 FORMAT(3I5)                                                            
  555 FORMAT( 10A4 )
  808 FORMAT(8A10)                                                           
  810 FORMAT( / 8A10 / 8A10 )                                                
  812 FORMAT( /  10X,35HEqually spaced data beginning Month,I3,4H Day,I3     
     1, 5H Year,I5, 12H  Julian Day,I4 )                                     
  813 FORMAT( /  10X,35H Random spaced data beginning Month,I3,4H Day,I3     
     1, 5H Year,I5, 12H  Julian Day,I4 )                                     
  814 FORMAT( 10X,46HNumber of Julian Days  to  beginning of series,I7)
  815 FORMAT( 10X,46HNumber of Julian Days from beginning of series,I7)
  817 FORMAT( 10X,22HStart time of data set, I3,29H from beginning of se
     1ries is ,F9.2, 7H  hours )
!  816 FORMAT(// 10X,73HPLEASE START YOUR ANALYSIS AT A TIME EQUIVALENT
!     1TO THE LOCAL TIME FOR T.M.,F7.2 /10X,65HTHE COMPUTER WILL ADJUST
!     2THE INPUT TIMES TO THE LOCAL MERIDIAN    )
 2037 FORMAT(5F7.3,F5.0)                                                    
 2038 FORMAT( F5.0, I5, F5.0, F5.2, 2F6.2, 2I5, F11.8 )
 2039 FORMAT(10A4,10A4)
 4040 FORMAT(1H0, 10X, F6.0, 6(1X, F9.2))
 4041 FORMAT(1H0, 4(1X, F9.5))
 7358 FORMAT(1H0,23HJOB TERMINATED........,24H NSPH OR VFAC = 0 (ZERO))     
 7359 FORMAT(6X, F9.1, 5X, F9.1, 4X, 2(1X, F9.1))
 7360 FORMAT(1H0)                                                           
 7361 FORMAT(10X,5H PAGE,I3,7H  OF  3 )                                     
 8128 FORMAT(3I5, 2I10, 2I6, 18X, A10, A5)                                  
 8229 FORMAT( I10, 5I10, 5X, A10, A5 )                                      
 8240 FORMAT(1X,36H *** JOB TERMINATED ** VFAC = 0.0 ** )                   
 8484 FORMAT(5X,66H**** SERIES LESS THAN 29 DAYS ARE NOT ALLOWED ON THIS
     1 PROGRAM ****  )
 8485 FORMAT(/ 10X,29HFirst data point of data set ,I3,' is ',F8.2,F8.0)    
 8486 FORMAT(/// 5X,I6,21H DATA VALUES READ-IN  )                           
 8487 FORMAT(   10X,I6,27H data values to be analyzed )                     
 8488 FORMAT(10X,15HAzimuth used = ,F8.2, 9H  degrees  )
 8489 FORMAT( 9I5, 6A4 )
 8490 FORMAT( /// 10X, 17HMEAN OF SERIES = , F9.5)
 8491 FORMAT( / 10X,52HGreenwich ( V(0)+U ) for the beginning of the ser
     1ies, I5 /  10X,45HNode factors are for the middle of the series )
 8492 FORMAT( /  16X, 10( A10 ))                                            
 8493 FORMAT( 1X, 10H( V(0)+U ),     10(1X,F9.2))
 8494 FORMAT( 4X, 4HNode, 3X,     10(1X,F9.4))
 8495 FORMAT(  1X,10H(AVE.)NODE,     10(1X,F9.4))
 8496 FORMAT(/// 10X,30HTHE AVERAGE NODE FACTOR IS FOR, I3, 6H YEARS  )     
 8497 FORMAT(/// 10X,47HADJUST AVERAGE NODE FACTOR TO MATCH NO. OF YRS.)      
 8498 FORMAT(1H1, 9X,49HLOCAL ( V(0)+U ) FOR THE BEGINNING OF THE SERIES
     1 , I5 )
 8499 FORMAT(// 5X,45H ** SYNODIC CYCLING DEFAULTS FOR CONSTITUENT ,A10)
 8500 FORMAT(// 10X,79HEQUILIBRIUM ARGUMENT ( V(0)+U ) FOR MERIDIAN OF G
     1REENWICH FOR BEGINNING OF YEAR, I5)
 8501 FORMAT(// 10X,45HNODE FACTORS ARE FOR THE MIDDLE OF THE YEAR  )
 8502 FORMAT(// 10X,45HNODE FACTORS ARE FOR THE MIDDLE OF THE SERIES)
 8503 FORMAT(/// 10X,60H**** START (TIME),(DATE),(YEAR) ERROR FOR THIS D
     1ATA SET ****,48H  IT MUST BE GREATER THAN PREVIOUS DATA SET **** )
 8504 FORMAT(/// 10X,34H***ERROR IN START DATE OF DATA SET, I3 )
 8505 FORMAT(/// F8.2, F8.0, 19H    LAST DATA POINT   )
 8506 FORMAT(10X,20HAverage of data set ,2F9.4,17H     Data values ,I6/)
      PIRAD = 3.14159265358979D0/180.0D0
      IF(ISET.NE.0) GO TO 888
c     NJOBX = 0                                             
      NUMX = 0                                              
      IALL = 0                                              
      JUDAA = 0                                             
      LPYER = 0                                             
      LNDAY = 0                                             
      LYDAY = 0
      LZDAY = 0
      LJDAA = 0                                             
      LJUDAA = 0                                            
      NNDAY = 0                                             
      TIMEX = 0.0
      FNODE = 0.0
      DO 886 JK = 1,180                                                       
  886 YNODE(JK) = 0.0
c 887 READ 20,(IDENS(JK),JK=1,16)
c     WRITE(6,810)(IDENS(JK),JK=1,16)
  888 continue
c     READ(5,NML=NML1)
c     READ(5,NML=NML2)
c     READ(5,NML=NML3)
c     WRITE(6,NML=NML1)
c     WRITE(6,NML=NML2)
c     WRITE(6,NML=NML3)
      READ(25,*) AZI,N,NSPH,CVAR,UMEAN,VMEAN
      READ(25,*) XYER,MONTH,DAY,STT,STTM,TM,GONL
      READ(25,*) INDATA,KINDAT,VFAC,NCON,ITYPE
      READ(25,'(a40)') AQ
      print*,' azi = ',azi,' n = ',n,' nsph = ',nsph
      print*,' cvar = ',cvar,' umean = ',umean, ' vmean = ',vmean
      print*,' xyer = ',xyer,' month = ',month,' day = ',day,' stt = ',
     #stt,'sttm = ',sttm
      print*,' tm = ',tm,' gonl = ',gonl
      print*,' indata = ',indata,' kindat = ',kindat,' vfac = ',vfac,
     #' ncon = ',ncon,' itype = ',itype
      THBAR = UMEAN
      TEBAR = VMEAN
      STRT = AINT(STT)
      MODATA = MORE
      JHH = IFIX(STRT)
      JMIN = IFIX(STTM)
      NDELTT = NSPH
      XNSPH = DFLOAT(NSPH)
      XLONG = GONL
      TIMEZ = TM
      STT = STRT + (STTM)/60.0
      TMX = TM
      JYEAR = IFIX(XYER)
      KAZI = 3
      NYYYY = JYEAR
      JMM = MONTH
      JDD = IFIX(DAY)
      NCONST = NCON
      JYY = JYEAR
      ISET = ISET + 1                                                     
      IF(ISET.NE.1) GO TO 7640
 7638 LXYERE = IFIX(XYER)
      LDAYY  = IFIX(DAY)
      MMONTH = MONTH
      FST  = STT
          yearb=LXYERE
          monthb=MMONTH
          dayb=LDAYY
          hourb=STT
      JDAYY=INT(JULIAN(yearb,monthb,dayb,hourb)-
     1       JULIAN(yearb,1.0d0,1.0d0,0.0d0) )+1  !! JUlianday=1 on Jan. 1

      IDENS(9)  = 'Least Squa'
      IDENS(10) = 'res H.A.  '
      IDENS(11) = 'Beginning '
      write(IDENST,3011) mmonth,ldayy,lxyere,fst
 3011 FORMAT(I2,1H-,I2,1H-,I4,10H  at Hour ,F5.2,2H     )
      idens(12)(1:10)=idenst(1:10)
      idens(13)(1:10)=idenst(11:20)
      idens(14)(1:10)=idenst(21:30)
 7640 continue
c7640 NJOBX = NJOBX + JOBX
      IF(MODATA.NE.0) GO TO 7645
      IF(MODATA.EQ.0.AND.NJOBX.GE.29) GO TO 7645
      PRINT 8484   ! zaj 12/12/2006 commentted out to perform HA for shorter time series with less constituents
      CALL PEXIT
 7645 LDAA = IFIX(DAY)
      LYERE = IFIX(XYER)                                                  
          yearb=LYERE
          monthb=MONTH
          dayb=LDAA
          hourb=0.0
      JDAYI=INT(JULIAN(yearb,monthb,dayb,hourb)-
     1       JULIAN(yearb,1.0d0,1.0d0,0.0d0) )+1  !! JUlianday=1 on Jan. 1

      JDAZ = JDAYI                                                        
      JDAYI = JDAYI - 1
 7648 IF(IALL.EQ.0) GO TO 7660
      IF(LYERE.EQ.LPYER) GO TO 7660
 7649 LNDAY = LYDAY + LZDAY
 7650 LZDAY = LNDAY
      LGYER = LYERE - LPYER
      IF(LGYER.LE.1) GO TO 7660
      LPYER = LPYER + 1
      CALL YTEST(LPYER,NODAYS,MMDAY)
      LNDAY = LNDAY + MMDAY
      GO TO 7650
 7660 JDAYI = JDAYI + LNDAY
      JUDAA = JDAYI + 1
      IF(JDAYI.LT.LJDAA) THEN
      PRINT 8504, ISET
      CALL PEXIT
      ENDIF
      IF(IALL.GT.0) GO TO 7676                                            
      GO TO (7672,7674,7675),ITYPE                                        
c7672 READ(5,NML=TSPEED)
 7672 READ(5,*)(SPEED(i),i=1,NCON)
      GO TO 7675                                                          
c7674 READ(5,NML=TLABEL)
C 7674 READ(5,*)(LABEL(i),i=1,NCON)
7674  CONTINUE
      DO I=1,NCON
        LABEL(i)=LABEL2(i)
      ENDDO
 7675 DO 2222 III = 1,NCON
      GO TO (195,197,198),ITYPE
  195 SPEEDX = SPEED(III)
      CALL NNAME(SPEEDX,LABL,ISUB,INUMX,1)
      LABEL(III) = LABL
      GO TO 2222
  197 LABL = LABEL(III)
      CALL NNAME(SPEEDX,LABL,ISUB,INUMX,2)
      SPEED(III) = SPEEDX
      GO TO 2222
  198 INUMX = III
      CALL POSITX(SPEEDX,IPOS,INUMX,2)
      SPEED(III) = SPEEDX
      CALL NNAME (SPEEDX,LABL,ISUB,INUMX,1)
      LABEL(III) = LABL
 2222 P(III) = FLOAT(ISUB)
      DO 2225 II = 1,NCON
 2225 SPEEDF(II) = (SPEED(II)*PIRAD)
 7676 CONTINUE                                                             
      IF(IALL.GT.0) GO TO 7677                                             
      TIMX = FLOAT(JDAYI)*24.00d0 + dble(STT)  
      TIMEX = FLOAT(JDAYI)*24.00 + STT  
      JDAAYI =  JDAYI
      MNYEAR = LYERE                                                       
      JDAAZ = JDAZ
      LLDAYY = LDAYY
c     PRINT 1                                                              
      PRINT 810,(IDENS(KJ),KJ = 1,16)                                      
 7677 IF(NSPH.EQ.0.OR.VFAC.EQ.0.0D0) GO TO 8356
      IF(KINDAT.EQ.1) GO TO 102
  101 IF(INDATA.EQ.1)THEN 
	do kz=1,N                                
          VIN(kz)=RDATA(kz,1)
	ENDDO
      ENDIF
!      IF(INDATA.EQ.1) READ (9,AQ)(VIN(I),I=1,N)
!      IF(INDATA.EQ.2) READ (9,AQ)(T(I),VIN(I),I=1,N)
      CALL TEMP(VIN,VN,N)                                                  
      GO TO 108                                                            
  102 continue 
	do kz=1,N                                
          VIN(kz)=RDATA(kz,1)
          DIN(kz)=RDATA(kz,2)
	ENDDO
!      IF(INDATA.EQ.1) READ(9,AQ) (VIN(I),DIN(I),I = 1,N)
!      IF(INDATA.EQ.2) READ(9,AQ) (T(I),VIN(I),DIN(I),I = 1,N)
      PRINT 8485, ISET,VIN(1),DIN(1)                                       
 8355 CALL AZIM(DIN,VIN,AZI,VN,VE,KAZI,N,CVAR)
      PRINT 8488, AZI                                                      
  108 PRINT 8487, N
      IF(VFAC.NE.0.0D0) GO TO 1083                                         
      PRINT 8240                                                           
      CALL PEXIT                                                           
 1083 IF(VFAC.GE.0.9999D0.AND.VFAC.LE.1.0D0) GO TO 7680
      DO 1084 IZP = 1,N                                                    
      VE(IZP) = VE(IZP)/VFAC
 1084 VN(IZP) = VN(IZP)/VFAC
      GO TO 7680
 8356 PRINT 7358
      CALL PEXIT
 7680 IF(MODATA.NE.0) GO TO 358
      ZYERE = FLOAT(MNYEAR)
      CALL YTEST(MNYEAR,NODAYS,MZDAY)
      BJOBX = FLOAT(NJOBX)
  995 CALL DATEX(BJOBX,FST,TM,JDAAZ,DAYB,DAYM,GRBS,GRMS,CP)
      TML = CP - GONL
      NJOBX = IFIX(BJOBX)
      CALL ASTRO(ZYERE,DAYB,DAYM,GRBS,GRMS,NJOBX)
      DO 200 III = 1,NCON
      SPEEDX = SPEED(III)
      CALL VANDUF(SPEEDX,EPOXX,FNODE,2)
  200 YNODE(III) = FNODE
      DJOBX = FLOAT(MZDAY)
      CALL DATEX(DJOBX,0.0,TM,1,DAYBB,DAYMM,GRBSS,GRMSS,CP)
      TML = CP - GONL
      CALL ASTRO(ZYERE,DAYBB,DAYMM,GRBSS,GRMSS,MZDAY)
      DO 2200 III = 1,NCON
      SPEEDX = SPEED(III)
      CALL VANDUF(SPEEDX,EPOXX,FNODE,2)
      ZEPOCH(III) = DBLE(EPOXX) + DBLE(TIMEX)*SPEED(III)
C     ZEPOCH(III) IS THE  GR(V(0)+U) FOR BEGINNING OF SERIES
 2200 YEPOCH(III) = ZEPOCH(III)
C     EPOXX = GR(V(0)+U) FOR BEGINNING OF YEAR
      CALL TWOPI(ZEPOCH,NCON)
      CALL TWOPI(YEPOCH,NCON)
      PRINT 8491, LXYERE
      DO 997 IV = 1,NCON,10                                                 
      NIVE = IV + 9                                                         
      IF(NIVE.GT.NCON) NIVE = NCON                                          
      PRINT 8492, ( LABEL(LL),LL = IV,NIVE)                                 
      PRINT 8493, (YEPOCH(LL),LL = IV,NIVE)                                 
      PRINT 8494, ( YNODE(LL),LL = IV,NIVE)                                 
      IF(NIVE.EQ.NCON) GO TO 358                                            
  997 CONTINUE                                                              
c 998 PRINT 8502
  358 NUM = N                                                               
      NOBS = N
      INUM(ISET) = NUM                                                      
      IDAY = IFIX(DAY)                                                      
      JDAZZ = JUDAA - 1                                                     
      JDAYII = JDAZZ - JDAAYI
      STTX = FLOAT(JDAZZ)*24.00d0 + dble(STT) - TIMX
      TIMX = STTX
      JDAZZX = JDAZZ
      STTTX = STTX
      GO TO (3587,3588),INDATA                                              
 3587 PRINT 812, MONTH, IDAY, JYEAR, JDAZ                                   
      PRINT 814, JDAZZ
      PRINT 815, JDAYII
      PRINT 817, ISET, STTX
      IF(STTX.GE.0.0D0) GO TO 359
 6008 PRINT 8503
      CALL PEXIT
 3588 PRINT 813, MONTH, IDAY, JYEAR, JDAZ                                   
      PRINT 814, JDAZZ
      PRINT 815, JDAYII
      PRINT 817, ISET, STTX
      IF(STTX.LT.0.0D0) GO TO 6008
 3592 TXX = STTX - (T(1) - 1.)*24.
      DO 3600 KK = 1,NUM
      T(KK) = (T(KK) - 1.)*24. + TXX
 3600 CONTINUE
  359 GO TO (360,403),INDATA                                              
  360 DO 400 III = 1,NUM                                                  
      T(III) = TIMX + DFLOAT(III-1)/XNSPH
  400 CONTINUE                                                            
  403 SUMN(ISET) = 0.0
      SUME(ISET) = 0.0
      DO 4444 J = 1,NUM
      SUMN(ISET) = SUMN(ISET) + VN(J)
      SUME(ISET) = SUME(ISET) + VE(J)
 4444 CONTINUE
      AVSUMN = SUMN(ISET)/FLOAT(NUM)
      AVSUME = SUME(ISET)/FLOAT(NUM)
 4460 PRINT 8506, AVSUMN,AVSUME,NUM
      LPYER = LYERE
      LYDAY = NNDAY
      LJDAA = JUDAA
      IALL = IALL + 1
      RETURN
      END
      SUBROUTINE FITAN( AUS,AUC,RTA,SPED,JMAP)                            
      SPED = 0.0                                                          
      IF(AUC.EQ.0.0) GO TO 14
      BXG = AUS/AUC
      GO TO 15
   14 IF(AUS)16,44,17                                                     
   16 RTA = 270.0                                                         
      GO TO 88                                                            
   17 RTA = 90.0                                                          
      GO  TO 88                                                           
   15 RTA = ATAN(BXG)*57.29578
      GO TO (88,38),JMAP                                                  
   38 IF(AUS)32,32,34                                                    
   32 IF(AUC)33,44,37                                                    
   33 RTA = 180.0 + RTA                                                  
      GO TO 88                                                           
   34 IF(AUC)35,17,88                                                    
   35 RTA = RTA + 180.0                                                  
      GO TO 88                                                           
   37 RTA = RTA + 360.0                                                  
      GO TO 88                                                           
   44 RTA = 0.0                                                          
   88 IF(AUC.EQ.0.0.AND.AUS.EQ.0.0) GO TO 55                             
      SPED = SQRT(AUC**2 + AUS**2)                                       
   55 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      SUBROUTINE HEADER(NECOM,K,AZI,KINDAT,TKIND,SKIND)                  
      CHARACTER*10 SKIND,TKIND
      IF(KINDAT.EQ.2) GO TO 7760                                         
      IF(AZI.EQ.0.0) GO TO 7750                                          
      IF(K.EQ.1) TKIND = ' (Major Ax'                                    
      IF(K.EQ.2) TKIND = ' (Minor Ax'                                    
      SKIND = 'is)'                                                      
      GO TO 7760                                                         
 7750 IF(K.EQ.1) TKIND = ' (North Co'                                    
      IF(K.EQ.2) TKIND = '  (East Co'                                    
      SKIND = 'mponent)'                                                 
 7760 IF(KINDAT.EQ.2) TKIND = ' ( Tides )'
      IF(KINDAT.EQ.2) SKIND = '      '
      RETURN                                                             
      END                                                                
      SUBROUTINE MATCOR(PP,CUTOFF,JR,M2,NCOL,NDEXY,NDEXP1,NROW,NSTEP,
     * MATYP)
C_______________________________________________________________________
C
C     FILE NAME:     MATCORTT
C
C     WRITTEN:     JANUARY 1984
C
C     PURPOSE:
C
C     SUBROUTINE MATCOR IS CALLED BY SUBROUTINE SCREEN TO DISPLAY THE
C     ORIGINAL CORRELATION MATRIX, AND ALSO BY SCREEN TO DISPLAY
C     THE MATRIX AS MODIFIED BY THE SOLUTION ALGORITHM.  CORRELATION
C     COEFFICIENTS ARE MULTIPLIED BY 1,000,000 TO ELIMINATE UNNEEDED
C     ZEROS AND DECIMALS.  OUTPUT CAN BE LIMITED TO NXN COLUMNS AND
C     ROWS, WHERE N IS SPECIFIED BY NCOL.  IT IS ANTICIPATED THAT THIS
C     SUBROUTINE WILL BE USED ONLY FOR TRAINING AND UNSTANDARD ANALYSES.
C
C_______________________________________________________________________
C
      REAL    PP(M2,M2)
      DIMENSION JR(1)
C
C_________________________________________________________________________
C
C*********************************************************************
C
C              STEP 1
C
C***********************************************************************
C
      NDEXP2 = NDEXP1 + 1
      IF (NCOL .GT. NDEXY) NCOL = NDEXY
      NPAGE = NCOL/15
      IF (NPAGE*15 .LT. NCOL) NPAGE=NPAGE+1
      J2 = 0
      IF (NPAGE .EQ. 0 .OR. NCOL .LE. 14) GO TO 170
C
C***********************************************************************
C
C              STEP 2 - A
C
C*********************************************************************
C
      DO 160 I=1,NPAGE
      J1 = J2 + 1
      J2 = J1 + 14
      IF (J2 .GE. NCOL) J2 = NCOL
C
C     PRINT PAGE HEADER
C
      PRINT 110, NSTEP,MATYP,I
  110 FORMAT (1H1, ' CORRELATION MATRIX X 1,000,000, STEP ',I4,
     * ' MATRIX TYPE',I4,' PAGE NO.',I4,/)
C
C      LABEL COLUMNS
C
      PRINT 115, (J,J=J1,J2)
  115 FORMAT (1H0,'LINE',123X,'VAR.',/,' NO.',1X,15(I5,3X),
     * 3X,'NO.',/)
C
C*********************************************************************
C
C              STEP 3 - A
C
C*********************************************************************
C
      DO 140 K = 1, NROW-1
C
C     SET UP ONE LINE
      DO 120 J = J1,J2
C
      IF (ABS(PP(K,J)) .LT. CUTOFF) THEN
      JR(J) = 0
      ELSE
      JR(J) = INT(PP(K,J)*1000000)
      END IF
C
      IF (J .NE. NCOL) GO TO 120
C
      IF (ABS(PP(K,NDEXY)) .LT. CUTOFF) THEN
      JR(J2) = 0
      ELSE
      JR(J2) = INT(PP(K,NDEXY)*1000000)
      END IF
C
  120 CONTINUE
C
      JR(J2+1) = INT(PP(K,NDEXP1))
C
      PRINT 130, K,(JR(J),J=J1,J2+1)
  130 FORMAT (1X,I3,15I8,I4)
C
  140 CONTINUE
C
      DO 150 J=J1,J2
C
      IF (ABS(PP(NDEXY,J)) .LT. CUTOFF) THEN
      JR(J) = 0
      ELSE
      JR(J) = INT(PP(NDEXY,J)*1000000)
      END IF
C
      IF (J2 .NE. NCOL) GO TO 150
C
      IF (ABS(PP(NDEXY,NDEXY)) .LT. CUTOFF) THEN
      JR(J2) = 0
      ELSE
      JR(J2) = INT(PP(NDEXY,NDEXY)*1000000)
      END IF
C
  150 CONTINUE
C
C     SET UP LAST LINE
C
      JR(J2+1) = INT(PP(NDEXY,NDEXP1))
      PRINT 130, NDEXY,(JR(J), J=J1,J2+1)
C
  160 CONTINUE
C
C*********************************************************************
C
C              STEP 2 - B
C
C*********************************************************************
C
C     SET UP FINAL PAGE
  170 J1 = J2 + 1
      J2 = NCOL-1
      PRINT 110, NSTEP,MATYP,NPAGE+1
      PRINT 185, (J,J=J1,J2),NDEXY
  185 FORMAT(1H0, 'LINE',123X,'VAR',/,1X,'NO.',1X,15(I5,3X),3X,'NO.',/)
      DO 195 K=1,NROW-1
C
C*********************************************************************
C
C              STEP 3 - B
C
C*********************************************************************
C
C     SET UP ONE LINE
      DO 205 J=J1,J2
      IF (ABS(PP(K,J)) .LT. CUTOFF) THEN
      JR(J) = 0
      ELSE
      JR(J) = INT(PP(K,J)*1000000)
      END IF
  205 CONTINUE
C
      IF (ABS(PP(K,NDEXY)) .LT. CUTOFF) THEN
      JR(J2+1) = 0
      ELSE
      JR(J2+1) = INT(PP(K,NDEXY)*1000000)
      END IF
C
      PRINT 220, K,(JR(J),J=J1,J2+1),INT(PP(K,NDEXP1)),INT(PP(K,NDEXP2))
  220 FORMAT (1X,I3,15I8,2I3)
C
  195 CONTINUE
C
      DO 230 J=J1,J2
      IF (ABS(PP(NDEXY,J)) .LT. CUTOFF) THEN
      JR(J) = 0
      ELSE
      JR(J) = INT(PP(NDEXY,J)*1000000)
      END IF
C
  230 CONTINUE
      JR(J2+1) = INT(PP(NDEXY,NDEXY)*1000000)
      PRINT 240, NDEXY,(JR(J),J=J1,J2+1)
  240 FORMAT (1X,I3,15I8)
      RETURN
      END
      SUBROUTINE MODT1(PP,FREQ,FHOUR,SX,TH,THBAR,V,I1,
     * M2,NDEXY,NOBS,NV,TIM,NDR)
C_______________________________________________________________________
C
C     1 DEC 83
C     THIS SUBROUTINE WAS TRANSFERED FROM JACK'S FILE UNDER THE NAME
C     -MODT1.0-.  THE NAME HAS BEEN CHANGED TO -MODT1P- FOR PERMANENT.
C
C_____________________________________________________________________________
C
C    WRITTEN MARCH 1983
C  PURPOSE:
C     SUBROUTINE MODT1 (MODIFY DATA, OPTION 1) IS DESIGNED FOR THE
C  ROUTINE ANALYSIS FOR ASTRONOMICAL TIDES OF HOURLY HEIGHTS, OR
C  REGULAR SPACING OF OBSERVATIONS OTHER THAN HOURLY.
C     REGRESSION VECTORS ARE CONSTRUCTED IN THE FORM
C        X(1), X(2), X(3), . . ., X(NV), Y
C  WHERE THE X(I) ARE TRIGONOMETRIC FUNCTIONS AND Y IS THE OBSERVED
C  TIDE HEIGHT.  AFTER THE CONSTRUCTION OF EACH VALUE OF THE VECTOR
C  SUMS AND SUMS-OF-PRODUCTS ARE ACCUMULATED. ONE OR MORE CALLS CAN
C  BE MADE TO SUBROUTINE MODT1 IN THE ANALYSIS OF ANY TIDE RECORD.
C  THIS PROVISION PERMITS BYPASSING DATA GAPS AND LARGE STORM TIDES.
C
C_______________________________________________________________________
C
C               DEFINITIONS
C
C     JTEST NOT EQUAL TO ZERO.  INDICATES A NEGATIVE PHASE ANGLE.
C            THIS IS IMPOSSIBLE IF CODE IS CORRECT.
C  +  NOBS    = NUMBER OF OBSERVATIONS IN THIS SET.  COMPUTED BY
C              SUBROUTINE REDATA.
C     SETS    = TOTAL NUMBER OF OBSERVATIONS TO BE CONSIDERED.
C                 IF THERE IS ONE SET OF INPUT DATA WITHOUT BREAKS
C                 SETS WILL BE EQUAL TO NOBS.
C     SX(I)   = SUM OF X(I) OVER ALL -I- AFTER ALL DATA ARE INCLUDED.
C               SX(I) IS REPLACED BY THE MEAN VALUE OF X(I).
C
C     JYY     = YEAR
C     JDAY(I) = NUMBER OF JULIAN DAYS -1 COMPLETED.  ?BLOCK DATA?
C     JMONTH  = INDEX FOR JDAY TABLE.
C  *  JMM     = MONTH OF THE YEAR.
C  *  JDD     = DAY OF THE MONTH.
C  *  JHH     = HOUR OF THE DAY (FROM 0 TO 23).
C     JMIN    = MINUTE (IN TENTHS OF HOURS)
C     HOUR    = (JULIAN HOUR - 1) OF THE FIRST OBSERVATION.
C **  NCONST  = NUMBER OF CONSTITUENTS CONSIDERED.
C     THETA(J) = THE PHASE ANGLE IN COMPUTATIONAL UNITS.
C  *  TH(I)   = OBSERVED TIDE HEIGHT.
C     THBAR   = ESTIMATED MEAN TIDE LEVEL.
C     V(J)    = REGRESSION VECTOR.
C     FREQ(N) = FREQUENCY IN COMPUTATIONAL UNITS, THE CHANGE IN
C                ONE HOUR.
C     JTHETA  = INTEGER PORTION OF THETA - LOCATION OF THE DESIRED SINE
C     TWOPI   = 2 TIMES PI.   ?BLOCK DATA?
C
C   * INPUT DATA SUPPLIED BY SUBROUTINE REDATA.
C  ** INPUT DATA SUPPLIED BY LSQHA.
C   + COMPUTED BY SUBROUTINE REDATA2.
C  ++ COMPUTED BY LSQHA.
C
C_______________________________________________________________________
C
C
      REAL*8 FREQ,FHOUR,TIM,DTHETA,HOUR,DELTT,SHOUR,DELTAX
      REAL   PP(M2,M2)
      DIMENSION TH(1), FREQ(1),V(1),SX(1),FHOUR(1),TIM(1)
      COMMON/BGROUP/DELTT,SHOUR(20),DELTAX
C
C___________________________________________________________________________
C
C  INITIALIZE TEST FOR NEGATIVE ANGLES.  THIS TEST TO BE REMOVED
C  AT THE END OF CHECKOUT.
C     JTEST = 0
C
C*************************************************************************
C
C     . . . STEP 1 . . .
C
C INITIALIZE JULIAN HOUR COMPUTATION.
C
C*****************************************************************************
C
C     FOR HOURLY DATA, HOUR IS CALCULATED AS (JULIAN HOUR-1)
C     TO INITIALIZE FOR ROUTINE INCREMENTATION BY THE LOOP DO 100.
C     HOUR = ((JDAY(JMONTH)+JDD) * 24+JHH+JMIN/60.-1)/DELTT
      HOUR = DINT(TIM(1)*DELTAX + 0.001D0)
      SHOUR(NDR) = HOUR
      FHOUR(NDR) = DINT(TIM(NOBS)*DELTAX + 0.001D0)
C  ONE TIME PERIOD IS SUBTRACTED TO INITIALIZE THE PHASE ANGLE
C  ESTABLISH REGRESSION VECTOR.
c     PRINT 5, I1,NOBS
 5    FORMAT(3X,'I1=',I4,'  NOBS=',I6,/)
C
C***************************************************************************
C
C     . . . STEP 2 . . .
C
C***************************************************************************
C
      DO 125 I=I1,NOBS
C
      V(NDEXY) = TH(I) - THBAR
C  ESTABLISH PERIODIC INDEPENDENT OBSERVATIONS IN REGRESSION VECTOR.
      DO 100 J = 1,NV,2
C  UPDATE PHASE ANGLE FOR EACH CONSTITUENT.
      N = (J + 1)/2
      DTHETA = FREQ(N)*TIM(I)
      CALL TWOPIR(DTHETA,1)
C  ESTABLISH INDEPENDENT VARIABLES IN THE REGRESSION VECTOR,
C  V(2J) = V((J + 1)/2) = SINE OF PHASE ANGLE, V(2J) = COSINE OF
C  PHASE ANGLE.
C     V(J) = SIN(THETA)
C     V(J+1) = COS(THETA)
C
C     THE SYSTEM FUNCTIONS WILL BE USED IN PLACE OF THE SINE TABLE.
C
      V(J) = DSIN(DTHETA)
      V(J+1) = DCOS(DTHETA)
C
  100 CONTINUE
C
C******************************************************************************
C
C     . . . STEP 3 . . .
C
C******************************************************************************
C
C
C*****************************************************************************
C
C     . . . STEP 4 . . .
C
C     FORM THE SUMS OF ALL NON-PERIODIC VARIABLES.
C
C******************************************************************************
C
  110 DO 162 J = NV+1,NDEXY
      SX(J) = SX(J) + V(J)
C     FORM THE SUM OF CROSS PRODUCTS OF DEPENDENT AND
C     INDEPENDENT VARIABLES.
      DO 161 K=1,J
  161 PP(K,J) = PP(K,J) + V(K) * V(J)
  162 CONTINUE
C
C  THE FOLLOWING IF - PRINT STATEMENT IS USED ONLY FOR TESTING.
C     IF ( I .LE. I1+3 ) PRINT 165, (V(J),J=1,NDEXY)
C 165 FORMAT (1H0,'TEST V(J):',/,(1X,10F12.5,/))
C
  125 CONTINUE
C
C
      RETURN
      END
      SUBROUTINE NNAME(SPDD,ITAG,ISUB,INUM,ICODE)
C     THIS SUBROUTINE IDENTIFIES THE CONSTITUENT BY ITS SPEED               1618
C     AND MAKES IT AVAILABLE FOR LABELING                                   1619
C     OR IDENTIFIES THE CONSTITUENT BY LABEL AND MAKES AVAILABLE ITS        1620
C     CONSTITUENT SPEED                                                     1621
C     IT ALSO DETERMINE THE SUBSCRIPT OF THE CONSTITUENT                    1622
C     SUBROUTINE PREPARED BY -- E. E. LONG  (NOAA/NOS) 1979
C        ORDER OF CONSTITUENT SPEEDS***  M(2),N(2),S(2),O(1),K(1),K(2)      1623
C      L(2),2N(2)R(2),T(2),LAMBDA(2),MU(2),NU(2),J(1),M(1),OO(1),P(1)       1624
C      Q(1),2Q(1),RHO(1),M(4),M(6),M(8),S(4),S(6),M(3),S(1),MK(3),2MK(3)    1625
C      MN(4),MS(4),2SM(2),MF,MSF,MM,SA,SSA                                  1626
      REAL*8 SPD,SPDD
      CHARACTER*10 LABLE(180),ITAG
      DIMENSION SPD(180),IP(180)
      DATA(IP(N),N = 1,37)/3*2,2*1,8*2,7*1,4,6,8,4,6,3,1,2*3,2*4,2,5*0/
      DATA IP(38)/3/
      DATA(IP(N),N=39,114)/5*1,13*2,2*3,7*4,8*5,12*6,4*7,10*8,4*9,6*10,
     111,4*12/
      DATA(IP(N),N=115,120)/3,4,2*6,2*8/
      DATA(IP(N),N=121,128)/4*1,2*2,3,4/                                  
      DATA(IP(N),N=129,140)/3*2,4,3*6,8,10,12,2*1/                        
      DATA(IP(N),N=141,150)/2*1,4*2,4*4/
      DATA(IP(N),N=151,164)/2,11*0,3,5/
      DATA(IP(N),N=165,175)/2*2,3*3,3*4,7,8,10/
      DATA(SPD(I),I=  1, 37)/ 28.9841042D0, 28.4397295D0, 30.0000000D0,   
     1    13.9430356D0,   15.0410686D0,   30.0821373D0,   29.5284789D0,   
     2    27.8953548D0,   30.0410667D0,   29.9589333D0,   29.4556253D0,   
     3    27.9682084D0,   28.5125831D0,   15.5854433D0,   14.4966939D0,   
     4    16.1391017D0,   14.9589314D0,   13.3986609D0,   12.8542862D0,   
     5    13.4715145D0,   57.9682084D0,   86.9523127D0,  115.9364169D0,   
     6    60.0000000D0,   90.0000000D0,   43.4761563D0,   15.0000000D0,   
     7    44.0251729D0,   42.9271398D0,   57.4238337D0,   58.9841042D0,   
     8    31.0158958D0,   01.0980331D0,   01.0158958D0,   00.5443747D0,   
     9    00.0410686D0,   00.0821373D0/                                   
      DATA SPD(38)/ 45.0410686D0/                                         
      DATA(SPD(I),I= 39, 77)/ 12.9271398D0, 14.0251729D0, 14.5695476D0,   
     A    15.9748272D0,   16.0569644D0,   26.8794590D0,   26.9615963D0,   
     B    27.4238337D0,   27.5059710D0,   27.8039338D0,   28.6040041D0,   
     C    28.9019669D0,   29.0662415D0,   29.1483788D0,   29.3734880D0,   
     D    30.5443747D0,   30.7086493D0,   31.0980331D0,   42.3827651D0,   
     E    43.9430356D0,   56.8794590D0,   56.9523127D0,   57.5059710D0,   
     F    58.4397295D0,   58.5218668D0,   59.0662415D0,   59.5284789D0,   
     G    71.3668693D0,   71.9112440D0,   71.9933813D0,   72.4649023D0,   
     H    72.9271398D0,   73.0092770D0,   74.0251728D0,   74.1073100D0,   
     I    85.4013258D0,   85.8635632D0,   85.9457005D0,   86.4079380D0/   
      DATA(SPD(I),I= 78,116)/ 86.4900752D0, 87.4238337D0, 87.5059710D0,   
     J    87.9682084D0,   88.0503457D0,   88.5218668D0,   88.9841042D0,   
     K    89.0662415D0,  100.3509735D0,  100.9046318D0,  101.9112440D0,     
     L   103.0092771D0,  114.8476674D0,  115.3920422D0,  115.4741794D0,     
     M   116.4079380D0,  116.4900752D0,  116.9523127D0,  117.0344500D0,     
     N   117.5059710D0,  117.9682084D0,  118.0503457D0,  129.8887360D0,     
     O   130.4331108D0,  130.9774855D0,  131.9933813D0,  144.3761464D0,     
     P   144.9205211D0,  145.3920422D0,  145.9364169D0,  146.4900752D0,     
     Q   146.9523127D0,  160.9774855D0,  174.3761464D0,  174.9205211D0,     
     R   175.4741794D0,  175.9364169D0,   41.9205276D0,   58.5125831D0/     
      DATA(SPD(I),I=117,140)/ 87.4966873D0, 88.5125831D0,117.4966873D0,     
     S   116.4807916D0,   14.9178647D0,   15.0821353D0,   15.1232059D0,     
     T    15.5125897D0,   30.6265119D0,   27.3416965D0,   42.9271397D0,     
     U    60.0821373D0,   26.9523126D0,   27.4966873D0,   28.5947204D0,     
     V    57.4966873D0,   85.3920421D0,   85.9364168D0,   86.4807916D0,     
     W   115.4648958D0,  146.4807916D0,  175.4648958D0,   16.1391016D0,     
     X    12.8450026D0/                                                     
      DATA(SPD(I),I=141,150)/ 15.1232058D0, 14.8767942D0, 30.0000001D0,
     Y    29.9178627D0,   30.1642746D0,   29.9178666D0,   59.9589333D0,
     Z    59.9178627D0,   60.2464119D0,   59.8767999D0/
      DATA(SPD(I),I=151,164)/ 28.9430356D0, 01.0569644D0, 00.5490165D0,
     1          00.5079479D0, 00.0410667D0, 00.1232059D0, 00.1642746D0,
     2          00.2464118D0, 00.3285841D0, 00.4106864D0, 00.4928237D0,
     3          00.9856473D0, 45.0000000D0, 75.0000000D0/
      DATA(SPD(I),I=165,175)/ 27.8860712D0, 30.0410686D0, 43.4807981D0,
     4          44.9589314D0, 45.1232059D0, 56.3258007D0, 56.8701754D0,
     5          57.8860712D0,105.0000000D0,120.0000000D0,150.0000000D0/
      DATA(LABLE(M),M = 1,37)/'M(2)'     ,'N(2)'     ,'S(2)'     ,
     1          'O(1)'     ,'K(1)'     ,'K(2)'     ,'L(2)'     ,
     2          '2N(2)'    ,'R(2)'     ,'T(2)'     ,'LAMBDA(2)',
     3          'MU(2)'    ,'NU(2)'    ,'J(1)'     ,'M(1)'     ,
     4          'OO(1)'    ,'P(1)'     ,'Q(1)'     ,'2Q(1)'    ,
     5          'RHO(1)'   ,'M(4)'     ,'M(6)'     ,'M(8)'     ,
     6          'S(4)'      ,'S(6)'      ,'M(3)'     ,'S(1)'     ,
     7          'MK(3)'    ,'2MK(3)'   ,'MN(4)'    ,'MS(4)'    ,
     8          '2SM(2)'    ,'MF'        ,'MSF'      ,'MM'        ,
     9          'SA'        ,'SSA'       /
      DATA LABLE(38)/'SK(3)'     /
      DATA(LABLE(M),M =39,77)/'SIGMA(1)'  ,'MP(1)'     ,'CHI(1)'    ,
     A          '2PO(1)'    ,'SO(1)'     ,'2NS(2)'    ,'2NK2S(2)'  ,
     B          'MNS(2)'    ,'MNK2S(2)'  ,'2MS2K(2)'  ,'2KN2S(2)'  ,
     C          'OP(2)'     ,'MKS(2)'    ,'M2(KS)(2)' ,'2SN(MK)(2)',
     D          'MSN(2)'    ,'2KM(SN)(2)','SKM(2)'    ,'NO(3)'     ,
     E          'SO(3)'     ,'N(4)'      ,'3MS(4)'    ,'MNKS(4)'   ,
     F          'SN(4)'     ,'KN(4)'     ,'MK(4)'     ,'SL(4)'     ,
     G          'MNO(5)'    ,'2MO(5)'    ,'3MP(5)'    ,'MNK(5)'    ,
     H          '2MP(5)'    ,'2MK(5)'    ,'MSK(5)'    ,'3KM(5)'    ,
     I          '3NKS(6)'   ,'2NM(6)'    ,'2NMKS(6)'  ,'2MN(6)'    /
      DATA(LABLE(M),M=78,114)/'2MNKS(6)'  ,'MSN(6)'    ,'MKN(6)'    ,
     J          '2MS(6)'    ,'2MK(6)'    ,'NSK(6)'    ,'2SM(6)'    ,
     K          'MSK(6)'    ,'2MNO(7)'   ,'2NMK(7)'   ,'2MSO(7)'   ,
     L          'MSKO(7)'   ,'2(MN)(8)'  ,'3MN(8)'    ,'3MNKS(8)'  ,
     M          '2MSN(8)'   ,'2MNK(8)'   ,'3MS(8)'    ,'3MK(8)'    ,
     N          'MSNK(8)'   ,'2(MS)(8)'  ,'2MSK(8)'   ,'2M2NK(9)'  ,
     O          '3MNK(9)'   ,'4MK(9)'    ,'3MSK(9)'   ,'4MN(10)'   ,
     P          'M(10)'     ,'3MNS(10)'  ,'4MS(10)'   ,'2MNSK(10)' ,
     Q          '3M2S(10)'  ,'4MSK(11)'  ,'4MNS(12)'  ,'5MS(12)'   ,
     R          '3MNKS(12)' ,'4M2S(12)'  /
      DATA(LABLE(M),M=115,120)/'2NP(3)'    ,'ML(4)'     ,'2ML(6)'
     1,          'MSL(6)'    ,'2MSL(8)'   ,'3ML(8)'    /
      DATA(LABLE(M),M=121,128)/'TK(1)'     ,'RP(1)'     ,'KP(1)'
     S,          'THETA(1)'  ,'KJ(2)'     ,'OO(2)'     ,'MO(3)'
     T,          'SK(4)'     /                                              
      DATA(LABLE(M),M=129,138)/'MLN2S(2)'  ,'2ML2S(2)'  ,'MKL2S(2)'
     U,          '2MLS(4)'   ,'2NMLS(6)'  ,'2MLNS(6)'  ,'3MLS(6)'
     V,          '4MLS(8)'   ,'3MSL(10)'  ,'4MSL(12)'  /
      DATA(LABLE(M),M=139,140)/'2KO(1)'    ,'2OK(1)'    /                   
      DATA(LABLE(M),M=141,150)/'2KP(1)'   ,'2PK(1)'   ,'KP(2)'
     W,          '2SK(2)'   ,'2KS(2)'   ,'2TS(2)'   ,'ST(4)'
     X,          '3SK(4)'    ,'3KS(4)'    ,'3TS(4)'     /
      DATA(LABLE(M),M=151,164)/'SO(2)'    ,'SO(0)'    ,'.5MF'
     1,          '.5MSF'     ,'ST(0)'     ,'3SA'       ,'4SA'
     2,          '6SA'       ,'8SA'       ,'10SA'      ,'12SA'
     3,          '24SA'      ,'HS(3)'      ,'HS(5)'      /
      DATA(LABLE(M),M=165,175)/'O(2)'     ,'SK(2)'    ,'NK(3)'
     4,          'SP(3)'    ,'K(3)'     ,'NO(4)'    ,'MO(4)'
     5,         'SO(4)'     ,'S(7)'      ,'S(8)'      ,'S(10)'     /
    1 FORMAT(1H1)                                                         
    2 FORMAT( 10X,25H CONSTITUENT OF SPEED   ,F12.7,15H  NOT IN LIST   )  
    3 FORMAT(/// 10X, 30H**** EXECUTION TERMINATED ****)                  
    4 FORMAT(10X,13HCONSTITUENT  , A10, 15HNOT IN THE LIST  )             
    5 FORMAT(10X,15HCONSTITUENT NO., I5, 13H  NOT IN LIST )               
      GO TO (20,30,40),ICODE                                              
   20 DO 100 J = 1,175                                                    
      IF(SPDD.NE.SPD(J)) GO TO 100                                        
      ITAG = LABLE(J)                                                     
      ISUB = IP(J)                                                        
      INUM = J                                                            
      GO TO 101                                                           
  100 CONTINUE                                                            
      PRINT 1                                                             
      PRINT 2,SPDD                                                        
      PRINT 3                                                             
      CALL PEXIT                                                          
   30 DO 200 I = 1,175                                                    
      IF(ITAG.NE.LABLE(I)) GO TO 200                                      
      SPDD = SPD(I)                                                       
      ISUB = IP(I)                                                        
      INUM = I                                                            
      GO TO 101                                                           
  200 CONTINUE                                                            
      PRINT 1                                                             
      PRINT 4, ITAG                                                       
      PRINT 3                                                             
      CALL PEXIT                                                          
   40 DO 300 K = 1,175                                                    
      IF(INUM.NE.K) GO TO 300                                             
      ITAG = LABLE(K)                                                     
      SPDD = SPD(K)                                                       
      ISUB = IP(K)                                                        
      GO TO 101                                                           
  300 CONTINUE                                                            
      PRINT 1                                                             
      PRINT 5, INUM                                                       
      PRINT 3                                                             
      CALL PEXIT                                                          
  101 RETURN                                                              
      END                                                                 
C  DETERMINATION OF ORBITAL ELEMENTS FOR BEGINNING OF CENTURY             
      SUBROUTINE ORBIT(XCEN,XSX,XPX,XHX,XP1X,XNX,OEX,T,XYER,NNN)          
      IMPLICIT REAL*4(A-H,O-Z)
      DIMENSION OEX(NNN)                                                  
      S = 13.1763968
      P = 0.1114040
      XH = 0.9856473
      P1 = 0.0000471
      XN = -.0529539
      XCAN = XYER*0.01 + 0.001                                            
      XCEN = AINT(XCAN)*100.0                                             
      T = -3.0                                                            
      YR = 2.5                                                            
      GAT = 1600.0                                                        
      DO 10 JK = 1,30                                                     
      GP = GAT/400.0 + 0.00001
      COL = GP - AINT(GP)
      IF(COL.LT.0.010) GO TO 11                                           
      IF(GAT.EQ.XCEN) GO TO 12                                            
      YR = YR - 1.0                                                       
      GO TO 9                                                             
   11 IF( GAT.EQ.XCEN) GO TO 12                                           
    9 GAT = GAT + 100.0                                                   
   10 CONTINUE                                                            
   12 T = (GAT - 1900.0)*0.01                                             
      OEX(1) = 270.437422 + 307.892*T + 0.002525*T**2 + .00000189*T**3 +
     1 YR*S
      OEX(2) = 334.328019 + 109.032206*T - 0.01034444*T**2 - .0000125*T*
     1*3 + YR*P
      OEX(3) = 279.696678 + 0.768925*T + .0003205*T**2 + YR*XH
      OEX(4) = 281.220833 + 1.719175*T + 0.0004528*T**2 + .00000333*T**3
     1 + YR*P1
      OEX(5) = 259.182533 - 134.142397*T + .00210556*T**2 + .00000222*T*
     1*3 + YR*XN
      CALL TWOPI(OEX,5)
      DO 100 I = 1,5
  100 OEX(I) = FLOAT(IFIX(OEX(I)*100.0 + 0.5))*0.01
      XSX = OEX(1)
      XPX = OEX(2)
      XHX = OEX(3)
      XP1X = OEX(4)
      XNX = OEX(5)
      RETURN
       END
      SUBROUTINE OUTPUT (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * EPOCH,HSUBN,ICNTL,                    INOBS,    K,KAPPA,
     * KSTORE,MCNST,   M2,      NCONST,NDEXM1,NDEXP1,NDEXP2,
     * NDEXY,NDR,NORDER,NYYYY,OMEAN,OSD,OVAR,  RT,R1,SETS,SIG,
     * SX,THBAR,TIMEZ,TZERO,XLONG)
C
C_______________________________________________________________________
C
C     1 DEC 83
C
C     LABEL:  -OUTPUTT-
C
C     THIS IS A TEMPORARY FILE IN WHICH THE INDICATED CHANGES PERMIT
C     CALCULATION OF REGRESSION COEFFICIENTS FROM THE CORRELATION
C     MATRIX.  THIS IS EXPECTED TO IMPROVE STABILITY BY DECREASING THE
C     FREQUENCY OF UNDERFLOW.
C
C_______________________________________________________________________
C
C  CREATED JULY 1983
C
C  PURPOSE:
C        THIS SUBROUTINE -OUTPUT- SETS UP OUTPUT FIELDS FOR PRINTING.
C  IT BEGINS WITH THE CALCULATIONS OF VALUES NEEDED FOR THE HEADER
C  LINE.
C
C        SUBROUTINE
C        VARIABLES
C
C        ( ADD THEM HERE )
C
C
C        - - - DIMENSION STATEMENTS - - -
C
      REAL*8 TEMP,A,B,FPOCH,CONV,SPEED,DELTT
      REAL   PP(M2,M2)
      REAL SX(1),HSUBN(1),KAPPA(1),
     * EPOCH(1),AHSUBN(1),AEPOCH(1),COLA(1),
     *      R1(1),RT(1),SIG(1)
      INTEGER INOBS(1),MCNST(1),NORDER(1),ICNTL(1)
      CHARACTER*10 LABEL,IDENS,LABLE,TKIND,SKIND,LABL
      CHARACTER*4 BANAME(8)*4
      COMMON/PARAMC/RAMP(180),RKAPP(180),LABLE(180),RKAP(180)
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(200000,2)
     1       ,RATIO                                             
      COMMON/PARAM4/IYY(20),IMM(20),IDD(20),IHH(20),IMIN(20)
      COMMON/PARAM5/YEPOCH(180),YNODE(180),SPEED(180),P(180)
      COMMON/PARAM6/IDENS(16),LABEL(180),TKIND,SKIND
C
C        SUBROUTINE BEGINS HERE
C
C---------------------------------------------------------
C
C       STEP 1
C
C------------------------------------------------------------
C
      DO 100 II = 1,180
      RAMP(II) = 0.0
      RKAP(II) = 0.0
  100 RKAPP(II) = 0.0
      DZERO=DBLE(TZERO)
      MVAR = PP(1,NDEXP2)
      EMEAN = THBAR + SX(NDEXY)
      CONV = 180.0D0/3.14159265358979D0
C     DO 150 N = 1, K
C     NX = PP(N,NDEXP1)
C     A = PP(N,NDEXY)*SX(NX)
C     IF(ICNTL(5) .EQ. 1) A=A*SIG(NDEXY)/SIG(N)
C 150 EMEAN = EMEAN - A
      IF(MVAR.NE.2) GO TO 151
C 150 EMEAN = EMEAN - PP(K+1,NDEXY)*SX(NX+1)
  151 TEMP=PP(NDEXY,NDEXY)
      IF(ICNTL(5) .EQ. 0) THEN
      CVAR = SNGL(TEMP)/SETS
      ELSE
      CVAR = SNGL(TEMP)*OVAR
      END IF
C
C      IF(CVAR.LT.TZERO) PRINT 154, K,N,PP(NDEXY,NDEXY),SETS,CVAR,
C     * PP(K+1,NDEXP1)
C 154  FORMAT(1H0,'CVAR < TZERO',2I5,4F10.5,/)
      IF(CVAR.LT.0.) CVAR = 0.
      ESD = SQRT(CVAR)
C
C     PRINT OUTPUT TABLE HEADER
      NSETS = SETS
c     PRINT 183
      PRINT 185, TKIND, SKIND
  185 FORMAT(/1X, 2A10 )
      PRINT 153, NSETS
  153 FORMAT (' Analysis based on',I6,' observation points as indicated 
     *below',/)
      PRINT 157, (INOBS(I),DELTT,IYY(I),IMM(I),IDD(I),IHH(I),IMIN(I),
     * I=1,NDR)
c     PRINT 183
      PRINT 185, TKIND, SKIND
      PRINT 190, (IDENS(I),I = 1,16)
  157 FORMAT (1X,I5,' observations at intervals of',F5.2,' hours  ',
     * ' Year=',I5,' Month=',I3,' Day=',I3,' Hour=',I2,' Min=',I3)
      PRINT 162, OMEAN,OVAR,OSD,EMEAN,CVAR,ESD
  162 FORMAT( ' Observed Mean = ', F8.3,' units  Observed Variance = ',
     * F9.3,' sq units  Observed Standard Deviation = ', F8.3,' units',/
     *,' Constant in regression = ', F8.3,' units  Residual Variance = '
     *, F9.3,' sq units  Residual Standard Deviation = ',F8.3,' units'/)
      J = 0
      MVAR = 1
      PRINT 172
  172 FORMAT('                 ---- Adjusted for a standard year ---- * 
     *', 'Constituent   R**2 this   R**2 thru   Selec-  * Unadjusted val
     *ues',/, 1X,56HConstituent        (H)      (K)     (K'- K)     (K')
     *   *, '    Speed     constituent  this level   tion   *    (R)    
     *  (Z)',/,' Num. Label       (units) (degrees) (degrees) (degrees) 
     **  (deg./hr.)    ','only (1)  screening(2) Number  *  (units) (deg
     *rees)' / )
C
C*********************************************************************
C
C      STEP 2
C
C*********************************************************************
C
      DO 212 N = 1,K
      IF(MVAR .EQ. 2) THEN
      MVAR = 1
      GO TO 212
      ELSE
      J = J + 1
      TEMP = PP(N,NDEXP1)
      TERMN=SNGL(TEMP)
      NTERM=INT(TERMN)
      IF(PP(N,NDEXP2).EQ.2.000) THEN
C     THE PREDICTOR JUST SELECTED IS PERIODIC
      MVAR = 2
C
      IF(ICNTL(5) .EQ. 0) THEN
      A = PP(N,NDEXY)
      B = PP(N+1,NDEXY)
      ELSE
C_______________________________________________________________________
C
C     THE FOLLOWING INSERT IS INTENDED TO ELIMINATE
C     DIVIDE CHECK AND OVERFLOW PROBLEMS IN THE FOLLOWING
C     CALCULATIONS.
C
      IF (ABS(PP(N,NDEXY)) .LT. DZERO .OR. ABS(PP(N+1,NDEXY))
     * .LT. DZERO .OR. SIG(N) .LT. TZERO .OR. SIG(N+1) .LT. TZERO)
     * THEN
      PRINT 400, J,K,N,NDEXY
  400 FORMAT(1H ,' J=',I4,' K=',I4,' N=',I4,' NDEXY=',I4,/)
C
      PRINT 410, PP(N,NDEXY),PP(N+1,NDEXY),SIG(N),SIG(N+1)
  410 FORMAT (1H ,' POTENTIAL ERROR IN SUBROUTINE OUTPUT ',/,
     * ' PP(N,NDEXY),PP(N+1,NDEXY),SIG(N),SIG(N+1) ARE ',/,
     * 4E14.6,/)
C
      A = 1.0D0
      B = 1.0D0
      GO TO 420
      END IF
C
C_______________________________________________________________________
      A1=SIG(NDEXY)/SIG(N)
      B1=SIG(NDEXY)/SIG(N+1)
      A=PP(N,NDEXY)*DBLE(A1)
      B=PP(N+1,NDEXY)*DBLE(B1)
      END IF
  420 TEMP = DSQRT(A*A + B*B)
      HSUBN(J)=SNGL(TEMP)
      FPOCH = DATAN2(A,B)
      EPOCH(J) = FPOCH*CONV
C
      TEMP = PP(N+1,NDEXP1)
      CNSTN=SNGL(TEMP)
      NCNST=INT(CNSTN)/2
      MCNST(J) = NCNST
      GEPOCH = EPOCH(J)
      CALL TWOPI(GEPOCH,1)
      EPOCH(J) = GEPOCH
C
C     MAKE ADJUSTMENTS FOR A STANDARD YEAR, IF THESE ARE WANTED.
C  THAT IS TO SAY, APPLY EQUILIBRIUM ARGUMENTS AND NODE FACTORS.
      IF ( NYYYY .EQ. 0 ) GO TO 196
      AHSUBN(J) = HSUBN(J) / YNODE(NCNST)
      AEPOCH(J) = EPOCH(J) + YEPOCH(NCNST)
C * * * * * * * * * * *
C  TEMPORARY PRINT STATEMENT
C
      IF (AEPOCH(J) .LT. 0.) PRINT 500,LABEL(NCNST),AEPOCH(J)
  500 FORMAT(1X,'TEMPORARY PRINT OF ',A10,' WITH AEPOCH =',
     *      1X,F9.2,' BEFORE ADDING 360 DEGREES.')
C
C * * * * * * * * * *
      ZEPOCH = AEPOCH(J)
      CALL TWOPI(ZEPOCH,1)
      AEPOCH(J) = ZEPOCH
196   TEMP=SPEED(NCNST)
      CALL NNAME(TEMP,LABL,ISUB,INUMX,1)
      LABLE(J) = LABL
      COLA(J) = P(NCNST)*XLONG - SNGL(TEMP)*TIMEZ/15.
      DOLA = COLA(J)
      COLA(J) = AMOD( COLA(J),360. )
      KAPPA(J) = AEPOCH(J) - COLA(J)
      TKAPP = KAPPA(J)
      CALL TWOPI(TKAPP,1)
      KAPPA(J) = TKAPP
      PRINT 175,NCNST,LABLE(J),    AHSUBN(J),KAPPA(J),COLA(J),AEPOCH(J),
     * SPEED(NCNST),R1(J),RT(J),NORDER(J),HSUBN(J),EPOCH(J)
  175 FORMAT( 1X,I3,2X,A10,F9.4, 3(1X,F9.2),' * ', F12.7,2X,F9.6,3X,
     * F9.6, I7, 4X, '*', F9.4, F8.2 )
      ELSE
C  THE PREDICTOR JUST SELECTED IS NON-PERIODIC.
      TEMP = PP(N,NDEXY)
      HSUBN(J)=SNGL(TEMP)
C
C     IF(ICNTL(5) .EQ. 1) HSUBN(J) = HSUBN(J)
C
      TEMP=PP(N,NDEXP1)
      CNST=SNGL(TEMP)
      NCNST=INT(CNST)-NCONST
      MCNST(J) = NCNST
      PRINT 177, NCNST,BANAME(NCNST),R1(J),RT(J),NORDER(J),
     * HSUBN(J),EPOCH(J)
  177 FORMAT( 1X,I3,2X,A4,11X,'NON-PERIODIC VARIABLE',12X,'*',12X,
     *2(4X,F9.6),I7,'     *', 3X,F9.4, F9.2 )
      END IF
      END IF
 212  CONTINUE
      PRINT 180
  180 FORMAT(' (1) R1 is the square of the multiple correlation coeffi',
     *'cient between the tide and this constituent.',/,' (2) RT is the',
     *' square of the multiple correlation coefficient between the tid',
     *'e and all constituents thus far selected.',/,'     Each number ',
     * 'is printed with at least one digit which is believed to be ',
     * 'insignificant.'/)
C     *********   TO BE CHECKED OUT LATER    **************
      IF(K .EQ. 1) GO TO 200
C
C*********************************************************************
C
C               STEP 3
C
C*********************************************************************
C
C     REARRANGE THE OUTPUT IN STANDARD ORDER.
C
      IF (KSTORE .LT. 27) GO TO 192
c     PRINT 183
  183 FORMAT(1H1)
C     PRINT 185, TKIND, SKIND
C 186 PRINT 190, (IDENS(I),I = 1,16)
  190 FORMAT( 1X, 8A10/1X,8A10 /)
      PRINT 153, NSETS
      PRINT 157, (INOBS(I),DELTT,IYY(I),IMM(I),IDD(I),IHH(I),IMIN(I),
     * I=1,NDR)
  192 DUMMY = 0.0
      IPEN = 0
C
      DO 195 N = 1,NDEXM1
      DO 193 L = 1,NDEXM1
      IF(MCNST(L).NE.N) GO TO 193
      IF(PP(N,NDEXP2) .EQ.2.) THEN
C     PREDICTOR IS PERIODIC
      NCNST = MCNST(L)
      TEMP = SPEED(NCNST)
      CALL NNAME(TEMP,LABL,ISUB,INUMX,1)
      LABEL(L) = LABL
      CALL POSITX(TEMP,LOCK,L,1)
      RKAP(LOCK) = KAPPA(L)
      RAMP(LOCK) = AHSUBN(L)
      RKAPP(LOCK) = AEPOCH(L)
      IF(IPEN.GT.0) GO TO 1920
      PRINT 183
      PRINT 185, TKIND, SKIND
      PRINT 190, (IDENS(I),I = 1,16)
      PRINT 162, OMEAN, OVAR, OSD, EMEAN, CVAR, ESD
      PRINT 172
1920  PRINT 175,MCNST(L),LABLE(L),    AHSUBN(L),KAPPA(L),COLA(L),
     * AEPOCH(L),SPEED(NCNST),R1(L),RT(L),NORDER(L),HSUBN(L),
     * EPOCH(L)
      IPEN = IPEN + 1
      IF(IPEN.EQ.37) IPEN = 0
      ELSE
C     PREDICTOR IS NON PERIODIC
      NCNST = PP(N,NDEXP1) - NCONST
      PRINT 177, NCNST,BANAME(NCNST),R1(L),RT(L),NORDER(L),HSUBN(L),
     * EPOCH(L)
      END IF
      GO TO 195
193   CONTINUE
  195 CONTINUE
      CALL CPUNCH(RAMP,RKAPP,175,IDENS,KAZI,AZI,KINDAT)
  200 RETURN
      END
      SUBROUTINE PEXIT
C
      STOP 'PREMATURE EXIT FROM LSQHA3 DURING EXECUTION'
      END
      SUBROUTINE PMATRX (PP,M2,NDEXM1,NDEXY,NDEXP1,NSTEP,MATYP)
C
C_________________________________________________________________________
C
C  WRITTEN OCTOBER 1983
C
C  PURPOSE:
C
C     SUBROUTINE -PMATRX- IS CALLED BY SUBROUTINE -SCREEN-.
C  IT PRINTS THE CROSS PRODUCT MATRIX FORMED BY -LSQHA2- OR AS RE-
C  ARRANGED BY -SCREEN-, 10 COLUMNS TO THE PAGE AS NEEDED TO DISPLAY THE
C  MATRIX IN STANDARD FORM.  IT IS ANTICIPATED THAT THIS SUBROUTINE
C  WILL ONLY BE USED FOR TRAINING OR TROUBLESHOOTING.
C
C          MATYP = 0 - ORIGINAL MATRIX
C                = 1 - REARRANGED MATRIX
C                = 2 - SOLVED FOR ONE VARIABLE
C                = 3 - SOLVED FOR THE SECOND VARIABLE
C_____________________________________________________________________________
C
C  .  .  .  DIMENSION AND DECLARATION STATEMENTS  .  .  .
C
      REAL    PP(M2,M2)

C___________________________________________________________________________
C
C  .  .  . SUBROUTINE BEGINS HERE  .  .  .
C__________________________________________________________________________
C
C     DETERMINE THE NUMBER OF PAGES NEEDED.
      NPAGE = NDEXY/10 + 1
      IF ( NDEXM1 .LT. NDEXY ) NPAGE = NDEXM1/10 + 1
      IF ( NDEXM1 .EQ. 10 ) NPAGE = 1
C
      J1 = -9
      DO 40 J = 1,NPAGE
      J1 = J1 + 10
      J2 = J1 + 9
      IF ( J2 .GT. NDEXM1 ) J2 = NDEXM1
C
C        PRINT PAGE HEADER
C
      PRINT 10, NSTEP,MATYP,J,NPAGE
   10 FORMAT (1H1,'CROSS PRODUCT MATRIX, STEP',I4,' MATRIX TYPE',I4,
     * '   PAGE NO.',I3,'  OF',I3,/)
C
C        LABEL COLUMNS.
C
      PRINT 15,(I,I=J1,J2)
   15 FORMAT (1H0,'LINE',120X,'VAR.',/,' NO.',5X,I4,9(8X,I4),
     * 5X,'NO.',/)
      DO 30 M = 1,NDEXY
      PRINT 20, M,(PP(M,N),N=J1,J2),PP(M,NDEXP1)
   20 FORMAT (1X,I3,10F12.5,F5.0)
   30 CONTINUE
   40 CONTINUE
      RETURN
      END

      SUBROUTINE POSITX(SPEDX,IPOS,NUMP,ICODE)
      REAL*8 A(180),SPEDX
      DATA (A(I),I=1,48)/             28.9841042D0,  30.0000000D0,    
     1  28.4397295D0,  15.0410686D0,  57.9682084D0,  13.9430356D0,    
     2  86.9523127D0,  44.0251729D0,  60.0000000D0,  57.4238337D0,    
     3  28.5125831D0,  90.0000000D0,  27.9682084D0,  27.8953548D0,    
     4  16.1391017D0,  29.4556253D0,  15.0000000D0,  14.4966939D0,    
     5  15.5854433D0,   0.5443747D0,   0.0821373D0,   0.0410686D0,    
     6   1.0158958D0,   1.0980331D0,  13.4715145D0,  13.3986609D0,    
     7  29.9589333D0,  30.0410667D0,  12.8542862D0,  14.9589314D0,    
     8  31.0158958D0,  43.4761563D0,  29.5284789D0,  42.9271398D0,    
     9  30.0821373D0, 115.9364169D0,  58.9841042D0,  12.9271398D0,    
     1  14.0251729D0,  14.5695476D0,  15.9748272D0,  16.0569644D0,    
     2  30.5443747D0,  27.4238337D0,  28.9019669D0,  29.0662415D0,    
     3  26.8794590D0,  26.9523126D0/                                  
      DATA (A(I),I=49,96)/            27.4966873D0,  31.0980331D0,    
     1  27.8039338D0,  28.5947204D0,  29.1483788D0,  29.3734880D0,    
     2  30.7086493D0,  43.9430356D0,  45.0410686D0,  42.3827651D0,    
     3  59.0662415D0,  58.4397295D0,  57.4966873D0,  56.9523127D0,    
     4  58.5125831D0,  56.8794590D0,  59.5284789D0,  71.3668693D0,    
     5  71.9112440D0,  73.0092770D0,  74.0251728D0,  74.1073100D0,    
     6  72.9271398D0,  71.9933813D0,  72.4649023D0,  88.9841042D0,    
     7  86.4079380D0,  87.4238337D0,  87.9682084D0,  85.3920421D0,    
     8  85.8635632D0,  88.5125831D0,  87.4966873D0,  89.0662415D0,    
     9  85.9364168D0,  86.4807916D0,  88.0503457D0, 100.3509735D0,    
     1 100.9046318D0, 101.9112440D0, 103.0092771D0, 116.4079380D0,    
     2 116.9523127D0, 117.9682084D0, 114.8476674D0, 115.3920422D0,    
     3 117.4966873D0, 115.4648958D0/                                  
      DATA (A(I),I=97,124)/          116.4807916D0, 117.0344500D0,    
     1 118.0503457D0, 129.8887360D0, 130.4331108D0, 130.9774855D0,    
     2 131.9933813D0, 144.3761464D0, 144.9205211D0, 145.3920422D0,    
     3 145.9364169D0, 146.4807916D0, 146.9523127D0, 160.9774855D0,   
     4 174.3761464D0, 174.9205211D0, 175.4648958D0, 175.9364169D0,   
     5  14.9178647D0,  15.0821353D0,  15.1232059D0,  15.5125897D0,   
     6  30.6265119D0,  27.3416965D0,  42.9271397D0,  60.0821373D0,   
     7  16.1391016D0,  12.8450026D0/                                 
      DATA(A(I),I=125,140)/           26.9615963D0,  27.5059710D0,
     A  28.6040041D0,  57.5059710D0,  58.5218668D0,  85.4013258D0,
     B  85.9457005D0,  86.4900752D0,  87.5059710D0,  88.5218668D0,
     C 115.4741794D0, 116.4900752D0, 117.5059710D0, 146.4900752D0,
     D 175.4741794D0,  41.9205276D0/
      DATA(A(I),I=141,150)/           15.1232058D0,  14.8767942D0,
     E  30.0000001D0,  29.9178627D0,  30.1642746D0,  29.9178666D0,
     F  59.9589333D0,  59.9178627D0,  60.2464119D0,  59.8767999D0/
      DATA(  A(I),I=151,164)/ 28.9430356D0, 01.0569644D0, 00.5490165D0,
     1          00.5079479D0, 00.0410667D0, 00.1232059D0, 00.1642746D0,
     2          00.2464118D0, 00.3285841D0, 00.4106864D0, 00.4928237D0,
     3          00.9856473D0, 45.0000000D0, 75.0000000D0/
      DATA(  A(I),I=165,175)/ 27.8860712D0, 30.0410686D0, 43.4807981D0,
     4          44.9589314D0, 45.1232059D0, 56.3258007D0, 56.8701754D0,
     5          57.8860712D0,105.0000000D0,120.0000000D0,150.0000000D0/
    1 FORMAT(// 10X,32HNO MATCHING CONSTITUENT OF SPEED, F14.7 )
      GO TO (20,30),ICODE
   20 DO 100 I = 1,175
      IPOS = I
      IF(SPEDX.EQ.A(I)) GO TO 105
  100 CONTINUE
      PRINT 1
      CALL PEXIT
   30 DO 200 II = 1,175
      SPEDX = A(II)
      IF(II.EQ.NUMP) GO TO 105
  200 CONTINUE
  105 RETURN
      END
      SUBROUTINE SCREEN (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * CUTOFF,EPOCH,HSUBN,ICNTL,                    INOBS,    JR,
     * K,KAPPA,MAXTRM,MCNST,MINTRM,   M2,      NCOL,NCONST,NDEXM1,
     * NDEXP1,NDEXP2,NDEXY,NDR,NORDER,NQP,     NV,NYYYY,OMEAN,OSD,
     * OVAR,  QQ,RT,R1,SETS,SIG,SX,THBAR,TIMEZ,TSTORE,TZERO,
     * XLONG)
C_______________________________________________________________________
C
C     DATE:      13 DEC 83
C
C     LABEL:     -SCREN3T- + -SCREN4T-
C
C     PURPOSE:   THIS IS A COMPLETELY REWRITTEN FORM OF -SCREEN- USING
C                -SCREN2T- AS A BASE. REWRITING AND CORRECTIONS BY DLH.
C
C_______________________________________________________________________
C
C  MAJOR MODIFICATION AND DOCUMENTATION, MAY 1983.
C
C     PURPOSE:
C     SUBROUTINE -SCREEN- USES THE CROSS PRODUCT MATRIX OR THE
C  CORRELATION MATRIX GENERATED BY -CSTAT2- TO GENERATE A REGRESSION
C  EQUATION RELATING AN ARBITRARY TABULAR FUNCTION TO THE SUM OF A
C  TRIGONOMETRIC SERIES OF SPECIFIED FREQUENCIES, AND OF OTHER
C  SPECIFIED PREDICTORS.
C
C     THE SUBPROGRAM OPERATES IN THREE BASIC STEPS.  1) THE SQUARE
C  OF THE MULTIPLE CORRELATION OF TWO SUCCESSIVE PREDICTORS (THE SINE
C  AND COSINE OF A COMMON FREQUENCY) AND THE DEPENDENT VARIABLE IS
C  COMPUTED FOR EACH FREQUENCY, AS THE PREDICTOR HAVING THE HIGHEST
C  CORRELATION TO THE DEPENDENT VARIABLE IS SELECTED FOR
C  PROCESSING.  2) THE CROSS PRODUCT MATRIX IS MODIFIED BY
C  INTERCHANGING THE ROW AND COLUMN VECTORS CORRESPONDING TO
C  THE SELECTED FREQUENCY WITH THE UPPER LEFT ROW AND COLUMN
C  VECTORS OF THE UNSELECTED PORTION OF THE MATRIX.  3) THE
C  REGRESSION COEFFICIENTS EXPRESSING EACH UNSELECTED PREDICTOR
C  TO THE SELECTED PREDICTORS IS COMPUTED, AND THE UNSELECTED
C  PREDICTORS ARE MADE ORTHOGONAL TO THE SELECTED PREDICTORS.
C     THE PROCESS IS REPEATED UNTIL ALL PREDICTORS HAVE BEEN CON-
C  SIDERED, OR UNTIL A PREDICTOR IS SELECTED WHICH DOES NOT MAKE
C  A SIGNIFICANT IMPROVEMENT IN THE PREDICTION.  THE CALCULATION
C  MAY BE TERMINATED AT THIS STAGE.  PARTIAL RESULTS MAY BE
C  WRITTEN ON DISK FOR STORAGE, OR PRINTED IF DESIRED.  A HARD
C  COPY PRINTOUT OF THE FINAL RESULTS MAY BE OBTAINED AT EACH
C  STEP, OR ONLY AT THE FINAL STEP AS DESIRED.
C     SUBROUTINE -SCREEN- IS QUITE GENERAL AND CAN BE USED FOR THE
C  TRIGONOMETRIC EXPANSION OF MANY TYPES OF DATA.
C
C
C
C   SUBROUTINE
C   VARIABLES
C
C     NQP(N)    = INDEX NUMBERS OF DEFECTIVE VARIABLES.
C     RMAX      = THE SQUARE OF THE HIGHEST CORRELATION COEFFICIENT
C                 OBTAINED AT ANY LEVEL IN THE SCREENING PROCESS.
C     PP(I,J)   = CROSS CORRELATION MATRIX.  PP(I,NDEXP1) CONTAINS
C                 THE INITIAL INDEX NUMBERS OF ALL VECTORS, WHERE
C                 I = 1,NDEXY.
C     RSQR      = SQUARE OF THE MULTIPLE CORRELATION BETWEEN THE
C                 RESIDUE OF X(NDEXY) AND ANY FREQUENCY AT ANY STAGE OF
C                 THE CALCULATIONS.
C     KEY=L     = THE INDEX NUMBER FOR A SELECTED CONSTITUENT.
C     CUTOFF    = THE MINIMUM FRACTIONAL REDUCTION OF THE VARIANCE
C                 CONSIDERED TO BE SIGNIFICANT.
C     TEMP      = A TEMPORARY STORAGE CELL.
C     F         = A DUMMY FACTOR, USUALLY A REGRESSION COEFFICIENT.
C     KSTORE    = THE NUMBER OF CONSTITUENTS EVALUATED AT ANY STAGE
C                 AND THE INDEX NUMBER OF STORAGE OF OUTPUT RESULTS.
C     MINTRM    = MINIMUM NUMBER OF TERMS WANTED IN AN OUTPUT EQUATION.
C     MAXTRM    = MAXIMUM NUMBER OF TERMS WANTED IN AN OUTPUT EQUATION.
C     OVAR      = OBSERVED VARIANCE.
C     IPRINT    = ZERO, PRINT UNSORTED RESULTS.
C               = 1, PRINT SORTED RESULTS.
C     EMEAN     = CONSTANT IN REGRESSION EQUATION.
C     NX        = A DUMMY INDEX.
C     ESD       = STANDARD DEVIATION OF THE RESIDUE.
C     SETS(NDEXP1) = TOTAL NUMBER OF OBSERVATIONS.  ONLY IF DATA SERIES
C                 IS ONE UNBROKEN SET IS SETS(NDEXP1) = NOBS.
C     INOBS     = NUMBER OF OBSERVATIONS READ AT ONE TIME.
C     IYY       = YEAR                         X
C     IMM       = MONTH                         X    FIRST OBS
C     IDD       = DAY                            >   IN EACH
C     IHH       = HOUR                          X    SERIES
C     IMIN      = MINUTE, IN TENTHS OF HOURS   X
C     NDR       = THE NUMBER OF PERIODS OF RECORD CONSIDERED IN
C                 AN ANALYSIS.
C     HSUBN(J)  = THE AMPLITUDE OF A TIDAL CONSTITUENT.
C     EPOCH(J)  = THE PHASE OF A TIDAL CONSTITUENT.
C     YNODE(J)  = NODE FACTOR TO MODIFY HSUBN TO A STANDARD YEAR.
C     YEPOCH(J) = EQUILIBRIUM ARGUMENTS TO MODIFY -EPOCH- TO A
C                 STANDARD YEAR.
C     COLA(J)   = CORRECTION FACTOR NEEDED TO CONVERT EPOCH TO ZERO
C                 LONGITUDE.
C     KAPPA(J)  = PHASE REFERRED TO ZERO LONGITUDE.
C     NAME(J)   = STANDARD SYMBOL FOR A CONSTITUENT.
C     SPEED(J)  = FREQUENCY IN DEGREES-PER-HOUR.
C     R1(J)     = MULTIPLE CORRELATION COEFFICIENT BETWEEN THE TIDE AND
C                 ANY ONE COMPONENT.
C     RT(J)     = THE MULTIPLE CORRELATION COEFFICIENT BETWEEN THE TIDE
C                 AND ALL CONSTITUENTS CONSIDERED IN THE -J- STAGE.
C     NORDER(J) = THE ORDER IN WHICH A CONSTITUENT IS SELECTED.
C     TSTORE(N) = TEMPORARY STORAGE OF A VECTOR, USED FOR
C                 TRANSFERING ROWS AND COLUMNS OF A MATRIX.
C_______________________________________________________________________
C
C     SUBROUTINE -SCREEN- STARTS HERE.
C_______________________________________________________________________
C
      PARAMETER(NZ=37)
      REAL*8 DELTT
      REAL   PP(M2,M2),QQ(M2,3)
      REAL AHSUBN(NZ),AEPOCH(NZ),KAPPA(NZ)
      CHARACTER*4 BANAME(8)*4
      DIMENSION ICNTL(NZ),NORDER(NZ)
      DIMENSION NQP(NZ),R1(NZ),RT(NZ),HSUBN(NZ),EPOCH(NZ),JR(NZ),
     * MCNST(NZ),COLA(NZ)
      DIMENSION TSTORE(NZ), SX(NZ), SIG(NZ)
      KP = 0
      KNP = 0
C_______________________________________________________________________
C
C      INITIALIZE CONSTANTS AND THE FIELD FOR STORING LIST OF
C      FUNCTIONS WITH ZERO VARIANCE.
C
      DO 10 K = 1,NV
        NQP(K) = 0
   10   CONTINUE
      CONV = 180./3.141592625
      NCHNG = 0
      MVAR = 1
      NNV = NV
      KRC = 0
      NSTOP = 0
C
      IF (ICNTL(7) .EQ. 0) GO TO 12
      IF (ICNTL(5) .EQ. 0) THEN
          CALL PMATRX(PP,M2,NDEXY,NDEXY,NDEXP1,0,0)
        ELSE
          CALL MATCOR(PP,CUTOFF,JR,M2,NDEXY,NDEXY,NDEXP1,NDEXY,0,0)
        END IF
   12 CONTINUE
C
      DO 500 K = 1,NDEXM1
      IF (MVAR .EQ. 2) THEN
C
C     THIS VALUE OF K PERTAINS TO THE SECOND MEMBER OF A PAIR.
C     CALCULATIONS FOR THIS K ARE COMPLETED.
C
      MVAR = 1
      GO TO 500
      END IF
C
C     FOR EACH VALUE OF K, THE CODE ENDING AT 500 WILL SELECT THE
C     REMAINING VARIABLE OR VARIABLE PAIR WHICH IS MOST HIGHLY
C     CORRELATED WITH THE RESIDUAL OF THE DEPENDENT VARIABLE.  THE
C     MULTIPLE REGRESSION EQUATION WILL BE COMPUTED AND PRINTED AT
C     SELECTED STEPS.
C
C*********************************************************************
C
C              STEP 2 - A
C
C*********************************************************************
C
      RMAX = 0
      IF (ICNTL(5) .EQ. 0) THEN
C     SUM OF PRODUCTS MATRIX IS BEING USED,
C     OTHERWISE, THE CORRELATION MATRIX IS BEING USED.
C
C*********************************************************************
C
C              STEP 2 - AP - 1
C
C*********************************************************************
C
      DO 60 L = K,NNV,2
      MVAR = PP(L,NDEXP2)
      IF (MVAR .LT. 0) GO TO 60
C     THIS VARIABLE HAS BEEN DISCARDED, IT ADDS NO NEW INFORMATION.
      IF (MVAR .NE. 2) GO TO 70
      IF (PP(L,L) .LT. TZERO .AND. PP(L+1,L+1) .LT. TZERO) THEN
C     THIS VARIABLE SHOULD BE DISCARDED, IT ADDS NO USEFUL INFORMATION.
      PP(L,NDEXP2) = -PP(L,NDEXP1)
      PP(L+1,NDEXP2) = -PP(L,NDEXP1)-1.
      GO TO 60
C
      ELSE
C
C     ELIMINATE ANY POSSIBLE ATTEMPT TO DIVIDE BY ZERO.
      IF (PP(L,L).GT.TZERO) GO TO 30
      NQP(NCHNG) = INT(PP(L,NDEXP1))
      PP(L,L) = TZERO
      NCHNG = NCHNG + 1
   30 IF (PP(L+1,L+1) .GT. TZERO) GO TO 40
      NQP(NCHNG) = INT(PP(L+1,NDEXP1))
      PP(L+1,L+1) = TZERO
      NCHNG = NCHNG + 1
   40 CONTINUE
C
C     FORM THE SQUARE OF THE JOINT CORRELATION COEFFICIENT -RSQR-
C     BETWEEN THE DEPENDENT VARIABLE AND A SINE AND COSINE PAIR.
C
c     RSQR = ((PP(L,NDEXY)*PP(L,NDEXY)) / (PP(L,L)*PP(NDEXY,NDEXY))
c    *  + (PP(L+1,NDEXY)*PP(L+1,NDEXY)) / (PP(L+1,L+1)*PP(NDEXY,NDEXY))
c    *  - 2.*PP(L,L+1) * PP(L,NDEXY) * PP(L+1,NDEXY) /( PP(L,L)
c    *  *PP(L+1,L+1) * PP(NDEXY,NDEXY)))/ (1. - (PP(L,L+1)*PP(L,L+1)
c    *  /( PP(L,L) * PP(L+1,L+1))))
      RSQR = ((PP(L,NDEXY)*PP(L,NDEXY)) / (SIG(L)*SIG(NDEXY)*SETS)**2
     *  + (PP(L+1,NDEXY)*PP(L+1,NDEXY)) / (SIG(L+1)*SIG(NDEXY)*SETS)**2
     *  - 2.*PP(L,L+1) * PP(L,NDEXY) * PP(L+1,NDEXY) /((SIG(L)
     *  *SIG(L+1) * SIG(NDEXY)*SETS)**2*SETS))/ (1.-(PP(L,L+1)*PP(L,L+1)
     *  /( SIG(L) * SIG(L+1)*SETS)**2))
C
      IF(K .GT. NCOL) GO TO 50
      J = (L+1)/2
      R1(J) = RSQR
 50   IF ( RSQR .LE. RMAX ) GO TO 60
      RMAX = RSQR
      KEY = L
      MVAR = 2
      END IF
C
   60 CONTINUE
C
C*********************************************************************
C
C              STEP 2 - AP - 2
C
C     FORM SQUARES OF THE SIMPLE CORRELATION COEFFICIENT -RSQR-
C     BETWEEN THE DEPENDENT VARIABLE AND THE OTHER VARIABLES.
C
C*********************************************************************
C
   70 DO 90 L = NNV+1,NDEXM1
C
      IF (PP(L,L) .LT. TZERO) THEN
C     THIS VARIABLE SHOULD BE DISCARDED.
C     IT CONTAINS NO USEFUL INFORMATION.
      PP(L,NDEXP2) = -PP(L,NDEXP1)
      ELSE
c     RSQR = (PP(L,NDEXY)*PP(L,NDEXY))/(PP(L,L)*PP(NDEXY,NDEXY))
      RSQR = (PP(L,NDEXY)*PP(L,NDEXY))/(SIG(L)*SIG(NDEXY)*SETS)**2
      END IF
C
      IF (K .GT. NCOL) GO TO 80
      J = L- NCONST
      R1(J) = RSQR
   80 IF (RSQR .LT. RMAX) GO TO 90
      RMAX = RSQR
      KEY = L
      MVAR = 1
   90 CONTINUE
C
      ELSE
C
C*********************************************************************
C
C              STEP 2 - AC - 1
C
C     THE CORRELATION MATRIX IS BEING USED
C
C**********************************************************************
C
      DO 120 L = K,NNV,2
      MVAR = PP(L,NDEXP2)
      IF (MVAR .LT. 0) GO TO 120
C     THIS VARIABLE HAS BEEN DISCARDED.
C     IT CONTAINS NO NEW INFORMATION.
C
      IF (MVAR .NE. 2) GO TO 130
C
C     THIS VARIABLE IS A MEMBER OF A PAIR.
C     ELIMINATE ANY POSSIBLE ATTEMPT TO DIVIDE BY ZERO
C
      IF (PP(L,L+1)**2 .GE.1.) THEN
C     THIS VARIABLE SHOULD BE DISCARDED.
C     BOTH MEMBERS OF THE PAIR ARE IDENTICAL
C
      PP(L,NDEXP2) = -PP(L,NDEXP1)
      PP(L+1,NDEXP2) = -PP(L,NDEXP1)-1.
      GO TO 120
C
      ELSE
C
      RSQR = (PP(L,NDEXY)**2+PP(L+1,NDEXY)**2-
     * 2*PP(L,L+1)*PP(L,NDEXY)*PP(L+1,NDEXY))/(1.0-
     * PP(L,L+1)**2)
C
      IF (K .GT. NCOL) GO TO 100
      J = (L+1)/2
      R1(J) = RSQR
C
  100 IF (RSQR .LE. RMAX) GO TO 120
      RMAX = RSQR
      KEY = L
      MVAR = 2
      END IF
C
  120 CONTINUE
C_______________________________________________________________________________
C
C     FORM SQUARE OF THE SIMPLE CORRELATION COEFFICIENT BETWEEN
C     THE DEPENDENT VARIABLE AND OTHER VARIABLES.
C
  130 DO 150 L = NNV+1,NDEXM1
      RSQR = PP(L,NDEXY)**2
      IF (K .GE. NCOL) GO TO 140
      J = L - NCONST
      R1(J) = RSQR
C
  140 IF (RSQR .LT. RMAX) GO TO 150
      RMAX = RSQR
      KEY = L
      MVAR = 1
  150 CONTINUE
C
      END IF
C
C*********************************************************************
C
C              STEP 2 - B
C
C*********************************************************************
C
      IF(MVAR.EQ.1) THEN
      KNP = KNP + 1
      ELSE
      KP = KP + 1
      END IF
C_______________________________________________________________________
C
C     RECORD ORDER OF SELECTION OF PREDICTOR
      KSTORE = KP + KNP
      NORDER(KSTORE) = KSTORE
C     IF(K .GT. NCOL) GO TO 170
C     PRINT 160, (J,R1(J),J=(J+1)/2,NDEXY-NCONST-1)
  160 FORMAT(1H0,' TABLE OF R SQUARES',/,10(I3,F9.5))
C
C 170 CONTINUE
      IF (RMAX .LT. CUTOFF) THEN
C     NO REMAINING VARIABLE MAKES A SIGNIFICANT
C     CONTRIBUTION TO PREDICTIVE SKILL
C
      IF (NSTOP .EQ. 0) THEN
C     THIS IS THE FIRST VARIABLE TO FAIL THE CUTOFF TEST
C
      IF (KSTORE-1 .LT. MINTRM .OR. KSTORE-1 .GT.MAXTRM) THEN
C     RESULTS OF THE LAST STEP WERE VALID BUT WERE
C     NOT PRINTED.  THE SOLUTION SHOULD BE PRINTED NOW.
C
      PRINT 180, KSTORE,RMAX,KSTORE-1
  180 FORMAT ('THE',I4,'TH CONSTITUENT ACCOUNTS FOR ONLY',F9.7,
     * ' OF THE VARIANCE.  CALCULATIONS FOR',I4,' CONSTITUENTS ARE',
     * ' GIVEN BELOW.',/,' RESULTS OBTAINED LATER IN THE PROGRAM ARE',
     * ' OF DOUBTFUL VALUE.')
C
C  GO DIRECTLY TO THE OUTPUT ROUTINE.  MATRIX OPERATIONS FOR THE
C  PRECEEDING STEP HAVE BEEN COMPLETED.
C
      CALL OUTPUT (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * EPOCH,HSUBN,ICNTL,                    INOBS,  K-1,KAPPA,
     * KSTORE,MCNST,   M2,      NCONST,NDEXM1,NDEXP1,NDEXP2,
     * NDEXY,NDR,NORDER,NYYYY,OMEAN,OSD,OVAR,  RT,R1,SETS,SIG,
     * SX,THBAR,TIMEZ,TZERO,XLONG )
C
      NSTOP = KSTORE
C
      ELSE
C     PRECEEDING RESULTS WERE PRINTED
C
      PRINT 190, KSTORE,RMAX,KSTORE-1
  190 FORMAT ( 'THE',I4,'TH CONSTITUENT ACCOUNTS FOR ONLY'
     * ,F9.7, 'OF THE VARIANCE.',/, ' RESULTS FOR THE',I4,'TH ',
     *'CONSTITUENTS HAVE BEEN PRINTED AND ARE THE FINAL RESULTS ',
     *'THAT ARE BELIEVED TO BE VALID',/)
      END IF
      ELSE
C
C**********************************************************************
C
C              STEP 2 - C
C
C*********************************************************************
C
      IF (ICNTL(9) .EQ. 0) GO TO 501
C     CALCULATIONS TERMINATED BY THIS TRANSFER, FOR RSQR .LT. CUTOFF
      END IF
C
      ELSE
C     RESULTS OF THIS STEP ARE VALID.  CONTINUE CALCULATIONS
C
      IF (ICNTL(10) .EQ. 1) GO TO 320
C     THIS TRANSFER AVOIDS REARRANGEMENT OF MATRIX
      IF (K .EQ. KEY) GO TO 320
C     REARRANGEMENT OF MATRIX NOT REQUIRED
C
      IF(MVAR .EQ. 2) THEN
C     REARRANGE THE MATRIX FOR THE SELECTION OF A  PERIODIC
C     PREDICTOR.
C
C     EXCHANGE COLUMN VECTORS
      DO 210 M = 1,NDEXY
      TEMP = PP(M,K)
      PP(M,K) = PP(M,KEY)
      PP(M,KEY) = TEMP
      TEMP = PP(M,K+1)
      PP(M,K+1) = PP(M,KEY+1)
      PP(M,KEY+1) = TEMP
 210  CONTINUE
C     MOVE ROW VECTORS.
      DO 220 M = 1,NDEXP2
      TEMP = PP(K,M)
      PP(K,M) = PP(KEY,M)
      PP(KEY,M) = TEMP
      TEMP = PP(K+1,M)
      PP(K+1,M) = PP(KEY+1,M)
      PP(KEY+1,M) = TEMP
  220 CONTINUE
C     EXCHANGE CORRELATION COEFFICIENTS
      TEMP = R1((K+1)/2)
      R1((K+1)/2) = R1((KEY+1)/2)
      R1((KEY+1)/2) = TEMP
      TEMP = SIG((K+1))
      SIG((K+1)) = SIG((KEY+1))
      SIG((KEY+1)) = TEMP
      TEMP = SIG((K))
      SIG((K)) = SIG((KEY))
      SIG((KEY)) = TEMP
      ELSE
C     REARRANGE THE MATRIX FOR THE SELECTION OF A
C     NON PERIODIC PREDICTOR.
      DO 230 N =1,NDEXP2
C     SAVE SELECTED ROW VECTOR
  230 TSTORE(N) = PP(KEY,N)
C     MOVE EACH ROW VECTOR FROM K TO L-1 DOWN ONE LINE
C     BEGINING WITH L-1.
      DO 250 M = KEY,K+1,-1
      DO 240 N = 1,NDEXP2
 240  PP(M,N) = PP(M-1,N)
 250  CONTINUE
C     PLACE VECTORS IN ROW K
      DO 270 N = 1,NDEXP2
 270  PP(K,N) = TSTORE(N)
C     SAVE SELECTED COLUMN VECTORS
      DO 280 N=1,NDEXY
 280  TSTORE(N) = PP(N,KEY)
C     MOVE EACH COLUMN VECTOR FROM K TO L-1 ONE VALUE
C     TO THE RIGHT BEGINING WITH L-1.
      DO 300 M = KEY,K+1,-1
      DO 290 N = 1,NDEXY
 290  PP(N,M) = PP(N,M-1)
  300 CONTINUE
C     PLACE SELECTED COLUMN VECTORS IN COLUMN K
      DO 310 N = 1,NDEXY
 310  PP(N,K) = TSTORE(N)
C     REARRANGE CORRELATION COEFFICIENTS
      TEMP = R1(KEY-NCONST)
      DO 311 M=KEY-NCONST,K+1-NCONST,-1
 311  R1(M) = R1(M-1)
      R1(K-NCONST) = TEMP
      TEMP = SIG(KEY)
      DO 312 M=KEY,K+1,-1
 312  SIG(M) = SIG(M-1)
      SIG(K) = TEMP
C     ESTABLISH NEW LOCATION FOR FINAL PERIODIC VARIABLES.
      NNV = NV + 1
      END IF
      END IF
C
C_________________________________________________________________________
      IF (ICNTL(7) .EQ. 0 .OR. K .GT. NCOL) GO TO 315
      IF (ICNTL(5) .EQ. 0) THEN
          CALL PMATRX(PP,M2,NCOL,NDEXY,NDEXP1,K,1)
        ELSE
          CALL MATCOR(PP,CUTOFF,JR,M2,NCOL,NDEXY,NDEXP1,NDEXY,K,1)
        END IF
  315 CONTINUE
C_______________________________________________________________________________
C
C     AT THIS POINT -KEY- INDICATES THE CURRENT POSITION OF THE SINE
C  FUNCTION CORRESPONDING TO THE SELECTED FREQUENCY ON THE POSITION
C  OF THE SELECTED NON-PERIODIC VARIABLE.
C
C     IN THE NEXT LOOP, THE CROSS PRODUCTS IN THE FIRST ROW OF THE
C  MATRIX ARE REPLACED BY THE REGRESSION COEFFICIENTS WHICH BEST EXPRESS
C  EACH OF THE UNSELECTED VARIABLES IN TERMS OF THE SELECTED VARIABLES.
C  IN ALL ROWS, THE PORTION OF THE VARIANCE EXPLAINED BY THE SELECTED
C  VARIABLE IS REMOVED FROM THE REMAINING VARIANCE.  THIS OPERATION IS
C  FIRST CARRIED OUT FOR THE SINE TERM AND THEN FOR THE COSINE TERM.
C
C_______________________________________________________________________
C
C     SOLVE FOR COEFFICIENT OF SINE OR NON PERIODIC FUNCTION.
 320  DO 340 L = 1,NDEXY
      IF (L.EQ.K) GO TO 340
      FKL = PP(L,K) / PP(K,K)
      DO  330 M = K,NDEXY
  330 PP(L,M) = PP(L,M) - FKL * PP(K,M)
  340 CONTINUE
      FKK =  PP(K,K)
      DO 350 M = K,NDEXY
 350  PP(K,M) = PP(K,M)/FKK
      IF(MVAR.NE.2) GO TO 390
C
C  FOLLOWING STATEMENT FOR TESTING.
C
      IF (ICNTL(7) .EQ. 0 .OR. K .GT. NCOL) GO TO 355
      IF (ICNTL(5) .EQ. 0) THEN
          CALL PMATRX(PP,M2,NCOL,NDEXY,NDEXP1,K,2)
        ELSE
          CALL MATCOR(PP,CUTOFF,JR,M2,NCOL,NDEXY,NDEXP1,NDEXY,K,2)
        END IF
  355 CONTINUE
C
C     SOLVE FOR THE COEFFICIENT OF THE COSINE TERM.
C
      DO 370 L = 1,NDEXY
      IF ( L.EQ.K+1 ) GO TO 370
      FKL = PP(L,K+1) / PP(K+1,K+1)
      DO 360 M = K+1,NDEXY
  360 PP(L,M) = PP(L,M) - FKL * PP(K+1,M)
  370 CONTINUE
      FKK =  PP(K+1,K+1)
      DO 380 M = K+1,NDEXY
  380 PP(K+1,M) = PP(K+1,M)/FKK
  390 CONTINUE
C
      IF (ICNTL(7) .EQ. 0 .OR. K .GT. NCOL) GO TO 405
      IF (ICNTL(5) .EQ. 0) THEN
          CALL PMATRX(PP,M2,NCOL,NDEXY,NDEXP1,K,3)
        ELSE
          CALL MATCOR(PP,CUTOFF,JR,M2,NCOL,NDEXY,NDEXP1,NDEXY,K,3)
        END IF
  405 CONTINUE
C________________________________________________________________________
C
C
C  END OF ABOVE LOOP
C_______________________________________________________________________
C
C     NOTE . . . NOTE . . . NOTE  THE ABOVE INSTRUCTIONS MUST BE
C  PRESERVED COMPLETELY AND IN ORDER
C
C_______________________________________________________________________
C
C     AT THE END OF THIS LOOP -PP(N,NDEXY)- CONTAINS THE REGRESSION
C  COEFFICIENT FOR THE FIRST (KSTORE) INDEPENDENT VARIABLES SELECTED
C  OR THE ASSOCIATED VALUE OF REGRESSION COEFFICIENT *
C  SIG(NDEXY)/SIG(K).
C
C  -PP(NV,NV)- CONTAINS THE RESIDUAL VARIATION OF THE DEPENDENT DATA,
C  I.E., THE PART OF THE VARIANCE NOT RELATED TO THE FIRST (KSTORE)
C  INDEPENDENT VARIABLES.
C
C_______________________________________________________________________
C
      IF (PP(NDEXY,NDEXY) .LT. TZERO) THEN
      PRINT 410, KSTORE
  410 FORMAT ( ' CALCULATIONS FOR THE ',I4,'TH CONSTITUENT ARE',
     * ' UNSTABLE, AS SIXSQR(NDEXY) .LT. ZERO.',/
     * ' CALCULATIONS ARE BEING TERMINATED. ',/)
C
      IF (KSTORE-NSTOP .EQ. 1) THEN
      PRINT 420, KSTORE-1
  420 FORMAT( ' THE FINAL VALID RESULTS WERE OBTAINED AT STEPS '
     * ,I4,' AND ARE GIVEN ABOVE. ')
      GO TO 501
      ELSE
      DO 430 MM = 1,NDEXY
C     RESTORES FIELDS NEEDED FOR PRINTING THE OUTPUT OF PREVIOUS STEP
      PP(MM,NDEXY) = QQ(MM,1)
      PP(MM,NDEXP1) = QQ(MM,2)
      PP(MM,NDEXP2) = QQ(MM,3)
  430 CONTINUE
C
      CALL OUTPUT (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * EPOCH,HSUBN,ICNTL,                    INOBS,    K,KAPPA,
     * KSTORE,MCNST,   M2,      NCONST,NDEXM1,NDEXP1,NDEXP2,
     * NDEXY,NDR,NORDER,NYYYY,OMEAN,OSD,OVAR,  RT,R1,SETS,SIG,
     * SX,THBAR,TIMEZ,TZERO,XLONG)
C
      NSTOP = KSTORE
      GO TO 501
C     THIS TRANSFER WILL TERMINATE CALCULATIONS
      END IF
      ELSE
C
C*********************************************************************
C
C              STEP 2 - E
C
C*********************************************************************
C
C     RESULTS OBTAINED AT THIS STEP ARE VALID. STORE DATA NEEDED
C     FOR RECOVERY IF NEXT STEP FAILS.
C
      DO 440 MM = 1,NDEXY
      QQ(MM,1) = PP(MM,NDEXY)
      QQ(MM,2) = PP(MM,NDEXP1)
      QQ(MM,3) = PP(MM,NDEXP2)
  440 CONTINUE
      END IF
C
C
C     STORE MULTIPLE CORRELATION COEFFICIENT FOR THIS STEP.
C
      IF(ICNTL(5) .EQ. 0) THEN
      RT(KSTORE) = 1.0 - (PP(NDEXY,NDEXY) / (OVAR*SETS))
      ELSE
      RT(KSTORE) = 1.0 - PP(NDEXY,NDEXY)
      END IF
C
C     PRINT 450, KSTORE, RT(KSTORE)
C 450 FORMAT(1H0,' KSTORE=',I3,' RT= ',F7.4,/)
C
C     TEST TO DETERMINE IF OUTPUT IS WANTED AT THIS STEP.  IF NOT,
C  GO TO THE END OF THE LOOP.
      IF(NSTOP.EQ.KSTORE.AND.ICNTL(9).EQ.0) GO TO 500
      IF ( K+2 .GE. NDEXM1 ) GO TO 460
C
C     A PRINTOUT OF RESULTS MAY BE OBTAINED FOR THE LAST
C  COMPUTATION STEP COMPLETED, BUT, THE RESULTS OF THE FIRST
C  (MINTRM-1) LOOPS MAY BE OMITTED IF COMPUTATIONS ARE CONTINUED.
C
      IF ( KSTORE .LT. MINTRM ) GO TO 500
      IF ( KSTORE .GT. MAXTRM .AND. K .LT. NDEXM1-1) GO TO 500
  460 CONTINUE
C
C     IF TRANSFER TO 500 IS NOT MADE, CALL OUTPUT
C
      CALL OUTPUT (PP,      BANAME,AEPOCH,AHSUBN,COLA,DELTT,
     * EPOCH,HSUBN,ICNTL,                    INOBS,    K,KAPPA,
     * KSTORE,MCNST,   M2,      NCONST,NDEXM1,NDEXP1,NDEXP2,
     * NDEXY,NDR,NORDER,NYYYY,OMEAN,OSD,OVAR,  RT,R1,SETS,SIG,
     * SX,THBAR,TIMEZ,TZERO,XLONG)
C
      KAY = K
  500 CONTINUE
C
C*********************************************************************
C
C              STEP 3
C
C*********************************************************************
C
C     SOLUTION COMPLETED
C     PRINT IDENTIFIERS FOR EACH DROPPED VARIABLE
C
  501 IF ( NCHNG .EQ. 0 ) GO TO 525
C
      IF (KSTORE .LT. NCONST) PRINT 503, CUTOFF
  503 FORMAT( ' VARIABLES FOR WHICH NUMBERS ARE LISTED BELOW WERE ',
     *'DROPPED TO AVOID INSTABILITY OR',/,'  BECAUSE EACH VARIABLE ACC',
     *'OUNTED FOR LESS THEN ',F9.7,' OF THE INITIAL VARIANCE. ',/)
      PRINT 507, (PP(I,NDEXP1),I=KAY,NDEXP1)
  507 FORMAT(1H ,20F4.0,/)
C
      J = 0
      DO 510 I = 1,NDEXM1
      IF ( NQP(I) .EQ. 0 ) GO TO 510
      J = J + 1
      NQP(J) = NQP(I)
  510 CONTINUE
C
      PRINT 520, (NQP(M),M=1,J)
  520 FORMAT (' THE FOLLOWING PREDICTORS WERE MODIFIED TO AVOID ',
     * 'DIVISION BY ZERO',/,'  CALCULATIONS INVOLVING THESE PREDIC',
     * 'TORS ARE OF DOUBTFUL VALUE.',/,(20I4))
C
  525 J = 0
      DO 530 I=1,NDEXM1
      IF (PP(I,NDEXP2) .GT. 0) GO TO 530
      J = J+1
      NQP(J) = INT(-PP(I,NDEXP2))
  530 CONTINUE
C
      PRINT 540, (NQP(I), I=1,J)
  540 FORMAT ( ' THE FOLLOWING VARIABLES WERE REJECTED TO AVOID ',
     * 'INSTABILITY',/,20I4,/)
C
      RETURN
      END
      SUBROUTINE SUMTRG(A1,FINAL,SINSUM,COSSUM)
      REAL*8 A1,A2,A3,COF,FINALL,FINAL
      FINALL = FINAL
      A2 = FINALL*A1
      A3 = A2 + A1
      CALL TWOPIR(A1,1)
      CALL TWOPIR(A2,1)
      CALL TWOPIR(A3,1)
      COF = DSIN(A3)/DSIN(A1)
      SINSUM = COF*DSIN(A2)
      COSSUM = COF*DCOS(A2)
      RETURN
      END
      SUBROUTINE SUMTRGRAN(A1,NOBS,SINSUM,COSSUM,TIM,DELTT)
      REAL*8 A1,A2,A3,COF,TIM(1),DELTT
      SINSUM=0.
      COSSUM=0.
      CALL TWOPIR(A1,1)
      IF(TIM(1).GE.DELTT)THEN
      A2 = A1 * NINT((TIM(1) - DELTT)/DELTT)
      A3 = A2 + A1
      CALL TWOPIR(A2,1)
      CALL TWOPIR(A3,1)
      COF = DSIN(A3)/DSIN(A1)
      SINSUM = SINSUM - COF*DSIN(A2)
      COSSUM = COSSUM - COF*DCOS(A2)
      END IF
      DO L=2,NOBS
      IF(TIM(L)-TIM(L-1).GT.1.5*DELTT)THEN
      A2 = A1 * NINT(TIM(L-1)/DELTT)
      A3 = A2 + A1
      CALL TWOPIR(A2,1)
      CALL TWOPIR(A3,1)
      COF = DSIN(A3)/DSIN(A1)
      SINSUM = SINSUM + COF*DSIN(A2)
      COSSUM = COSSUM + COF*DCOS(A2)
      A2 = A1 * NINT((TIM(L) - DELTT)/DELTT)
      A3 = A2 + A1
      CALL TWOPIR(A2,1)
      CALL TWOPIR(A3,1)
      COF = DSIN(A3)/DSIN(A1)
      SINSUM = SINSUM - COF*DSIN(A2)
      COSSUM = COSSUM - COF*DCOS(A2)
      END IF
      END DO
      A2 = A1 * NINT(TIM(NOBS)/DELTT)
      A3 = A2 + A1
      CALL TWOPIR(A2,1)
      CALL TWOPIR(A3,1)
      COF = DSIN(A3)/DSIN(A1)
      SINSUM = SINSUM + COF*DSIN(A2)
      COSSUM = COSSUM + COF*DCOS(A2)
      RETURN
      END
      SUBROUTINE TABLE6(VI,V,XI,VP,VPP,CIG,CVX,CEX,PVC,PVCP,ANG,AN,AT)     
      COMMON/BOXS/AW,AI,AE,AE1,ASP                                         
      V = 0.0                                                              
      XI = 0.0                                                             
      VP = 0.0                                                             
      VPP = 0.0                                                            
      AN = ANG*0.0174533                                                   
      AX = ANG                                                             
      EYE = COS(AI)*COS(AW) - SIN(AI)*SIN(AW)*COS(AN)                      
      C9 = ACOS(EYE)*57.2957795                                            
      VI = FLOAT(IFIX(C9*100.0 + 0.5))*0.01                                
      CIG = VI*0.0174533                                                   
      AT = CIG                                                             
      IF(CIG.EQ.0.0) GO TO 230                                             
      IF(AX.EQ.0.0.OR.AX.EQ.180.0) GO TO 230                               
      VXXE = SIN(AI)*SIN(AN)                                               
      VXXN = COS(AI)*SIN(AW) + SIN(AI)*COS(AW)*COS(AN)                     
      IF(VXXE.EQ.0.0.OR.VXXN.EQ.0.0) GO TO 201                             
      VXX = VXXE/VXXN                                                      
      C10 = ATAN(VXX)*57.2957795                                           
      V = FLOAT(IFIX(C10*100.0 + 0.5))*0.01                                
      IF(AX.GT.180.0.AND.V.GT.0.0) V = -1.0*V                              
  201 CVX = V*0.0174533                                                    
      TERM = SIN(AI)*(COS(AW)/SIN(AW))                                     
      EXX = TERM*(SIN(AN)/COS(AN)) + (COS(AI) - 1.0)*SIN(AN)               
      IF(EXX.EQ.0.0) GO TO 202                                             
      EZZ = TERM + COS(AI)*COS(AN) + (SIN(AN)**2/COS(AN))                  
      IF(EZZ.EQ.0.0) GO TO 202                                             
      EXEZ = EXX/EZZ                                                       
      IF(EXEZ.GT.3450.0) GO TO 202                                         
      C11 = ATAN(EXEZ)*57.2957795                                          
      XI = FLOAT(IFIX(C11*100.0 + 0.5))*0.01                               
      IF(AX.GT.180.0.AND.XI.GT.0.0) XI = -1.0*XI                           
  202 CEX = XI*0.0174533                                                   
      A22 = (0.5 + 0.75*AE**2)*SIN(2.0*CIG)                                
      B22 = (0.5 + 0.75*AE1**2)*SIN(2.0*AW)*ASP                            
      VPXE = A22*SIN(CVX)                                                  
      VPXN = A22*COS(CVX) + B22                                            
      IF(VPXE.EQ.0.0.OR.VPXN.EQ.0.0) GO TO 203                             
      VPX = VPXE/VPXN                                                      
      IF(VPX.GT.3450.0) GO TO 203                                          
      VP = ATAN(VPX)*57.2957795                                            
      IF(AX.GT.180.0.AND.VP.GT.0.0) VP = -1.0*VP                           
  203 PVC = VP*0.0174533                                                   
      A47 = (0.5 + 0.75*AE**2)*SIN(CIG)**2                                 
      B47 = (0.5 + 0.75*AE1**2)*ASP*SIN(AW)**2                             
      VPYE = A47*SIN(2.0*CVX)                                              
      VPYN = A47*COS(2.0*CVX) + B47                                        
      IF(VPYE.EQ.0.0.OR.VPYN.EQ.0.0) GO TO 204                             
      VPY = VPYE/VPYN                                                      
      IF(VPY.GT.3450.0) GO TO 204                                          
      VPP = ATAN(VPY)*57.2957795                                           
      IF(AX.GT.180.0.AND.VPP.GT.0.0) VPP = -1.0*VPP                        
  204 PVCP = VPP*0.0174533                                                 
  230 RETURN                                                               
      END                                                                  
      SUBROUTINE  TEMP(XSER,XTEMP,I )                                      
      DIMENSION   XSER(I) , XTEMP(I)                                       
      DO 100 J = 1, I                                                      
  100 XTEMP(J) = XSER(J)                                                   
      RETURN                                                               
      END                                                                  
      SUBROUTINE TRGSA ( PP,FREQ,FHOUR,SX,NDR,INOBS,M2,NV)
C
C___________________________________________________________________________
C
C  WRITTEN MAY 1983
C     TEST VERSION 7/11/83
C
C_____________________________________________________________________________
C
C  PURPOSE:
C
C    SUBROUTINE -TRGSA2- COMPUTES THE SUMS AND SUMS-OF-PRODUCTS OF TRIG-
C    ONOMETRIC FUNCTIONS.  THE FOLLOWING RELATIONS ARE USED FOR THE SUM-
C    MATION:   SIN(NX) OR COS(NX), N = 1,M  ARE EQUAL TO
C       ?SIN((M+1)X/2) / SIN(X/2)?  TIMES  SIN(MX/2) OR COS(MX/2).
C
C    SIN(NX)**2, COS(NX)**2 AND SIN(NX)COS(NX) ARE EXPRESSED AS
C    1/2? 1 +OR- COS(2NX)? OR 1/2? SIN(2NX)?.  OTHER PRODUCTS ARE
C    EXPRESSED AS FUNCTIONS OF (X-Y) AND (X+Y).  A PSEUDO SUBROUTINE IS
C    USED FOR THE SUMMATION.
C
C_____________________________________________________________________________
C
C . . . . DIMENSION AND DECLARATION STATEMENTS . . . .
C
      REAL*8 FREQ(1),FHOUR(1),A1,HORS,HOUR,DELTT,DELTN,SHOUR,FINNAL
      REAL*8 DFREQ(180),DELTAX,FINAL
      REAL   PP(M2,M2),SX(M2)
      INTEGER INOBS(1)
      COMMON/BGROUP/DELTT,SHOUR(20),DELTAX
C____________________________________________________________________________
C
C  .  .  SUBROUTINE STARTS HERE  .  .
C
C_____________________________________________________________________________
C
C
C******************************************************************************
C
C     . . . STEP 1 . . .
C
C******************************************************************************
C
      DO 444 KKK = 1,NV,2
      KK = (KKK + 1)/2
  444 DFREQ(KK) = FREQ(KK)/DELTAX
      HOUR=SHOUR(1)
C
C****************************************************************************
C
C     . . . STEP 2 . . .
C
C****************************************************************************
C
      DO 300 L=1,NDR
      HORS = FLOAT(INOBS(L))
      FINAL = FHOUR(L)
      FINNAL = FINAL
      PRINT 2, L
    2 FORMAT (' Data Block No.',I3)
      IF (SHOUR(L+1) .GT.FINNAL+1.0D0) PRINT 7,FINNAL, SHOUR(L+1)
C  THERE IS A DATA BREAK BETWEEN THIS DATA BLOCK AND THE NEXT.
    7 FORMAT (' From TRGSA: Gap in data between ',F14.5,' and ',F14.5)
      ALPHA = 1.0
C
C************************************************************************
C
C     . . . STEP 3 . . .
C
C  START OF BASIC SUB-PROGRAM
C
C***************************************************************************
C
   20 DO 200 I = 1,NV,2
C     COMPUTE CONSTITUENT INDEX NUMBER
      J = (I+1)/2
C     COMPUTE SUMS OF SINES AND COSINES
      A1 = (DFREQ(J)/2.0D0)
      KEY = 1
      CALL SUMTRG(A1,FINAL,SINSUM,COSSUM)
   21 SX(I) = SX(I) + SINSUM * ALPHA
      SX(I+1) = SX(I+1) + COSSUM * ALPHA
C
C*****************************************************************************
C
C     . . . STEP 4 . . .
C
C     COMPUTE SUMS OF SQUARES AND PRODUCTS SIN(X)COS(X)
C
C****************************************************************************
C
      A1 = DFREQ(J)
      KEY = 2
      CALL SUMTRG(A1,FINAL,SINSUM,COSSUM)
   31 PP(I,I) = PP(I,I) + ((FINAL +1. - COSSUM)/2.) * ALPHA
      PP(I,I+1) = PP(I,I+1) + (SINSUM/2.) * ALPHA
      PP(I+1,I+1) = PP(I+1,I+1) + ((FINAL +1. + COSSUM)/2.) * ALPHA
      IF(I.GE.(NV-1)) GO TO 200
C
C*******************************************************************************
C
C     . . . STEP 5 . . .
C
C     COMPUTE SUMS OF CROSS PRODUCTS OF TRIGONOMETRIC TERMS WITH UNLIKE
C     ARGUMENTS, TAKING THE DIFFERENCE TERM FIRST,
C
C****************************************************************************
C
      M = I + 2
      DO 100 K=M,NV,2
      LL = (K+1)/2
      A1 = ((DFREQ(J) -DFREQ(LL))/2.0D0)
      KEY = 3
      CALL SUMTRG(A1,FINAL,SINSUM,COSSUM)
   41 P1 = SINSUM
      P2 = COSSUM
C     TAKE THE SUM OF ARGUMENTS
      A1 = ((DFREQ(J) +DFREQ(LL))/2.0D0)
      KEY = 4
      CALL SUMTRG(A1,FINAL,SINSUM,COSSUM)
   51 PP(I,K) = PP(I,K) + ((P2 - COSSUM)/2.) * ALPHA
      PP(I+1,K+1) = PP(I+1,K+1) + ((P2 + COSSUM)/2.)*ALPHA
      PP(I,K+1) = PP(I,K+1) + ((SINSUM + P1)/2.)*ALPHA
      PP(I+1,K) = PP(I+1,K) + ((SINSUM - P1)/2.)*ALPHA
  100 CONTINUE
  200 CONTINUE
C
C
C****************************************************************************
C
C     . . . STEP 6 . . .
C
C  IF THE SERIES DID NOT START AT THE FIRST SCHEDULED OBSERVATION OF THE
C  YEAR, SUBTRACT THE PORTION OF THE SUM DUE TO THE PERIOD BEFORE THE
C  BEGINNING OF THE SERIES.
C
C******************************************************************************
      DELTN =  0.0D0
      IF(HOUR.NE.DELTN) THEN
      FINAL = SHOUR(L)
      HOUR = DELTN
      ALPHA = -1.
      GO TO 20
      ELSE
C  THERE IS NO DATA GAP AT THIS POINT.  MATRIX CALCULATIONS NOT
C  NEEDED HERE.
      END IF
C
      IF (L .NE. NDR) GO TO 290
      GO TO 300
  290 HOUR = SHOUR(L+1)
C
C*****************************************************************************
C
C     . . . STEP 7 . . .
C
C******************************************************************************
C
  300 CONTINUE
C
C
      RETURN
      END
      SUBROUTINE TRGSARAN ( PP,FREQ,SX,NOBS,M2,NV,TIM)
C
C___________________________________________________________________________
C
C  WRITTEN December 1995 by Chris Zervas
C_____________________________________________________________________________
C
C  PURPOSE:
C
C    SUBROUTINE -TRGSA2- COMPUTES THE SUMS AND SUMS-OF-PRODUCTS OF TRIG-
C    ONOMETRIC FUNCTIONS.
C_____________________________________________________________________________
C
C . . . . DIMENSION AND DECLARATION STATEMENTS . . . .
C
      REAL*8 FREQ(1),A1,DELTT,SHOUR,TIM(1),DELTAX,DFREQ(180)
      REAL   PP(M2,M2),SX(M2)
      COMMON/BGROUP/DELTT,SHOUR(20),DELTAX
C____________________________________________________________________________
C
C  .  .  SUBROUTINE STARTS HERE  .  .
C
C_____________________________________________________________________________
C
      DO 444 KKK = 1,NV,2
      KK = (KKK + 1)/2
  444 DFREQ(KK) = FREQ(KK)/DELTAX
   20 DO 200 I = 1,NV,2
C     COMPUTE CONSTITUENT INDEX NUMBER
      J = (I+1)/2
C     COMPUTE SUMS OF SINES AND COSINES
      A1 = (DFREQ(J)/2.0D0)
      CALL SUMTRGRAN(A1,NOBS,SINSUM,COSSUM,TIM,DELTT)
   21 SX(I) = SX(I) + SINSUM
      SX(I+1) = SX(I+1) + COSSUM
C
C*****************************************************************************
C
C     COMPUTE SUMS OF SQUARES AND PRODUCTS SIN(X)COS(X)
C
C****************************************************************************
C
      A1 = DFREQ(J)
      CALL SUMTRGRAN(A1,NOBS,SINSUM,COSSUM,TIM,DELTT)
   31 PP(I,I) = PP(I,I) + ((NOBS - COSSUM)/2.)
      PP(I,I+1) = PP(I,I+1) + (SINSUM/2.)
      PP(I+1,I+1) = PP(I+1,I+1) + ((NOBS + COSSUM)/2.)
      IF(I.GE.(NV-1)) GO TO 200
C
C*******************************************************************************
C
C     COMPUTE SUMS OF CROSS PRODUCTS OF TRIGONOMETRIC TERMS WITH UNLIKE
C     ARGUMENTS, TAKING THE DIFFERENCE TERM FIRST,
C
C****************************************************************************
C
      M = I + 2
      DO 100 K=M,NV,2
      LL = (K+1)/2
      A1 = ((DFREQ(J) -DFREQ(LL))/2.0D0)
      CALL SUMTRGRAN(A1,NOBS,SINSUM,COSSUM,TIM,DELTT)
   41 P1 = SINSUM
      P2 = COSSUM
C     TAKE THE SUM OF ARGUMENTS
      A1 = ((DFREQ(J) +DFREQ(LL))/2.0D0)
      CALL SUMTRGRAN(A1,NOBS,SINSUM,COSSUM,TIM,DELTT)
   51 PP(I,K) = PP(I,K) + ((P2 - COSSUM)/2.)
      PP(I+1,K+1) = PP(I+1,K+1) + ((P2 + COSSUM)/2.)
      PP(I,K+1) = PP(I,K+1) + ((SINSUM + P1)/2.)
      PP(I+1,K) = PP(I+1,K) + ((SINSUM - P1)/2.)
  100 CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE TWOPI(AUG,NUM)
      DIMENSION AUG(NUM)
      DO 100 I = 1,NUM
      AZT = AUG(I)/360.0
      IF(AZT) 77,88,88
   77 AUG(I) = ((AZT - AINT(AZT)) + 1.00)*360.0
      GO TO 100
   88 AUG(I) = (AZT - AINT(AZT))*360.0
  100 CONTINUE
      RETURN
      END
      SUBROUTINE TWOPIR(AUG,IO)
      REAL*8 AUG(IO),AZT,PI,TOOPI
      PI = 3.14159265358979D0
      TOOPI = 2.0D0*PI
      DO 100 MO = 1,IO
       AZT = AUG(MO)/TOOPI
      IF(AZT) 77,88,88
   77 AUG(MO) = ((AZT - DINT(AZT)) + 1.00)*TOOPI
      GO TO 100
   88 AUG(MO) =  (AZT - DINT(AZT))*TOOPI
  100 CONTINUE
      RETURN
      END
      SUBROUTINE VANDUF(SPEED,E,F,ITYPE)
C     THIS SUBROUTINE COMPUTES THE EQUILIBRIUM ARGUMENT V(0) + U  AND
C     NODE FACTOR (F) OF THE NATIONAL OCEAN SERVICE SPECIAL PUBLICATION
C     NUMBER 98  --  PREPERED BY E. E. LONG, 1979 --
C     ****   ORDER OF CONSTISUENTS   ****
C     M(2),N(2),S(2),O(1),K(1),K(2),L(2),2N(2),R(2),T(2),LAMBDA(2),MU(2)
C     NU(2),J(1),M(1),OO(1),P(1),Q(1),2Q(1),RHO(1),M(4),M(6),M(8),S(4)
C     S(6),M(3),S(1),MK(3),2MK(3),MN(4),MS(4),2SM(2),MF,MSF,MM,SA,SSA
      REAL*8 SPD,SPEED
      DIMENSION SPD(180),MS(180)
      COMMON/LOCAT/TM,GONL
      COMMON/VEE/TML,CON,U,Q,UI
      COMMON/BOXA/S,XL,PM,PL,SL,PS,PLM,SKYN,VI,V,XI,VPP
      COMMON/BOXB/VP,P,AUL,AUM,CRA,CQA                                     
      COMMON/BOXS/AW,AI,AE,AE1,ASP                                         
      DATA(MS(N),N = 1,37)/3*2,2*1,8*2,7*1,4,6,8,4,6,3,1,2*3,2*4,2,5*0/    
      DATA(MS(N),N = 38, 48)/ 5*1,6*2/                                     
      DATA(MS(N),N =  49, 96)/7*2,3*3,7*4,8*5,12*6,4*7,7*8/                
      DATA(MS(N),N =  97,124)/3*8,4*9,6*10,11,4*12,4*1,2*2,3,4,2*1/        
      DATA(MS(N),N=125,140)/3*2,2*4,5*6,3*8,10,12,3/                       
      DATA(MS(N),N=141,150)/2*1,4*2,4*4/
      DATA(MS(N),N=151,164)/2,11*0,3,5/
      DATA(MS(N),N=165,175)/2*2,3*3,3*4,7,8,10/
      DATA(SPD(I),I=  1, 37)/ 28.9841042D0, 28.4397295D0, 30.0000000D0,    
     1    13.9430356D0,   15.0410686D0,   30.0821373D0,   29.5284789D0,    
     2    27.8953548D0,   30.0410667D0,   29.9589333D0,   29.4556253D0,    
     3    27.9682084D0,   28.5125831D0,   15.5854433D0,   14.4966939D0,    
     4    16.1391017D0,   14.9589314D0,   13.3986609D0,   12.8542862D0,    
     5    13.4715145D0,   57.9682084D0,   86.9523127D0,  115.9364169D0,     
     6    60.0000000D0,   90.0000000D0,   43.4761563D0,   15.0000000D0,     
     7    44.0251729D0,   42.9271398D0,   57.4238337D0,   58.9841042D0,     
     8    31.0158958D0,   01.0980331D0,   01.0158958D0,   00.5443747D0,     
     9    00.0410686D0,   00.0821373D0/                                     
      DATA(SPD(I),I= 38, 48)/ 12.9271398D0, 14.0251729D0, 14.5695476D0,     
     A    15.9748272D0,   16.0569644D0,   30.5443747D0,   27.4238337D0,     
     B    28.9019669D0,   29.0662415D0,   26.8794590D0,   26.9523126D0/     
      DATA(SPD(I),I= 49, 87)/ 27.4966873D0, 31.0980331D0, 27.8039338D0,     
     C    28.5947204D0,   29.1483788D0,   29.3734880D0,   30.7086493D0,     
     D    43.9430356D0,   45.0410686D0,   42.3827651D0,   59.0662415D0,     
     E    58.4397295D0,   57.4966873D0,   56.9523127D0,   58.5125831D0,     
     F    56.8794590D0,   59.5284789D0,   71.3668693D0,   71.9112440D0,     
     G    73.0092770D0,   74.0251728D0,   74.1073100D0,   72.9271398D0,     
     H    71.9933813D0,   72.4649023D0,   88.9841042D0,   86.4079380D0,     
     I    87.4238337D0,   87.9682084D0,   85.3920421D0,   85.8635632D0,     
     J    88.5125831D0,   87.4966873D0,   89.0662415D0,   85.9364168D0,     
     K    86.4807916D0,   88.0503457D0,  100.3509735D0,  100.9046318D0/     
      DATA(SPD(I),I= 88,124)/101.9112440D0,103.0092771D0,116.4079380D0,     
     L   116.9523127D0,  117.9682084D0,  114.8476674D0,  115.3920422D0,     
     M   117.4966873D0,  115.4648958D0,  116.4807916D0,  117.0344500D0,     
     N   118.0503457D0,  129.8887360D0,  130.4331108D0,  130.9774855D0,     
     O   131.9933813D0,  144.3761464D0,  144.9205211D0,  145.3920422D0,     
     P   145.9364169D0,  146.4807916D0,  146.9523127D0,  160.9774855D0,     
     Q   174.3761464D0,  174.9205211D0,  175.4648958D0,  175.9364169D0,     
     R    14.9178647D0,   15.0821353D0,   15.1232059D0,  15.5125897D0,      
     S    30.6265119D0,   27.3416965D0,   42.9271397D0,   60.0821373D0,     
     T    16.1391016D0,   12.8450026D0/                                     
      DATA(SPD(I),I=125,140)/ 26.9615963D0, 27.5059710D0, 28.6040041D0,     
     U    57.5059710D0,   58.5218668D0,   85.4013258D0,   85.9457005D0,     
     V    86.4900752D0,   87.5059710D0,   88.5218668D0,  115.4741794D0,     
     W   116.4900752D0,  117.5059710D0,  146.4900752D0,  175.4741794D0,     
     X    41.9205276D0/                                                     
      DATA(SPD(I),I=141,150)/ 15.1232058D0, 14.8767942D0, 30.0000001D0,
     Y    29.9178627D0,   30.1642746D0,   29.9178666D0,   59.9589333D0,
     Z    59.9178627D0,   60.2464119D0,   59.8767999D0/
      DATA(SPD(I),I=151,164)/ 28.9430356D0, 01.0569644D0, 00.5490165D0,
     1          00.5079479D0, 00.0410667D0, 00.1232059D0, 00.1642746D0,
     2          00.2464118D0, 00.3285841D0, 00.4106864D0, 00.4928237D0,
     3          00.9856473D0, 45.0000000D0, 75.0000000D0/
      DATA(SPD(I),I=165,175)/ 27.8860712D0, 30.0410686D0, 43.4807981D0,
     4          44.9589314D0, 45.1232059D0, 56.3258007D0, 56.8701754D0,
     5          57.8860712D0,105.0000000D0,120.0000000D0,150.0000000D0/
  500 FORMAT(1H0,49H*** (V + U) NOT COMPUTED FOR CONSTITUENT OF SPEED,     
     1F12.7, 7H   ****,22H  EXECUTION TERMINATED  )                        
      CON = SL + TML                                                       
      DO 600 J = 1,175
      IPICK = J
      IF(SPEED.EQ.SPD(J)) GO TO 611                                        
  600 CONTINUE                                                             
      IPICK = 176
  611 IF(IPICK.GT.164) GO TO 620
      IF(IPICK.GT.124) GO TO 618
      IF(IPICK.GT.88) GO TO 616                                            
      IF(IPICK.GT.38) GO TO 614                                            
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23   
     1,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38),IPICK                 
  614 IPIK = IPICK - 38                                                    
      GO TO (39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58   
     1,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80   
     2,81,82,83,84,85,86,87,88),IPIK                                       
  616 IPIKK = IPICK - 88                                                   
      GO TO (89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,10   
     16,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,   
     2123,124),IPIKK                                                       
  618 IPIKKK = IPICK - 124                                                 
      GO TO (125,126,127,128,129,130,131,132,133,134,135,136,137,138,139   
     1,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,
     2156,157,158,159,160,161,162,163,164),IPIKKK
  620 IPIIKK = IPICK - 164
      GO TO (165,166,167,168,169,170,171,172,173,174,175,488),IPIIKK
    1 E = 2.0*(CON - PM + XI - V)                                          
  201 F   = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))             
      GO TO 888                                                            
    2 E =  2.0*(CON + XI - V) - 3.0*PM + PL                                
      GO TO 201                                                            
    3 E = 2.0*TML                                                          
  203 F = 1.0                                                              
      GO TO 888                                                            
    4 E = CON - V - 2.0*(PM - XI) - 90.0                                   
      GO TO 218                                                            
    5 E = CON - VP + 90.0                                                  
  205 F = 1.0/SQRT(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0   
     1.1006)                                                               
      GO TO 888                                                            
    6 E = 2.0*CON - VPP                                                    
  206 F = 1.0/SQRT(19.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.0*U) + 0   
     1.0981)                                                               
      GO TO 888                                                            
    7 E = 2.0*CON - PM - PL + AUL                                          
  207 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))*(1.0/CRA)     
      GO TO 888                                                            
    8 E = 2.0*(CON + XI - V + PL) - 4.0*PM                                 
      GO TO 201                                                            
    9 E = SL - PS + 180.0 + 2.0*TML                                        
      GO TO 203                                                            
   10 E = 2.0*TML - (SL - PS)                                              
      GO TO 203                                                            
   11 E = 2.0*(CON + XI - V - SL) - PM + PL + 180.0                        
      GO TO 201                                                            
   12 E = 2.0*(CON + XI - V + SL) - 4.0*PM                                 
      GO TO 201                                                            
   13 E = 2.0*(CON + XI - V + SL) - 3.0*PM - PL                            
      GO TO 201                                                            
   14 E = CON + PM - PL - V + 90.0                                         
  214 F = (SIN(2.0*AW)*(1.0 - 1.5*SIN(AI)**2))/SIN(2.0*UI)                
      GO TO 888                                                           
   15 E = CON - PM + AUM                                                  
      F =  ((SIN(AW)*COS(0.5*AW)**2*COS(0.5*AI)**4)/(SIN(UI)*COS(0.5*UI)  
     1**2))*(1.0/CQA)                                                     
      GO TO 888                                                           
   16 E = CON - V + 2.0*(PM - XI)+ 90.0                                   
  216 F = (SIN(AW)*SIN(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*SIN(.5*UI)**2)   
      GO TO 888                                                           
   17 E = TML + 270.0 - SL                                                
      GO TO 203                                                           
   18 E = CON - V - 3.0*PM + 2.0*XI + PL - 90.0                           
  218 F = (SIN(AW)*COS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2)   
      GO TO 888                                                           
   19 E = CON - V - 4.0*PM + 2.0*XI + 2.0*PL - 90.0                       
      GO TO 218                                                           
   20 E = CON - V - 3.0*PM + 2.0*XI - PL + 2.0*SL - 90.0                  
      GO TO 218                                                           
   21 E = 4.0*(CON - PM + XI - V)                                         
  221 F   = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**2         
      GO TO 888                                                           
   22 E = 6.0*(CON - PM + XI - V)                                         
  222 F   = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**3         
      GO TO 888                                                           
   23 E = 8.0*(CON - PM + XI - V)                                         
  223 F   = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**4         
      GO TO 888                                                           
   24 E = 4.0*TML                                                         
      GO TO 203                                                           
   25 E = 6.0*TML                                                         
      GO TO 203                                                           
   26 E = 3.0*(CON - PM + XI - V) + 180.0                                 
  226 F = (COS(0.5*AW )**6*COS(0.5*AI)**6)/COS(0.5*UI)**6                 
      GO TO 888                                                           
   27 E = TML + 180.0                                                     
      GO TO 203                                                           
   28 E = 2.0*(CON - PM + XI - V) + ( CON - VP + 90.0)                    
  228 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**1*(1./SQRT  
     1(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0.1006))       
      GO TO 888                                                           
   29 E = 4.0*(CON - PM + XI - V) - (CON - VP + 90.0)                     
  229 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**2*(1./SQRT  
     1(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0.1006))       
      GO TO 888                                                           
   30 E = 4.0*(CON + XI - V) + PL - 5.0*PM                                
      GO TO 221                                                           
   31 E = 2.0*(CON - PM + XI - V) + 2.0*TML                               
      GO TO 201                                                           
   32 E = 4.0*TML - 2.0*(CON - PM + XI - V)                               
      GO TO 201                                                           
   33 E = 2.0*(PM - XI)                                                   
  233 F = (SIN(AW)**2*COS(0.5*AI)**4)/SIN(UI)**2                          
      GO TO 888                                                           
   34 E = 2.0*TML - 2.0*(CON - PM + XI - V)                               
      GO TO 201                                                           
   35 E = PM - PL                                                         
  235 F = ((2./3.- SIN(AW)**2)*(1.- 1.5*SIN(AI)**2))/(2./3.- SIN(UI)**2)  
      GO TO 888                                                           
   36 E = SL                                                              
      GO TO 203                                                           
   37 E = 2.0*SL                                                          
      GO TO 203                                                           
   38 E = 2.0*(CON - V - 2.0*(PM - XI) - 90.0) - (TML+270.0 - SL)         
  238 F = ((SIN(AW)*COS(0.5*AW)**2*COS(0.5*AI)**4)/(SIN(UI)*COS(0.5*UI)*  
     1*2))**2                                                             
      GO TO 888                                                           
   39 E = 2.0*(CON - PM + XI - V) - (TML + 270.0 - SL)                    
      GO TO 201                                                           
   40 E = CON + 2.0*SL - V - PM - PL + 90.0                               
      GO TO 258                                                           
   41 E = 2.0*(TML + 270.0 - SL) - (CON - V - 2.0*(PM - XI) - 90.0)       
      GO TO 218                                                           
   42 E = 2.0*TML - (CON - V - 2.0*(PM - XI) - 90.0)                      
      GO TO 218                                                           
   43 E = 2.0*TML + PM - PL                                               
      GO TO 221                                                           
   44 E = 4.0*(CON + XI - V) - 5.0*PM - 2.0*TML  + PL                     
      GO TO 221                                                           
   45 E = CON - V - 2.0*(PM - XI) + TML - SL + 180.0                      
      GO TO 218                                                           
   46 E = 4.0*CON - 2.0*(PM - XI + V + TML) - VPP                         
      GO TO 259                                                           
   47 E = 4.0*(CON + XI - V) - 6.0*PM + 2.0*(PL - TML)                    
      GO TO 221                                                           
   48 E = 6.0*(CON - PM) + 4.0*(XI - V - TML) + AUL                       
      GO TO 261                                                           
   49 E = 6.0*CON - 5.0*PM + 4.0*(XI - V - TML) - PL + AUL                
      GO TO 261                                                           
   50 E = 2.0*(TML + PM - XI + V) - VPP                                   
      GO TO 259                                                           
   51 E = 4.0*(XI - PM - V) + 2.0*(TML + VPP)                             
  251 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**2*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981)**2)         
      GO TO 888                                                           
   52 E = 6.0*CON - 4.0*TML - 3.0*PM + 2.0*(XI - V) - PL - VPP + AUL      
  252 CRA = SQRT(1.0 - 12.0*(SIN(.5*UI)**2/COS(.5*UI)**2)*COS(2.0*Q) + 3  
     16.0*(SIN(.5*UI)**4/COS(.5*UI)**4))                                  
      F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**2*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))*(1./CRA)   
      GO TO 888                                                           
   53 E = 6.0*CON - 4.0*TML + 2.0*(XI - V - PM - VPP)                     
  253 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**1*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981)**2)         
      GO TO 888                                                           
   54 E = 4.0*TML - 2.0*CON + PL - PM + VPP                               
  254 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**2*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))            
      GO TO 888                                                           
   55 E = 4.0*CON - 2.0*(TML + VPP) + PM - PL                             
      GO TO 251                                                           
   56 E = 2.0*TML + CON - V - 2.0*(PM - XI) - 90.0                        
      GO TO 218                                                           
   57 E = 2.0*TML + CON - VP + 90.0                                       
      GO TO 205                                                           
   58 E = 3.0*(CON - V) + 4.0*XI - 5.0*PM + PL - 90.0                     
  258 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**1*((SIN(AW)*C  
     1OS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))                
      GO TO 888                                                           
   59 E = 4.0*CON - 2.0*(PM - XI + V) - VPP                               
  259 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**1*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))            
      GO TO 888                                                           
   60 E = 2.0*(TML + CON + XI - V) - 3.0*PM + PL                          
      GO TO 201                                                           
   61 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) - 2.0*TML - PL + AUL            
  261 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**3*(1./SQRT(1.  
     10 - 12.0*(SIN(.5*UI)**2/COS(.5*UI)**2)*COS(2.0*Q) + 36.0*(SIN(.5*U  
     2I)**4/COS(.5*UI)**4)))                                              
      GO TO 888                                                           
   62 E = 6.0*(CON - PM + XI - V) - 2.0*TML                               
      GO TO 222                                                           
   63 E = 4.0*CON - 3.0*PM + 2.0*(XI - V) - PL + AUL                      
  263 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**2*(1./SQRT(1.  
     10 - 12.0*(SIN(.5*UI)**2/COS(.5*UI)**2)*COS(2.0*Q) + 36.0*(SIN(.5*U  
     2I)**4/COS(.5*UI)**4)))                                              
      GO TO 888                                                           
   64 E = 4.0*(CON + XI - V) - 6.0*PM + 2.0*PL                            
      GO TO 221                                                           
   65 E = 2.0*(CON + TML) - PM - PL + AUL                                 
      GO TO 207                                                           
   66 E = 5.0*(CON - V) - 7.0*PM + 6.0*XI + PL - 90.0                     
  266 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**2*((SIN(AW)*C  
     1OS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))                
      GO TO 888                                                           
   67 E = 5.0*(CON - V) + 6.0*(XI - PM) - 90.0                            
      GO TO 266                                                           
   68 E = 5.0*CON + 4.0*(XI - PM - V) - VP + 90.0                         
      GO TO 229                                                           
   69 E = 3.0*CON + 2.0*(XI - PM - V + TML) - VP + 90.0                   
      GO TO 228                                                           
   70 E = 5.0*CON + 2.0*(XI - PM - V) - (VP + VPP) + 90.0                 
  270 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**1*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))*(1./SQRT(  
     2.8965*SIN(2.*UI)**2 + .6001*SIN(2.*UI)*COS(U) + 0.1006))            
      GO TO 888                                                           
   71 E = 4.0*(CON - PM + XI - V) + TML - SL + 270.0                      
      GO TO 221                                                           
   72 E = 6.0*(CON - PM + XI - V) - TML + SL - 270.0                      
      GO TO 222                                                           
   73 E = 5.0*(CON - PM) + 4.0*(XI - V) + PL - VP + 90.0                  
      GO TO 229                                                           
   74 E = 2.0*(CON - PM + XI - V) + 4.0*TML                               
      GO TO 201                                                           
   75 E = 6.0*(CON + XI - V) - 7.0*PM + PL                                
      GO TO 222                                                           
   76 E = 4.0*(CON + XI - V) - 5.0*PM + 2.0*TML + PL                      
      GO TO 221                                                           
   77 E = 4.0*(CON - PM + XI - V) + 2.0*TML                               
      GO TO 221                                                           
   78 E = 8.0*CON - 9.0*PM + 6.0*(XI - V) - 2.0*TML + PL + AUL            
  278 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**4*(1./SQRT(1.  
     10 - 12.0*(SIN(.5*UI)**2/COS(.5*UI)**2)*COS(2.0*Q) + 36.0*(SIN(.5*U  
     2I)**4/COS(.5*UI)**4)))                                              
      GO TO 888                                                           
   79 E = 6.0*(CON + XI - V) - 8.0*PM + 2.0*PL                            
      GO TO 222                                                           
   80 E = 4.0*CON - 3.0*PM + 2.0*(XI - V + TML) - PL + AUL                
      GO TO 263                                                           
   81 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) - PL + AUL                      
      GO TO 261                                                           
   82 E = 4.0*CON + 2.0*(XI - PM - V + TML) - VPP                         
      GO TO 259                                                           
   83 E = 8.0*(CON - PM) + 6.0*(XI - V) - 2.0*TML + AUL                   
      GO TO 278                                                           
   84 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) - 2.0*TML - PL + AUL            
      GO TO 278                                                           
   85 E = 6.0*CON + 4.0*(XI - PM - V) - VPP                               
      GO TO 254                                                           
   86 E = 7.0*(CON - V) - 9.0*PM + 8.0*XI + PL - 90.0                     
  286 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**3*((SIN(AW)*C  
     1OS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))                
      GO TO 888                                                           
   87 E = 7.0*CON + 6.0*(XI - V) - 8.0*PM + 2.0*PL - VP + 90.0            
  287 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**3*(1./SQRT  
     1(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0.1006))       
      GO TO 888                                                           
   88 E = 5.0*(CON - V) + 6.0*(XI - PM) + 2.0*TML - 90.0                  
      GO TO 266                                                           
   89 E = 5.0*CON + 4.0*(XI - PM) - 3.0*V + 2.0*TML - VPP - 90.0          
  289 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**1*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))*( (SIN(AW  
     2)*COS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))             
      GO TO 888                                                           
   90 E = 6.0*CON - 7.0*PM + 6.0*(XI - V) + 2.0*TML + PL                  
      GO TO 222                                                           
   91 E = 6.0*(CON - PM + XI - V) + 2.0*TML                               
      GO TO 222                                                           
   92 E = 4.0*(CON - PM + XI - V) + 4.0*TML                               
      GO TO 221                                                           
   93 E = 8.0*(CON + XI - V) - 10.0*PM + 2.0*PL                           
      GO TO 223                                                           
   94 E = 8.0*(CON + XI - V) - 9.0*PM + PL                                
      GO TO 223                                                           
   95 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) + 2.0*TML - PL + AUL            
      GO TO 261                                                           
   96 E = 10.0*CON - 9.0*PM + 8.0*(XI - V) - 2.0*TML - PL + AUL           
  296 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**5*(1./SQRT(1.  
     10 - 12.0*(SIN(.5*UI)**2/COS(.5*UI)**2)*COS(2.0*Q) + 36.0*(SIN(.5*U  
     2I)**4/COS(.5*UI)**4)))                                             
      GO TO 888                                                          
   97 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) - PL + AUL                     
      GO TO 278                                                          
   98 E = 8.0*CON + 6.0*(XI - PM - V) - VPP                              
  298 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**3*(1./SQRT(19 
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))           
      GO TO 888                                                          
   99 E = 6.0*CON + 4.0*(XI - PM - V) + 2.0*TML - VPP                    
      GO TO 254                                                          
  100 E = 9.0*CON - 10.0*PM + 8.0*(XI - V) + 2.0*PL - VP + 90.0          
  300 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))**4*(1./SQRT 
     1(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0.1006))      
      GO TO 888                                                          
  101 E = 9.0*(CON - PM) + 8.0*(XI - V) + PL - VP + 90.0                 
      GO TO 300                                                          
  102 E = 9.0*CON + 8.0*(XI - PM - V) - VP + 90.0                        
      GO TO 300                                                          
  103 E = 7.0*CON + 6.0*(XI - PM - V) + 2.0*TML - VP + 90.0              
      GO TO 287                                                          
  104 E = 10.0*(CON + XI - V) - 11.0*PM + PL                             
  304 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**5             
      GO TO 888                                                          
  105 E = 10.0*(CON - PM + XI - V)                                       
      GO TO 304                                                          
  106 E = 8.0*(CON + XI - V) - 9.0*PM + PL + 2.0*TML                     
      GO TO 223                                                          
  107 E = 8.0*(CON - PM + XI - V) + 2.0*TML                              
      GO TO 223                                                          
  108 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) - PL + AUL                     
      GO TO 278                                                          
  109 E = 6.0*(CON - PM + XI - V) + 4.0*TML                              
      GO TO 222                                                          
  110 E = 9.0*CON + 8.0*(XI - PM - V) + 2.0*TML - VP + 90.0              
      GO TO 300                                                          
  111 E = 10.0*(CON + XI - V) - 11.0*PM + PL                              
      GO TO 304                                                           
  112 E = 10.0*(CON - PM + XI - V) + 2.0*TML                              
      GO TO 304                                                           
  113 E = 10.0*CON - 9.0*PM + 8.0*(XI - V) + 2.0*TML - PL + AUL           
      GO TO 296                                                           
  114 E = 8.0*(CON - PM + XI - V) + 4.0*TML                               
      GO TO 223                                                           
  115 E = TML - 2.0*SL + PS - 90.0                                        
      GO TO 203                                                           
  116 E = CON + SL - PS + 90.0                                            
      GO TO 203                                                           
  117 E = CON + 2.0*SL + 90.0                                             
      GO TO 203                                                           
  118 E = TML + PM - SL + PL - V + 90.0                                   
      GO TO 258                                                           
  119 E = 2.0*CON + PM - V - VP - PL + 180.0                              
  319 F =  ((SIN(2.*AW)*(1.0 - 1.5*SIN(AI)**2))/SIN(2.0*UI))*(1./SQRT(.8  
     1965*SIN(2.*UI)**2 + .6011*SIN(2.*UI)*COS(U) + 0.1006))              
      GO TO 888                                                           
  120 E = 2.0*(CON - V) - 5.0*PM + 4.0*XI + PL + 180.0                    
      GO TO 238                                                           
  121 E = 3.0*(CON - V) - 4.0*(PM - XI) - 90.0                            
      GO TO 258                                                           
  122 E = 2.0*(CON + TML) - VPP                                           
      GO TO 206                                                           
  123 E = CON + 2.0*(PM - XI - VP) + V + 270.0                            
  323 F = ((SIN(AW)*COS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2)  
     1)**1*(1.0/SQRT(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) +  
     2 0.1006))**2                                                        
      GO TO 888                                                           
  124 E = CON - 2.0*V - 4.0*(PM - XI) + VP - 270.0                        
  324 F = ((SIN(AW)*COS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2)  
     1)**2*(1.0/SQRT(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) +  
     2 0.1006))                                                           
      GO TO 888                                                           
  125 E = 6.0*(CON - PM) + 4.0*(XI - V - TML) + 2.0*PL - VPP              
  325 GO TO 254                                                           
  126 E = 6.0*CON - 5.0*PM + 4.0*(XI - V - TML) + PL - VPP                
  326 GO TO 254                                                           
  127 E = 6.0*CON - 4.0*TML - 3.0*PM + 2.0*(XI - V - VPP) + PL            
  327 GO TO 253                                                           
  128 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) - 2.0*TML + PL - VPP            
  328 GO TO 254                                                           
  129 E = 4.0*CON - 3.0*PM + 2.0*(XI - V) + PL - VPP                      
  329 GO TO 259                                                           
  130 E = 8.0*CON - 9.0*PM + 6.0*(XI - V) + 3.0*PL - 2.0*TML - VPP        
  330 GO TO 298                                                           
  131 E = 8.0*(CON - PM) + 6.0*(XI - V) + 2.0*(PL - TML) - VPP            
  331 GO TO 298                                                           
  132 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) - 2.0*TML + PL - VPP            
  332 GO TO 298                                                           
  133 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) + PL - VPP                      
  333 GO TO 254                                                           
  134 E = 4.0*CON - 3.0*PM + 2.0*(XI - V + TML) + PL - VPP                
  334 GO TO 259                                                           
  135 E = 10.0*CON - 9.0*PM + 8.0*(XI - V) - 2.0*TML + PL - VPP           
  335 F = ((COS(.5*AW)**4*COS(.5*AI)**4)/(COS(.5*UI)**4))**4*(1./SQRT(19  
     1.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.*U) + .0981))            
      GO TO 888                                                           
  136 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) + PL - VPP                      
  336 GO TO 298                                                           
  137 E = 6.0*CON - 5.0*PM + 4.0*(XI - V) - 2.0*TML + PL - VPP            
  337 GO TO 254                                                           
  138 E = 8.0*CON - 7.0*PM + 6.0*(XI - V) + 2.0*TML + PL - VPP            
  338 GO TO 298                                                           
  139 E = 10.0*CON - 9.0*PM + 8.0*(XI - V) + 2.0*TML + PL - VPP           
  339 GO TO 335                                                           
  140 E = 4.0*(CON + XI - V) - 6.0*PM + 2.0*PL + SL - TML - 270.0         
  340 GO TO 221                                                           
  141 E = 3.0*SL - 2.0*VP - TML - 90.0
  341 F = 1.0/(.8965*SIN(2.0*UI)**2 + .6001*SIN(2.0*UI)*COS(U) + .1006)
      GO TO 888
  142 E = TML + VP - 3.0*SL + 90.0
  342 GO TO 205
  143 E = 2.0*TML - VP
  343 GO TO 205
  144 E = 2.0*(TML - SL) + VPP
  344 GO TO 206
  145 E = 2.0*(TML - VPP) + 4.0*SL
  345 F = (1.0/SQRT(19.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.0*U) +
     10.0981))**2
      GO TO 888
  146 E = 2.0*TML - 2.0*(SL - PS)
  346 GO TO 203
  147 E = 4.0*TML + PS - SL
  347 GO TO 203
  148 E = 4.0*TML - 2.0*SL + VPP
  348 GO TO 206
  149 E = 4.0*TML + 6.0*SL - 3.0*VPP
  349 F = (1.0/SQRT(19.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.0*U) +
     10.0981))**3
      GO TO 888
  150 E = 4.0*TML - 3.0*(SL - PS)
  350 GO TO 203
  151 E = CON - V - 2.0*(PM - XI) + TML + 90.0
  351 GO TO 218
  152 E = 2.0*(PM - XI) + V + TML - CON - 90.0
  352 GO TO 218
  153 E = PM - XI - 90.0
  353 F = 0.31920/(SIN(UI) - SIN(UI)**3)
      GO TO 888
  154 E = PM - SL
  354 GO TO 235
  155 E = SL - PS
  355 GO TO 203
  156 E = 3.0*SL
  356 F = 1.0
      GO TO 888
  157 E = 4.0*SL
  357 GO TO 356
  158 E = 6.0*SL
  358 GO TO 356
  159 E = 8.0*SL
  359 GO TO 356
  160 E = 10.0*SL
  360 GO TO 356
  161 E = 12.0*SL
  361 GO TO 356
  162 E = 24.0*SL
  362 GO TO 356
  163 E = 3.0*TML + 180.0
  363 GO TO 356
  164 E = 5.0*TML + 180.0
  364 GO TO 356
  165 E = 2.0*(CON - V) - 4.0*(PM - XI) - 180.0
  365 F =((SIN(AW)*COS(.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))
     1**2
      GO TO 888
  166 E = TML + CON - VP + 270.0
  366 GO TO 205
  167 E = 3.0*(CON - PM) + 2.0*(XI - V) + PL - VP + 90.0
  367 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/COS(0.5*UI)**4 )*(1.0/SQRT(0.
     18965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) + 0.1006))
      GO TO 888
  168 E = 3.0*TML + 270.0 - SL
  368 GO TO 203
  169 E = 3.0*CON - VPP - VP + 90.0
  369 F = 1.0/(SQRT(19.0444*SIN(UI)**4 + 2.7702*SIN(UI)**2*COS(2.0*U) +
     10.0981)*SQRT(0.8965*SIN(2.0*UI)**2 + 0.6001*SIN(2.0*UI)*COS(U) +
     20.1006))
      GO TO 888
  170 E = 4.0*(CON - V) + 6.0*XI - 7.0*PM + PL - 180.0
  370 F = ((COS(0.5*AW)**4*COS(0.5*AI)**4)/(COS(0.5*UI)**4))*((SIN(AW)*C
     1OS(0.5*AW)**2*COS(.5*AI)**4)/(SIN(UI)*COS(.5*UI)**2))**2
      GO TO 888
  171 E = 4.0*(CON - V) - 6.0*(PM - XI) - 180.0
  371 GO TO 370
  172 E = 2.0*(TML + CON - V) - 4.0*(PM - XI) - 180.0
  372 GO TO 365
  173 E = 7.0*TML + 180.0
  373 GO TO 356
  174 E = 8.0*TML
  374 GO TO 356
  175 E = 10.0*TML
  375 GO TO 356
  488 PRINT 500, SPEED
      CALL PEXIT                                                           
  888 GO TO (400,390),ITYPE                                                
  390 E = E + FLOAT(MS(IPICK))*GONL - SPD(IPICK)*(TM/15.0)
  400 CALL TWOPI(E,1)                                                      
      F = 1.0/F                                                            
      RETURN                                                               
      END                                                                  
      SUBROUTINE YTEST(IYEAR,NODAYS,NOYR)
      DIMENSION NODAYS(12)
      NOYR = 365
      NODAYS(2) = 28
      FX = FLOAT(IYEAR)
      FY = AINT(FX*0.01)*100.0
      DIFF = FX - FY
      DIF = DIFF/4.0 + 0.0001
      IF(DIFF.NE.0.0) GO TO 9
      DIF = FY/400.0 + 0.0001
    9 FID = DIF - AINT(DIF)
      IF(FID.GT.0.001) GO TO 10
      NOYR = 366
      NODAYS(2) = 29
   10 RETURN
      END
      SUBROUTINE prcmp(ndata,x,y,pcd,rxy,RATIO)
C principle component computation
      parameter (nmax=100000)
      dimension x(nmax),y(nmax)
      dimension e(2,2),a(2,nmax)
      character*80 filein
C statement function
      zeta (angle,xi,yi) = xi * cos (angle) + yi * sin(angle)
      eta (angle,xi,yi) = -xi * sin (angle) + yi * cos(angle)
      call cmsv(x,ndata,xmean,xmnu,xmnl,dum1,dum2,dum3)
c     write(6,801) xmean,xmnl,xmnu
801   format(/' U-compoment mean velocity & 95% confidence interval',
     #3f8.4)
      call cmsv(y,ndata,ymean,ymnu,ymnl,dum1,dum2,dum3)
c     write(6,802) ymean,ymnl,ymnu
802   format(' V-component mean velocity & 95% confidence interval',
     #3f8.4)
      call atan3(xmean,ymean,ampl,direc)
      write(6,803) ampl,direc
803   format(/' Mean velocity --- Speed ',f8.4,' Direction ',f6.2)
      if(xmnl*xmnu.lt.0.and.ymnl*ymnu.lt.0)write(6,*)'Mean velocity ',
     #'is not significantly different than zero'
        xbar = 0.0
        ybar = 0.0
        do 200 n = 1,ndata
          xbar = xbar + x(n)
          ybar = ybar + y(n)
  200   continue
        xbar = xbar / ndata
        ybar = ybar / ndata
c       write (12,*) ' xbar = ',xbar,' ybar = ',ybar
        do 300 n = 1, ndata
          x(n) = x(n) - xbar
          y(n) = y(n) - ybar
  300   continue
        sxx = 0.0
        syy = 0.0
        sxy = 0.0
        do 400 n = 1,ndata
          sxx = sxx + x(n)*x(n)
          syy = syy + y(n)*y(n)
          sxy = sxy + x(n)*y(n)
  400   continue
        sxx = sxx / (ndata-1.)
        syy = syy / (ndata-1.)
        sxy = sxy / (ndata-1.)
        rxy = sxy / (sqrt(sxx) * sqrt(syy))
	write(6,*)
	write(6,*) 'Correlation coefficient = ',rxy
C principle angles
        theta_prime = 0.5 * atan2(2.*sxy,sxx-syy)
        pi = 4.*atan(1.)
        theta_p = 180. * theta_prime / pi
        write(6,*) 'Principal current direction = ',90.-Theta_p
     &  ,' degrees (true)'
C principle variances
        s11 = 0.5*(sxx+syy+sqrt((sxx-syy)**2+4.*sxy**2))
        s22 = 0.5*(sxx+syy-sqrt((sxx-syy)**2+4.*sxy**2))
        s12 = 0.5*(syy-sxx)*sin(2.*theta_prime)+sxy*cos(2.*theta_prime)
c       write (12,*) ' Principle variances:'
c       write (12,*) ' s11 = ',s11,' s22 = ',s22,' s12 = ',s12
	write(6,520) s11,s11 * 100. / (s11+s22)
	write(6,521) s22,s22 * 100. / (s11+s22)
        RATIO=s22/S11
        WRITE(*,*)'RATIO=',RATIO
520     format(/' Major axis variance = ',f10.4,5x,f8.4,' %')
521     format(' Minor axis variance = ',f10.4,5x,f8.4,' %')
        write(6,525) sqrt(s11)
        write(6,526) sqrt(s22)
525     format(/' Major axis standard deviation = ',f10.4)
526     format(' Minor axis standard deviation = ',f10.4)
C principle directions
        e(1,1) = cos(theta_prime)
        e(1,2) = sin(theta_prime)
        e(2,1) = -sin(theta_prime)
        e(2,2) = cos(theta_prime)
c       write (12,*) ' Principle directions: '
c       write (12,*) ' e(1) = (',e(1,1),',',e(1,2),')'
c       write (12,*) ' e(2) = (',e(2,1),',',e(2,2),')'
C Principle components
        do 500 n =1, ndata
          xr = x(n)
          yr = y(n)
          a(1,n) = zeta (theta_prime,xr,yr)
          a(2,n) = eta (theta_prime,xr,yr)
  500   continue
C Original data vs basis space
        rms_error = 0.0
c       write (12,*) ' Data vs. basis representation: '
        do 600 n = 1, ndata
          xbasis = a(1,n)*cos(theta_prime)-a(2,n)*sin(theta_prime)
          ybasis = a(1,n)*sin(theta_prime)+a(2,n)*cos(theta_prime)
c         write (12,*) ' d: (',x(n),',',y(n),')',
c    &                 ' b: (',xbasis,',',ybasis,')'
        rms_error = rms_error+(x(n)-xbasis)**2+(y(n)-ybasis)**2
  600   continue
        rms_error = sqrt (rms_error/ndata)
c       write (12,*) ' rms error = ',rms_error
!      write(12,*) 'Number of data points = ',ndata
!      write(12,801) xmean,xmnl,xmnu
!      write(12,802) ymean,ymnl,ymnu
!      write(12,803) ampl,direc
!      if(xmnl*xmnu.lt.0.and.ymnl*ymnu.lt.0)write(12,*)'Mean velocity ',
!     #'is not significantly different than zero'
!	write(12,*)
!	write(12,*) 'Correlation coefficient = ',rxy
!        write(12,*) 'Principal current direction = ',90.-Theta_p
!     &  ,' degrees (true)'
!	write(12,520) s11,s11 * 100. / (s11+s22)
!	write(12,521) s22,s22 * 100. / (s11+s22)
!        write(12,525) sqrt(s11)
!        write(12,526) sqrt(s22)
        pcd=90.-Theta_p
        RETURN
      end
