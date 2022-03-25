CCCCCC    This program is modified from pred.f of Chris Zervas so that it can make prediction 
CCCCCC    of crossing year, also call equarg.f to calculate XODE and VPU, instead of reading 
CCCCCC    from data file 'yr'
CCCCCC    modified on Feb. 13, 2005 by Aijun Zhang 
CCCCCC    lf95 pred.f -o pred.x
C         run at unix command line
C    pred_ngofs.x "$BEGINDATE" "$ENDDATE" $KINDAT $DELT $filein $fileout
C         BEGINDATE="2005 01 01 12 30"
C         ENDDATE=  "2005 12 31 12 30"
C         KINDAT=1, for current prediction; =2 for water level prediction
C         DELT is time interval of output time series in minutes 
C         XMAJOR is principle current direction in degrees
C         filein is input file name which includes tide constituents
C         fileout is output name which contains predicted water level or current time series
C                              wl
C    2003 01 02 00 00 00     0.5085
C    2003 01 02 00 06 00     0.5169
C                            sp        dir  
C   2003 01 02 00 00 00     0.7091   258.0500
C   2003 01 02 00 06 00     0.7237   258.3601
	   
C    1    2    1.   0.    0              ! nsta ipredk conv tconv il2
C   tss.out                                 ! Output time series file
C   0  4  15 0  6 30 0 0.1 1998 1998 106.0   IEL,IMMS,IDDS,TIME,IMME,IDDE,TIMEL,DELT,IYRS,IYRE,XMAJOR
C   Harmonic Analysis of Data in  325j4b05.dat                                      
C   29-Day H.A.  Beginning  4-15-1998  at Hour 17.30  along 106 degrees             
C   12718
C       1931621828140022115186641641117973376 68733224 81743495 5868 163
C       2    0   0  804  92    0   0 36211666  4431590    0   0 24821455
C       3  3513257  6521961    0   0  5803436  6463317    0   0    0   0
C       4    0   0    0   0    0   0  3113546 15863554  8262104  1122127
C       5  213  13 39053385    0   0    0   0 26691641    0   0 38082138
C       6 28911704    0   0
C   Harmonic Analysis of Data in  325j4b05.dat                                      
C   29-Day H.A.  Beginning  4-15-1998  at Hour 17.30  along 196 degrees             
C -5820
C       1 57222694 14073452 22192976 1203 154 47261638 1292 478 26243445
C       2    0   0  658 796    0   0  4302938  6232349    0   0  2953258
C       3   563430   403046    0   0   92 316  1023593    0   0    0   0
C       4    0   0    0   0    0   0   49 617  251 639   833422   113482
C       5   34 800  398 178    0   0    0   0  3172976    0   0  3833513
C       6 12661394    0   0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C****          Purpose:
C                 1. Calculates residuals
C                 2. Calculates a predicted series
C **************************************************************
C
C               UNIT=10 is the output file written in ASCII format.
C               UNIT=11 is the input observations in CDF or ASCII format.
C
C
C
C      UNIT=         PATH NAME=                    INPUT      OUTPUT
C
C       5       Redirected std. input (< pathname)   X
C       6                                                       X
C      10                                                       X
C      11                                            X
C************************************************************************************
C      NSTA =     NUMBER OF STATIONS TO PREDICT
C      CONV =     FACTOR FOR CONVERTING PREDICTED TIME SERIES TO NEW UNITS
C**** Conversion options for time ***
C
C     TCONV = 0   NO CONVERSION OF PREDICTED TIMES
C     TCONV = TIME MERIDIAN FOR WHICH THE KAPPA PRIMES IN THE HARMONIC
C             CONSTANTS WERE DERIVED. THIS OPTION IS USED TO CONVERT
C             THE TIMES TO GREENWICH IF THE CONSTANTS WERE CALCULATED
C             FOR A LOCAL TIME MERIDIAN. A REASON FOR USING THIS OPTION
C             IS IF COMPARISONS WITH ACTUAL DATA IS REQUIRED WHICH IS
C             IN GREENWICH TIME AND THE HARMONIC CONSTANTS WERE OBTAINED
C             FROM THE PREDICTION BRANCH WHERE HARMONIC CONSTANTS ARE
C             FOR LOCAL MERIDIANS ALWAYS
C
C     TCONV changed to hours to shift predicted time series 
C     (positive = later) --- Chris Zervas (7/97)
C
c**** conversion option for using 2mn2 in the harmonic constants  
c     file versus the standard L2
c
c     ***note*** it is node (33) which is re-calculated
c
c     IL2 = 0 --- use the standard <L2> harmonic constants
c           1 --- use the <2MN2> harmonic constants
c
C
C*****************************************************************
C
C     **NOTE(1)** USE FORMAT NO. 531 FOR HARMONIC CONSTANTS
C        TEMPORARY FORMULATION FORMAT 532 IS USED
C
C******************************************************************
C
C
C     IPREDK = 0 -- USED TO CALCULATE THE DIFFERENCES BETWEEN A SET OF
C                   PREDICTED AND OBSERVED SERIES. THIS CALCULATION IS
C                   DONE USING AN INPUT OF OBSERVED VALUES IN ASCII
C                   FORMAT AND OUTPUT IS WRITTEN IN AN ASCII FILE
C
C     IPREDK = 1 -- USED TO CALCULATE THE DIFFERENCES BETWEEN A SET OF
C                   PREDICTED AND OBSERVED SERIES. THIS CALCULATION IS
C                   DONE USING AN INPUT OF OBSERVED VALUES IN CDF
C                   FORMAT AND OUTPUT IS WRITTEN IN AN ASCII FILE
C
C     **NOTE(1)** STARTING TIME FOR PREDICTIONS IS SET EQUAL TO
C                 THE STARTING TIME DEFINED BY THE VALUE OF ISTART
C     **NOTE(2)** NOS1 IS CALCULATED BY THE TIMES OF THE OBSERVED
C                 SERIES
C
C
C      IPRED = 2 -- USED TO CALCULATE A PREDICTED SERIES
C
C     **NOTE(3)** THE PREDICTIONS AND DIFFERENCES CAN NOT BE PERFORMED
C                 ACROSS 2 YEARS
C     **NOTE(4)** THE YEARS FOR EACH NSTA MUST BE THE SAME
C     **NOTE(5)** USE FORMAT NO. 532 FOR HARMONIC CONSTANTS
C ***************************** modified IEL to match KINDAT
C     IEL =    THE ELEMENT OF THE DATA SERIES TO PERFORM CALCULATIONS
C              1   MAJOR/MINOR COMPONENTS OF VECTOR VARIABLE (I.E. CURRENT)
C              2   SCALAR VARIABLE (I.E. TIDAL HEIGHT)
C              3   TEMPERATURE (CDF INPUT FIELD)
C              4   CONDUCTIVITY (CDF INPUT FIELD)
C              6   PRESSURE (CDF INPUT FIELD)
C     IYRS=YEAR OF THE FIRST DATA POINT(CAN NOT GO ACROSS YRS)
C     IMMS=MONTH OF FIRST DATA POINT
C     IDDS=DAY OF FIRST DATA POINT
C     TIME=TIME OF FIRST DATA POINT
C     MON=MONTH OF LAST DATA POINT
C     IDDE=DAY OF LAST DATA POINT
C     TIMEL=TIME OF LAST DATA POINT
C     DELT= desired time interval in hours, delt=0.1 for 6 minutes data
C     NOS1=NUMBER OF DATA POINTS PER HOUR
C     XMAJOR ---- AXIS FOR MAJOR/MINOR COMPONENTS (MAKE SURE HARMONIC
C              CONSTANTS WERE DERIVED ALONG THIS AXIS)
C                       IF A 0.0 IS READ --- 0 DEGRESS TRUE IS
C                       ASSUMED. ALWAYS READ IN THE CONSTITUENTS
C                       FOR THE MAJOR AXIS FIRST.
C
C Modified by Lianyuan Zheng on 03/01/2017

      PROGRAM PREDICTION
      PARAMETER (MXDIM=200000)
      PARAMETER (XMISS=999.0)
      character*100 cdfin,cdfout,BUFFER
      CHARACTER*10 ALIST
      CHARACTER*100 HEAD(2)
      CHARACTER*10 LABLE(180)
      DIMENSION TSTART(4)
      DIMENSION XDATA(14,50),ALIST(37),IHEAD(34),DAT(6),SUM(6),SMN(6),
     1 TIM(MXDIM),SPEED(MXDIM),DIREC(MXDIM),STORX(MXDIM),STORN(MXDIM)
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE,JULIAN
      REAL*8 YEARB,MONTHB,DAYB,HOURB
      REAL*8 DAYJ,DELT,TIMEE
      REAL NOS1
      REAL*8 SPD(180),FFF(180),VAU(180)
      INTEGER ORDER(37),IAMP(37),IEPOC(37)
      LOGICAL FEXIST

      common/virt1/
     1      A(37),AMP(37),EPOC(37),XODE(114),VPU(114),
     2      XCOS(1025),ARG(37),TABHR(24),ANG(37),SPD0(37),
     3      EPOCH(37),AMPA(37),IYR(15),NUM(15),ISTA(6),NO(6),
     4      JODAYS(12),C(37)
      common /XDATA/XDATA
      common /speeds/spd
      common /names/lable

      DATA order/1,3,2,5,21,4,22,28,24,30,13,25,12,8,16,11,27,15,14,
     1  35,37,36,34,33,20,18,10,9,19,17,32,26,7,29,6,23,31/
      DATA (ALIST(I),I=1,37) /'M(2)      ','S(2)      ','N(2)      ',
     1                        'K(1)      ','M(4)      ','O(1)      ',
     2                        'M(6)      ','MK(3)     ','S(4)      ',
     3                        'MN(4)     ','NU(2)     ','S(6)      ',
     4                        'MU(2)     ','2N(2)     ','OO(1)     ',
     5                        'LAMBDA(2) ','S(1)      ','M(1)      ',
     6                        'J(1)      ','MM        ','SSA       ',
     7                        'SA        ','MSF       ','MF        ',
     8                        'RHO(1)    ','Q(1)      ','T(2)      ',
     9                        'R(2)      ','2Q(1)     ','P(1)      ',
     1                        '2SM(2)    ','M(3)      ','L(2)      ',
     2                        '2MK(3)    ','K(2)      ','M(8)      ',
     3                        'MS(4)     '/

      LIN  = 5
      LOUT = 6

      NSTA   = 1
      IPREDK = 2
      CONV   = 1.0
      TCONV  = 0.0
      IL2    = 0

      IF(IPREDK .GT. 2 .OR. IPREDK .LT. 0) THEN
        WRITE(*,*) 'Error in IPREDK value'
        STOP
      ENDIF

C  Develop Cosine Table
      H = 0.00153398078789
      DO 35 I = 1, 1024
        XCOS(I) = COS(H*(I-1))
35    CONTINUE
      XCOS(1025) = 0.0
      MS0 = 1
      CON = 1024.0 / 90.0

      DO 612 JOB = 1, NSTA
        CALL GETARG(1,BUFFER)
        READ(BUFFER,*) IYRS, IMMS, IDDS, IHHS, MNS
        TIME = IHHS + MNS/60.0
        CALL GETARG(2,BUFFER)
        READ(BUFFER,*) IYRE, IMME, IDDE, IHHE, MNE
        TIMEL = IHHE + MNE/60.0
        CALL GETARG(3,BUFFER)
        READ(BUFFER,*) IEL
        CALL GETARG(4,BUFFER)
        READ(BUFFER,*) DELT
        DELT = DELT/60.0  !! convert to hours from minutes 11/30/2006
        CALL GETARG(5,CDFIN)
        CALL GETARG(6,CDFOUT)

        YEARB  = dble(IYRS*1.0)
        MONTHB = dble(1.0)
        DAYB   = dble(1.0)
        HOURB  = dble(0.0)
        JBASE_DATE=JULIAN(YEARB,MONTHB,DAYB,HOURB)

        NOS1 = 1.0 / DELT   !!!  data points per hour
C        WRITE(*,*) 'Run pred.f from ',IYRS,IMMS,IDDS,IHHS, ' to ',
C     1   IYRE,IMME,IDDE,IHHE

        YEARB  = DBLE(IYRS*1.0)
        MONTHB = DBLE(IMMS*1.0)
        DAYB   = DBLE(IDDS*1.0)
        HOURB  = DBLE(IHHS*1.0 + MNS/60.0)
        JDAY0  = JULIAN(YEARB,MONTHB,DAYB,HOURB)

        YEARB  = DBLE(IYRE*1.0)
        MONTHB = DBLE(IMME*1.0)
        DAYB   = DBLE(IDDE*1.0)
        HOURB  = DBLE(IHHE*1.0 + MNE/60.0)
        JDAY1  = JULIAN(YEARB,MONTHB,DAYB,HOURB)

        TSTART(4) = IYRS*1.0
        TSTART(3) = (JDAY0-JBASE_DATE)*1.0
        TSTART(2) = IHHS*1.0
        TSTART(1) = MNS*1.0
        NPTS = INT((JDAY1-JDAY0)*24/DELT+1+0.1)

        INQUIRE(FILE = TRIM(CDFIN), EXIST = FEXIST)
        IF(.NOT. FEXIST) THEN
          WRITE(*,*)
     1    'Tide Constituent File of '//TRIM(CDFIN)//' does not exist'
          OPEN(10,FILE=TRIM(CDFOUT),STATUS='UNKNOWN',FORM='FORMATTED')
	  DO N = 1, NPTS
            JDAY = JDAY0 + dble((N-1)*DELT/24.0) 
            CALL GREGORIAN(jday,yearb,monthb,dayb,hourb)
            IYEAR = INT(yearb)
            ICM   = INT(monthb)
            ICD   = INT(dayb)
            IHR   = INT(hourb)
            IMN   = INT((hourb-IHR)*60)
            ISEC  = 0
            DAYJ  = JDAY - JBASE_DATE + dble(1.0)
            WRITE(10,571) DAYJ,IYEAR,ICM,ICD,IHR,IMN,-99.99
	  ENDDO
	  CLOSE(10)  
          STOP
        ENDIF

        OPEN(30,file = trim(CDFIN),status='old')
        IZAJ = 0
119     CIND = NOS1

111     READ(30,'(A100)') HEAD(1) 
        READ(30,'(A100)') HEAD(2)
        READ(30,'(I6)')  IDATUM
        DATUM = IDATUM / 1000.0

        DO N = 1, 5
          NN = 7*(N-1)
          READ(30,'(I8,7(I5,I4))') 
     1         NO(N),(IAMP(NN+JJ),IEPOC(NN+JJ),JJ = 1, 7)
        END DO
        READ(30,'(I8,2(I5,I4))') NO(6),IAMP(36),IEPOC(36),
     2       IAMP(37),IEPOC(37)
        DO N = 1, 37
          AMP(N) = IAMP(N) / 1000.0
          EPOC(N) = IEPOC(N) / 10.0
        END DO

C  READ IN PRINCIPLE DURRENT DIRECTION FROM HARMONIC CONSTANTS DATA FILE
        IF((IZAJ .EQ. 0) .and. (IEL .eq. 1)) THEN
          IND = INDEX(HEAD(2),'along') + 5
          READ(HEAD(2)(IND:IND+4),'(I4)') IPCD
          XMAJOR = IPCD
          IZAJ = 1
          WRITE(*,*) 'XMAJOR: ',XMAJOR
        ENDIF

        DO 114 L = 1, 6
          IF(NO(L) .NE. L) GO TO 450
114     CONTINUE

C     CONVERT CONSTANTS IF TCONV IS NOT EQUAL TO ZERO
C     TCONV IS THE TIME MERIDIAN
        IF(TCONV .EQ. 0.0) GO TO 120
        WRITE(*,*) '    '
        WRITE(*,'(A100)') HEAD(1)
        WRITE(*,'(A100)') HEAD(2)
        WRITE(*,*) '    '
        WRITE(LOUT,994) TCONV
994     FORMAT('Values of the Epochs before ',F6.2,' hour time shift')

        DO 650 J = 1, 37
          IF(AMP(J) .NE. 0.0) THEN
            EPOC(J) = EPOC(J) + A(J) * TCONV
            EPOC(J) = MOD(EPOC(J), 360.0)
          END IF
650     CONTINUE
        WRITE(*,*) 'AMP: ',AMP
        WRITE(*,*) 'EPOC: ',EPOC

120     IF(MS0 .EQ. 2) GOTO 150

CCCCCCCC  using equarg replace yrcrds.f
C The start time (IYRS,IMMS,IDDS) must be the same as the base date since  
C the initial phase is computed in equarg rechoned from IYRS,IMMS,IDDS 
        LENGTH = INT(JDAY1-JDAY0) + 1
        CALL EQUARG(37,IYRS,1,1,365,LABLE(1),FFF(1),VAU(1))
        DO J = 1, 37
          VPU(J) = VAU(ORDER(J))
          XODE(J) = FFF(ORDER(J))
        END DO

150     DO 155 J = 1, 37
          C(J) = A(J) * (CON/CIND)
155     CONTINUE

C     SET UP TABLES FOR NON-ZERO CONSTITUENTS
        K = 0
        DO 180 J = 1, 37
          IF(AMP(J) .NE. 0.0) THEN
            K = K + 1
            AMPA(K) = AMP(J) * XODE(J)
            TEMX = VPU(J) - EPOC(J)
            IF(TEMX .LT. 0.0) THEN
              TEMX = TEMX + 360.0
            END IF
            EPOCH(K) = TEMX * CON
            SPD0(K) = C(J)
          END IF
180     CONTINUE
        NOCON = K

C     OPERATING TABLES NOW STORED AS AMPA(K),EPOCH(K),SPD0(K)
C**** CHECK LENGTH OF SERIES

        IF(NPTS .LE. MXDIM) GO TO 191
        NPTS = MXDIM
        PRINT 100
        PRINT 1190, NPTS
191     WRITE (LOUT,1002) NPTS
1002    FORMAT(' Total number of prediction times = ', I10)

C**** DETERMINE FIRST HOUR OF TIME PERIOD  at 00:00 of Jan. 1, first=0.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC modified by zaj on May 13, 2004
        FIRST = (JDAY0 - JBASE_DATE) * 24.0 * CIND
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc      
        DO 220 J = 1, MXDIM
          STORX(J) = 0.0
220     CONTINUE
        KOUNT = 0
        KT = 0

        DO 380 K = 1, NPTS
C**** THE PREDICTION/TIME STEP IS STORED IN VARIABLE PRED
          PRED = 0.0
C
          IF(KOUNT.GT.0) GO TO 260
          KOUNT = 1
231       IF(NOCON.EQ.0) GO TO 375

          DO 250 J = 1, NOCON
            ARGU = SPD0(J) * FIRST + EPOCH(J)
            ARG(J) = AMOD(ARGU,4096.0)
250       CONTINUE
          GO TO 290
260       IF(NOCON .EQ. 0) GO TO 375

          DO 280 J = 1, NOCON
            ARG(J) = ARG(J) + SPD0(J)
            ARG(J) = MOD(ARG(J),4096.0)
280       CONTINUE

290       DO 374 J = 1, NOCON
            IF(ARG(J) - 1024.0) 320, 320, 300
300         IF(ARG(J) - 2048.0) 350, 350, 310
310         IF(ARG(J) - 3072.0) 360, 360, 330
320         ANG(J) = ARG(J)
            GO TO 340
330         ANG(J) = 4096.0 - ARG(J)
340         NP = ANG(J) + 1.5
            PRED = PRED + AMPA(J) * XCOS(NP)
            GO TO 374
350         ANG(J) = 2048.0 - ARG(J)
            GO TO 370
360         ANG(J) = ARG(J) - 2048.0
370         NP = ANG(J) + 1.5
            PRED = PRED - AMPA(J) * XCOS(NP)
374       CONTINUE

375       IF(K .NE. NPTS) GO TO 381
          IF(IPREDK .EQ. 0) GO TO 381
          IF(KT .EQ. 1) GO TO 378
          FIRST = FIRST + NPTS - 1.0
          KT = 1
          CHECK = PRED
          PRED = 0.0
          GO TO 231
378       CKSUM = CHECK - PRED

C**** CONVERT RESULTS ACCORDING TO CONV
381       PRED = (PRED + DATUM) * CONV

C**** STORE THE RESULTS
621       STORX(K) = PRED
380     CONTINUE

        NOS2 = 0
        IF(IEL .EQ. 1) NOS2 = 1
        IF(IEL .EQ. 1 .AND. XMAJOR .NE. 0) NOS2 = 2
        IF(NOS2 .EQ. 0) GO TO 5004
        IF(ms0 .EQ. 2) GO TO 4996
        IF(NOS2 .EQ. 1) GO TO 4995
        PRINT 5006
        GO TO 5004
4995    PRINT 5008
        GO TO 5004
4996    IF(NOS2.EQ.1) GO TO 5002
        PRINT 5011
        GO TO 5004
5002    PRINT 5009
5004    PRINT 5902

        IF(NOS2) 400, 400, 395
395     IF(MS0 .EQ. 2) GO TO 410
        MS0 = 2

C**** STORES THE MAJOR AXIS
400     DO 1395 I = 1, NPTS
          STORN(I) = STORX(I)
1395    CONTINUE

        IF(NOS2) 410,410,119
410     CONTINUE

C**** FORM NEW DATA ARRAY
        OPEN(10,FILE=TRIM(CDFOUT),STATUS='UNKNOWN',FORM='FORMATTED')
        IDET = 0
        IPRED = 1
        JDAY = JDAY0
        TIMEE = 0.0
      
427     CONTINUE
        JDAY = JDAY0+TIMEE
        CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
        IYEAR = INT(YEARB)
        ICM   = INT(MONTHB+0.001)
        ICD   = INT(DAYB+0.001)
        IHR   = INT(HOURB+0.001)
        IMN   = INT((HOURB-IHR)*60+0.1)
        DAYJ  = JDAY - JBASE_DATE + dble(1.0)

        IF(IEL .EQ. 1) GO TO 418
        WRITE(10,571) DAYJ,IYEAR,ICM,ICD,IHR,IMN,STORN(IPRED)
        GO TO 421

418     CONTINUE 
        CALL VELDIR(STORN(IPRED),STORX(IPRED),DR,SP)
        DR = DR + XMAJOR
        IF(DR .GT. 360.0) DR = DR - 360.0
        u = sp * sin(dr*3.1415926/180.0)
        v = sp * cos(dr*3.1415926/180.0)
        WRITE(10,571) DAYJ,IYEAR,ICM,ICD,IHR,IMN,SP,DR,U,V
571     FORMAT(F13.8,I5,4I3,4F10.4)

421     IPRED = IPRED + 1
        IF(IPRED .LT. NPTS) GO TO 424
        IF(IPRED .EQ. NPTS) GO TO 424
        GO TO 602

424     CONTINUE      
        TIMEE = TIMEE + dble(DELT/24.0)
        GO TO 427

602     MS0 = 1
612   CONTINUE
      CLOSE(30)
      CLOSE(10)
      STOP

450   WRITE(*,*) 'STATION CARDS OUT OF ORDER'
      CLOSE(30)
      CLOSE(10)
      STOP

100   FORMAT(1H0)
5902  FORMAT(1X,4H    )
5006  FORMAT(1X,38HHarmonic Constants (Major Axis) ------)
5011  FORMAT(1X,38HHarmonic Constants (Minor Axis) ------)
5008  FORMAT(1X,26HHarmonic Constants (North))
5009  FORMAT(1X,25HHarmonic Constants (East))
1190  FORMAT(' LENGTH OF PREDICTIONS SHORTENED TO',I10)
      END


      BLOCK DATA
      PARAMETER (MXDIM=200000)

      CHARACTER*10 LABLE(180)
      DOUBLE PRECISION SPD(180)
      INTEGER MS(180)

      common /speeds/spd
      common /names/lable
      common /mmss/ms
      common /datec/idtbc(12,2),iltbc(12,2)
      common /datej/idtbj(12),iltbj(12)
      common/virt1/
     1      A(37),AMP(37),EPOC(37),XODE(114),VPU(114),
     2      XCOS(1025),ARG(37),TABHR(24),ANG(37),SPD0(37),
     3      EPOCH(37),AMPA(37),IYR(15),NUM(15),ISTA(6),NO(6),
     4      JODAYS(12),C(37)

      DATA ((IDTBC(I,J),J=1,2),I=1,12)/1,31,32,59,60,90,91,120,121,151,
     1                                  152,181,182,212,213,243,244,273,
     2                                  274,304,305,334,335,365/
      DATA  ((ILTBC(I,J),J=1,2),I=1,12)/1,31,32,60,61,91,92,121,122,
     1                              152,153,182,183,213,214,244,245,
     2                                    274,275,305,306,335,336,366/
      DATA (IDTBJ(I),I=1,12)/1,32,60,91,121,152,182,213,244,274,305,335/
      DATA (ILTBJ(I),I=1,12)/1,32,61,92,122,153,183,214,245,275,306,336/
      DATA(TABHR(I), I=1,24)/  -24.,  720., 1392., 2136., 2856., 3600.,
     1    4320., 5064., 5808., 6528., 7272., 7992.,  -24.,  720., 1416.,
     2    2160., 2880., 3624., 4344., 5088., 5832., 6552., 7296., 8016./
      DATA (A(I), I=1,37)/        28.9841042,  30.0000000,  28.4397295,
     1   15.0410686,  57.9682084,  13.9430356,  86.9523127,  44.0251729,
     2   60.0000000,  57.4238337,  28.5125831,  90.0000000,  27.9682084,
     3   27.8953548,  16.1391017,  29.4556253,  15.0000000,  14.4966939,
     4   15.5854433,   0.5443747,   0.0821373,   0.0410686,   1.0158958,
     5    1.0980331,  13.4715145,  13.3986609,  29.9589333,  30.0410667,
     6   12.8542862,  14.9589314,  31.0158958,  43.4761563,  29.5284789,
     7   42.9271398,  30.0821373, 115.9364169,  58.9841042/
      DATA (JODAYS(LL),LL=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA(spd(i),i=1,37)/28.9841042d0,   28.4397295d0,   30.0000000d0,
     1    13.9430356d0,   15.0410686d0,   30.0821373d0,   29.5284789d0,
     2    27.8953548d0,   30.0410667d0,   29.9589333d0,   29.4556253d0,
     3    27.9682084d0,   28.5125831d0,   15.5854433d0,   14.4966939d0,
     4    16.1391017d0,   14.9589314d0,   13.3986609d0,   12.8542862d0,
     5    13.4715145d0,   57.9682084d0,   86.9523127d0,  115.9364169d0,
     6    60.0000000d0,   90.0000000d0,   43.4761563d0,   15.0000000d0,
     7    44.0251729d0,   42.9271398d0,   57.4238337d0,   58.9841042d0,
     8    31.0158958d0,   01.0980331d0,   01.0158958d0,   00.5443747d0,
     9     0.0410686d0,    0.0821373d0/
      DATA(spd(i),i= 38,48)/12.9271398d0, 14.0251729d0,   14.5695476d0,
     1    15.9748272d0,   16.0569644d0,   30.5443747d0,   27.4238337d0,
     2    28.9019669d0,   29.0662415d0,   26.8794590d0,   26.9523126d0/
      DATA(spd(i),i= 49,87)/27.4966873d0, 31.0980331d0,   27.8039338d0,
     1    28.5947204d0,   29.1483788d0,   29.3734880d0,   30.7086493d0,
     2    43.9430356d0,   45.0410686d0,   42.3827651d0,   59.0662415d0,
     3    58.4397295d0,   57.4966873d0,   56.9523127d0,   58.5125831d0,
     4    56.8794590d0,   59.5284789d0,   71.3668693d0,   71.9112440d0,
     5    73.0092770d0,   74.0251728d0,   74.1073100d0,   72.9271398d0,
     6    71.9933813d0,   72.4649023d0,   88.9841042d0,   86.4079380d0,
     7    87.4238337d0,   87.9682084d0,   85.3920421d0,   85.8635632d0,
     8    88.5125831d0,   87.4966873d0,   89.0662415d0,   85.9364168d0,
     9    86.4807916d0,   88.0503457d0,  100.3509735d0,  100.9046318d0/
      DATA(spd(i),i= 88,124)/101.9112440d0,103.0092771d0,116.4079380d0,
     1   116.9523127d0,  117.9682084d0,  114.8476674d0,  115.3920422d0,
     2   117.4966873d0,  115.4648958d0,  116.4807916d0,  117.0344500d0,
     3   118.0503457d0,  129.8887360d0,  130.4331108d0,  130.9774855d0,
     4   131.9933813d0,  144.3761464d0,  144.9205211d0,  145.3920422d0,
     5   145.9364169d0,  146.4807916d0,  146.9523127d0,  160.9774855d0,
     6   174.3761464d0,  174.9205211d0,  175.4648958d0,  175.9364169d0,
     7    14.9178647d0,   15.0821353d0,   15.1232059d0,  15.5125897d0,
     8    30.6265119d0,   27.3416965d0,   42.9271397d0,   60.0821373d0,
     9    16.1391016d0,   12.8450026d0/
      DATA(spd(i),i=125,140)/ 26.9615963d0, 27.5059710d0, 28.6040041d0,
     1    57.5059710d0,   58.5218668d0,   85.4013258d0,   85.9457005d0,
     2    86.4900752d0,   87.5059710d0,   88.5218668d0,  115.4741794d0,
     3   116.4900752d0,  117.5059710d0,  146.4900752d0,  175.4741794d0,
     4    41.9205276d0/
      DATA(spd(i),i=141,150)/ 15.1232058d0, 14.8767942d0, 30.0000001d0,
     1    29.9178627d0,   30.1642746d0,   29.9178666d0,   59.9589333d0,
     2    59.9178627d0,   60.2464119d0,   59.8767999d0/
      DATA(spd(i),i=151,164)/ 28.9430356d0, 01.0569644d0, 00.5490165d0,
     1          00.5079479d0, 00.0410667d0, 00.1232059d0, 00.1642746d0,
     2          00.2464118d0, 00.3285841d0, 00.4106864d0, 00.4928237d0,
     3          00.9856473d0, 45.0000000d0, 75.0000000d0/
      DATA(spd(i),i=165,175)/ 27.8860712d0, 30.0410686d0, 43.4807981d0,
     1          44.9589314d0, 45.1232059d0, 56.3258007d0, 56.8701754d0,
     2          57.8860712d0,105.0000000d0,120.0000000d0,150.0000000d0/

      DATA(ms(n),n=1,  37)/3*2,2*1,8*2,7*1,4,6,8,4,6,3,1,2*3,2*4,2,5*0/
      DATA(ms(n),n=38, 48)/ 5*1,6*2/
      DATA(ms(n),n=49, 96)/7*2,3*3,7*4,8*5,12*6,4*7,7*8/
      DATA(ms(n),n=97, 124)/3*8,4*9,6*10,11,4*12,4*1,2*2,3,4,2*1/
      DATA(ms(n),n=125,140)/3*2,2*4,5*6,3*8,10,12,3/
      DATA(ms(n),n=141,150)/2*1,4*2,4*4/
      DATA(ms(n),n=151,164)/2,11*0,3,5/
      DATA(ms(n),n=165,175)/2*2,3*3,3*4,7,8,10/

      DATA(lable(m),m = 1,37)/ 'M(2)      ','N(2)      ','S(2)      ',
     1            'O(1)      ','K(1)      ','K(2)      ','L(2)      ',
     2            '2N(2)     ','R(2)      ','T(2)      ','Lambda(2) ',
     3            'Mu(2)     ','Nu(2)     ','J(1)      ','M(1)      ',
     4            'OO(1)     ','P(1)      ','Q(1)      ','2Q(1)     ',
     5            'Rho(1)    ','M(4)      ','M(6)      ','M(8)      ',
     6            'S(4)      ','S(6)      ','M(3)      ','S(1)      ',
     7            'MK(3)     ','2MK(3)    ','MN(4)     ','MS(4)     ',
     8            '2SM(2)    ','Mf        ','Msf       ','Mm        ',
     9            'Sa        ','Ssa       '/
      DATA lable(38)/'SK(3)     '/
      DATA(lable(m),m =39,77)/ 'Sigma(1)  ','MP(1)     ','Chi(1)    ',
     1            '2PO(1)    ','SO(1)     ','2NS(2)    ','2NK2S(2)  ',
     2            'MNS(2)    ','MNK2S(2)  ','2MS2K(2)  ','2KN2S(2)  ',
     3            'OP(2)     ','MKS(2)    ','M2(KS)(2) ','2SN(MK)(2)',
     4            'MSN(2)    ','2KM(SN)(2)','SKM(2)    ','NO(3)     ',
     5            'SO(3)     ','N(4)      ','3MS(4)    ','MNKS(4)   ',
     6            'SN(4)     ','KN(4)     ','MK(4)     ','SL(4)     ',
     7            'MNO(5)    ','2MO(5)    ','3MP(5)    ','MNK(5)    ',
     8            '2MP(5)    ','2MK(5)    ','MSK(5)    ','3KM(5)    ',
     9            '3NKS(6)   ','2NM(6)    ','2NMKS(6)  ','2MN(6)    '/
      DATA(lable(m),m=78,114)/ '2MNKS(6)  ','MSN(6)    ','MKN(6)    ',
     1            '2MS(6)    ','2MK(6)    ','NSK(6)    ','2SM(6)    ',
     2            'MSK(6)    ','2MNO(7)   ','2NMK(7)   ','2MSO(7)   ',
     3            'MSKO(7)   ','2(MN)(8)  ','3MN(8)    ','3MNKS(8)  ',
     4            '2MSN(8)   ','2MNK(8)   ','3MS(8)    ','3MK(8)    ',
     5            'MSNK(8)   ','2(MS)(8)  ','2MSK(8)   ','2M2NK(9)  ',
     6            '3MNK(9)   ','4MK(9)    ','3MSK(9)   ','4MN(10)   ',
     7            'M(10)     ','3MNS(10)  ','4MS(10)   ','2MNSK(10) ',
     8            '3M2S(10)  ','4MSK(11)  ','4MNS(12)  ','5MS(12)   ',
     9            '3MNKS(12) ','4M2S(12)  '/
      DATA(lable(m),m=115,120)/'2NP(3)    ','ML(4)     ','2ML(6)    ',
     1            'MSL(6)    ','2MSL(8)   ','3ML(8)    '/
      DATA(lable(m),m=121,128)/'TK(1)     ','RP(1)     ','KP(1)     ',
     1            'THETA(1)  ','KJ(2)     ','OO(2)     ','MO(3)     ',
     2            'SK(4)     '/
      DATA(lable(m),m=129,138)/'MLN2S(2)  ','2ML2S(2)  ','MKL2S(2)  ',
     1            '2MLS(4)   ','2NMLS(6)  ','2MLNS(6)  ','3MLS(6)   ',
     2            '4MLS(8)   ','3MSL(10)  ','4MSL(12)  '/
      DATA(lable(m),m=139,140)/'2KO(1)    ','2OK(1)    '/
      DATA(lable(m),m=141,150)/'2KP(1)    ','2PK(1)    ','KP(2)     ',
     1            '2SK(2)    ','2KS(2)    ','2TS(2)    ','ST(4)     ',
     2            '3SK(4)    ','3KS(4)    ','3TS(4)    '/
      DATA(lable(m),m=151,164)/'SO(2)     ','SO(0)     ','.5Mf      ',
     1            '.5Msf     ','ST(0)     ','3Sa       ','4Sa       ',
     2            '6Sa       ','8Sa       ','10Sa      ','12Sa      ',
     3            '24Sa      ','S(3)      ','S(5)      '/
      DATA(lable(m),m=165,175)/'O(2)      ','SK(2)     ','NK(3)     ',
     1            'SP(3)     ','K(3)      ','NO(4)     ','MO(4)     ',
     2            'SO(4)     ','S(7)      ','S(8)      ','S(10)     '/
      END


      SUBROUTINE EQUARG(NSPED,IYR,ICM,ICD,LENGTH,LABEL,FFF,VAU)
c     program to calculate equilibrium arguments and node factors
c     by e.e.long      (slight revision by b.b.parker)
c       major revisions by len hickman 6/11/86 to make this a subroutine
c       of lsqha.  v is calculated for the beginning of the series.
c       u and f are adjusted to the midpoint of the series (regardless
c       of length).
c     More revisions by Geoff French 9/1/86

      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE,JULIAN,YEARB
      REAL*8 MONTHB,DAYB,HOURB
      REAL*8 SPEED,SPD,FFF,VAU
      CHARACTER*10 LNAME,LABEL
      DIMENSION LABEL(180),VAU(180),FFF(180)
      DIMENSION SPD(180),VUU(180),VOU(180)
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp

C     ******************************************************************
C     *   nsped = number of constituents to be calculated              *
C     *                                                                *
C     *     ICM  = month of first data point                        *
C     *     ICD    = day  of first data point                         *
C     *     IYR   = year of first data point                         *
C     *     length  = length of series to be analyzed (in days)        *
C     ******************************************************************
      YEARB  = IYR
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JBASE_DATE=JULIAN(YEARB,MONTHB,DAYB,HOURB)

      IYER = IYR
      DAYB = 0.0
      DAYM = 183.0
      GRBS = 0.0
      GRMS = 12.0

      IF(MOD(IYER, 4) .EQ. 0) THEN
        DAYM = 184.0
        GRMS = 0.0
      END IF
      XYER = IYER
      DAYBB = REAL(DAYB)

C     Compute basic astronimical constants for the time period wanted
      CALL ASTRO(XYER, DAYBB, DAYM, GRBS, GRMS)
      TM   = 0.0
      GONL = 0.0
      TML  = 0.0

C     Look up constituent parameters by matching the name
      DO 100 J = 1, NSPED
        LNAME = LABEL(J)
        SPEED = 0.0
        CALL NAME(SPEED, LNAME, ISUB, INUM, 2)
        SPD(J) = SPEED
        CALL VANDUF(SPEED,E,F,2)
        FFF(J) = F
        VOU(J) = E
100   CONTINUE

      YEARB  = IYR
      MONTHB = ICM
      DAYB   = ICD
      HOURB  = 0.0
      JDAY   = JULIAN(YEARB, MONTHB, DAYB, HOURB)
      JUDAY  = JDAY - JBASE_DATE + 1
      DO 200 J = 1, NSPED
        VUU(J) = VOU(J) + SPD(J)*(FLOAT(JUDAY-1)*24.0)
200   CONTINUE
      CALL TWOPI(VUU,NSPED)
      DAYBB = JUDAY
      DAYMM = JUDAY
      GRBSS = 0.0
      GRMSS = 0.0
      DAYMM = JUDAY + LENGTH/2
      CALL ASTRO(XYER, DAYBB, DAYMM, GRBSS, GRMSS)

      DO 277 IZ = 1, NSPED
        SPEED = SPD(IZ)
        CALL VANDUF(SPEED, E, F, 2)
        VAU(IZ) = E
        FFF(IZ) = F
277   CONTINUE

C     Round node to 3 decimals, vo+u to 1 place
      DO 222 MT = 1, NSPED
        FFF(MT) = ANINT(FFF(MT)*1000.0) * 0.001
        VAU(MT) = ANINT(VAU(MT)*10.0) * 0.1
222   CONTINUE

      RETURN
      END


      SUBROUTINE FITAN(AUS, AUC, RTA, SPDX, JMAP)
      IF(AUC .EQ. 0.0) THEN
        IF(AUS .LT. 0.0) THEN
          RTA = 270.0
        ELSEIF(AUS .EQ. 0.0) THEN
          RTA = 0.0
        ELSE
          RTA = 90.0
        END IF
      ELSE
        BXG = AUS / AUC
        RTA = ATAN(BXG) * 57.2957795
        IF(JMAP .GT. 0) THEN
          IF(AUS .LE. 0.0) THEN
            IF(AUC .LT. 0.0) THEN
              RTA = 180.0 +RTA
            ELSEIF(AUC .EQ. 0.0) THEN
              RTA = 0.0
            ELSE
              RTA = RTA + 360.0
            END IF
          ELSE
            IF(AUC .LT. 0.0) THEN
              RTA = 180.0 +RTA
            ELSEIF(AUC .EQ. 0.0) THEN
              RTA = 90.0
            END IF
          END IF
        END IF
      END IF
      SPDX = 0.0

      RETURN
      END


      SUBROUTINE TWOPI(AUG,IO)
      DIMENSION AUG(IO)

      DO 114 MO = 1, IO
        ZAT = AUG(MO)/360.0
        IF(ZAT.LT. 0.0) THEN
          AUG(MO) = ((ZAT - AINT(ZAT)) + 1.00)*360.0
        ELSE
          AUG(MO) = (ZAT - AINT(ZAT))*360.0
        END IF
114   CONTINUE

      RETURN
      END


      SUBROUTINE TWOPI0(AUG)

      ZAT = AUG/360.0
      IF(ZAT.LT. 0.0) THEN
        AUG = ((ZAT - AINT(ZAT)) + 1.00)*360.0
      ELSE
        AUG = (ZAT - AINT(ZAT))*360.0
      END IF

      RETURN
      END


      subroutine astro(xyer,dayb,daym,grbs,grms)
      implicit real*4(a-h,o-z)
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/vee/tml,con,u,q,ui
      common/boxs/aw,ai,ae,ae1,asp
      common/boxxs/vib,vb,xib,vpb,vppb,cxsb,cxpb,cxhb,cxp1b

      pinv = 57.29578
      nyear = ifix(xyer)
      call orbit(xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,5)
      xw = 23.4522944 - .0130125*t - .00000164*t**2 + .000000503*t**3
      xi = 5.14537628
      aw = xw*0.0174533
      ai = xi*0.0174533
      ae = 0.0548997
      ae1 = 0.01675104 - 0.0000418*t - .000000126*t**2
      asp = .46022931
      do 30 noe = 1,30
   30 cxx(noe) = 0.0
      if(dayb.gt.0.0) dayb = dayb - 1.0
      if(daym.gt.0.0) daym = daym - 1.0
      doby = 0.0
      amit = 0.0
      ami = xyer - xcen
      cplex = xcen/400.0 + 0.0001
      dicf = cplex - aint(cplex)
      if(ami.eq.0.0) go to 32
      xcet = xcen + 1.0
      cdif = xyer - xcet
      doby = cdif/4.0 + 0.0001
      amit = aint(doby)
      if(dicf.lt.0.001) amit = amit + 1.0
   32 farm = 0.25*ami
      farx = farm - aint(farm)
      cxx(1) = xsx + 129.384820*ami + 13.1763968*(dayb + amit) + 0.54901
     16532*grbs
      cxx(2) = xpx + 40.6624658*ami + 0.111404016*(dayb + amit) + 0.0046
     141834*grbs
      cxx(3) = xhx - 0.238724988*ami + 0.985647329*(dayb + amit) + 0.041
     1068639*grbs
      cxx(4) = xp1x + 0.01717836*ami + 0.000047064*(dayb + amit) + 0.000
     1001961*grbs
      cxx(5) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grms
      cxx(6) = xnx - 19.3281858*ami - 0.052953934*(daym + amit) - 0.0022
     106414*grms
   40 cxx(7) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grbs
      cxx(8) = xnx - 19.328185764*ami - 0.0529539336*(dayb + amit) - 0.0
     106414*grbs
      call twopi(cxx, 8)
   41 do 100 ii = 1,8
  100 cxx(ii) = float(ifix(cxx(ii)*100.0 + 0.5))*0.01
      ang = cxx(8)
      call table6(vib,vb,xib,vpb,vppb,xx,xx,xx,xx,xx,ang,anb,atb)
      cxx(26) = vib
      cxx(27) = vb
      cxx(28) = xib
      cxx(29) = vpb
      cxx(30) = vppb
      cxsb = cxx(1)
      cxpb = cxx(2)
      cxhb = cxx(3)
      cxp1b= cxx(4)
      ang = cxx(6)
      call table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
      cxx(9 ) = vi
      cxx(10) = v
      cxx(11) = xi
      cxx(12) = vp
      cxx(13) = vpp
  230 do 333 ii = 9,13
  333 cxx(ii) = float(ifix(cxx(ii)*100.0 + 0.5))*0.01
      pgx = cxx(5) - cxx(11)
      pgx = float(ifix(pgx*100.0 + 0.5))*0.01
      call twopi0(pgx)
      xpg = pgx*0.0174533
      cxx(14) = pgx
      raxe = sin(2.0*xpg)
      raxn = (cos(0.5*at)**2/(6.0*sin(0.5*at)**2)) - cos(2.0*xpg)
      rxx = 0.0
      if(raxe.eq.0.0.or.raxn.eq.0.0) go to 232
      rax = raxe/raxn
      if(rax.gt.3450.0) go to 232
        rxx   = atan(rax )*pinv
      cxx(22) = rxx
  232 cra = sqrt(1.0 - 12.0*(sin(0.5*at)**2/cos(0.5*at)**2)*cos(2.0*xpg)
     1 + 36.0*(sin(0.5*at)**4/cos(0.5*at)**4))
      um2 = 2.0*(cxx(11) - cxx(10))
      cxx(21) = um2
      cxx(24) = cra
      ul2 = um2 - rxx
  404 ul2 = ul2 + 180.0
  405 cxx(15) = ul2
      zes = (5.0*cos(aw) - 1.0)*sin(xpg)
      zec = (7.0*cos(aw) + 1.0)*cos(xpg)
      call fitan(zes,zec,qxx,spxx,2)
      cxx(23) = qxx
      crav = 0.5*um2 + qxx + 090.0
      cxx(16) = crav
      cqa = sqrt(0.25 + 1.5*((cos(aw)/cos(0.5*aw)**2)*cos(2.0*xpg)) + 2.
     125*(cos(aw)**2/cos(0.5*aw)**4))
      cxx(25) = cqa
      do 444 iii = 14,23
  444 cxx(iii) = float(ifix(cxx(iii)*100.0 + 0.5 ))*0.01
      pm   = cxx(1)
      pl   = cxx(2)
      sl   = cxx(3)
      ps   = cxx(4)
      plm  = cxx(5)
      skyn = cxx(6)
      vi   = cxx(9)
      v    = cxx(10)
      xi   = cxx(11)
      vp   = cxx(12)
      vpp  = cxx(13)
      p    = cxx(14)
      aul  = cxx(15)
      aum  = cxx(16)
      cra = cxx(24)
      cqa = cxx(25)
      u = v*0.0174533
      q = p*0.0174533
      ui = vi*0.0174533

      return
      end


C     Order of constituents is same as in NAMES common
      subroutine vanduf(speed,e,f,itype)
      real*8 spd
      double precision speed
      dimension spd(180),ms(180)

      common/locat/tm,gonl
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp
      common /speeds/spd
      common /mmss/ms

  500 format('0*** (v + u) not computed for constituent of speed',
     1 f12.7,'  ****',' Execution terminated')
      con = sl + tml
      do 600 j = 1,175
      ipick = j
      if(speed.eq.spd(j)) go to 611
  600 continue
      ipick = 176
  611 if(ipick.gt.164) go to 620
      if(ipick.gt.124) go to 618
      if(ipick.gt.88) go to 616
      if(ipick.gt.38) go to 614
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
     1,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38),ipick
  614 ipik = ipick - 38
      go to (39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58
     1,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80
     2,81,82,83,84,85,86,87,88),ipik
  616 ipikk = ipick - 88
      go to (89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,10
     16,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,
     2123,124),ipikk
  618 ipikkk = ipick - 124
      go to (125,126,127,128,129,130,131,132,133,134,135,136,137,138,139
     1,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,
     2156,157,158,159,160,161,162,163,164),ipikkk
  620 ipiikk = ipick - 164
      go to (165,166,167,168,169,170,171,172,173,174,175,488),ipiikk

    1 e = 2.0*(con - pm + xi - v)
  201 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))
      go to 888

    2 e =  2.0*(con + xi - v) - 3.0*pm + pl
      go to 201

    3 e = 2.0*tml
  203 f = 1.0

      go to 888
    4 e = con - v - 2.0*(pm - xi) - 90.0
      go to 218

    5 e = con - vp + 90.0
  205 f = 1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0
     1.1006)
      go to 888

    6 e = 2.0*con - vpp
  206 f = 1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) + 0
     1.0981)
      go to 888

    7 e = 2.0*con - pm - pl + aul
  207 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*(1.0/cra)
      go to 888

    8 e = 2.0*(con + xi - v + pl) - 4.0*pm
      go to 201

    9 e = sl - ps + 180.0 + 2.0*tml
      go to 203

   10 e = 2.0*tml - (sl - ps)
      go to 203

   11 e = 2.0*(con + xi - v - sl) - pm + pl + 180.0
      go to 201

   12 e = 2.0*(con + xi - v + sl) - 4.0*pm
      go to 201

   13 e = 2.0*(con + xi - v + sl) - 3.0*pm - pl
      go to 201

   14 e = con + pm - pl - v + 90.0
  214 f = (sin(2.0*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui)
      go to 888

   15 e = con - pm + aum
      f =  ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)
     1**2))*(1.0/cqa)
      go to 888
   16 e = con - v + 2.0*(pm - xi)+ 90.0
  216 f = (sin(aw)*sin(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*sin(.5*ui)**2)
      go to 888
   17 e = tml + 270.0 - sl
      go to 203
   18 e = con - v - 3.0*pm + 2.0*xi + pl - 90.0
  218 f = (sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
      go to 888
   19 e = con - v - 4.0*pm + 2.0*xi + 2.0*pl - 90.0
      go to 218
   20 e = con - v - 3.0*pm + 2.0*xi - pl + 2.0*sl - 90.0
      go to 218
   21 e = 4.0*(con - pm + xi - v)
  221 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2
      go to 888
   22 e = 6.0*(con - pm + xi - v)
  222 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3
      go to 888
   23 e = 8.0*(con - pm + xi - v)
  223 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4
      go to 888
   24 e = 4.0*tml
      go to 203
   25 e = 6.0*tml
      go to 203
   26 e = 3.0*(con - pm + xi - v) + 180.0
  226 f = (cos(0.5*aw )**6*cos(0.5*ai)**6)/cos(0.5*ui)**6
      go to 888
   27 e = tml + 180.0
      go to 203
   28 e = 2.0*(con - pm + xi - v) + ( con - vp + 90.0)
  228 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**1*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   29 e = 4.0*(con - pm + xi - v) - (con - vp + 90.0)
  229 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   30 e = 4.0*(con + xi - v) + pl - 5.0*pm
      go to 221
   31 e = 2.0*(con - pm + xi - v) + 2.0*tml
      go to 201
   32 e = 4.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   33 e = 2.0*(pm - xi)
  233 f = (sin(aw)**2*cos(0.5*ai)**4)/sin(ui)**2
      go to 888
   34 e = 2.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   35 e = pm - pl
  235 f = ((2./3.- sin(aw)**2)*(1.- 1.5*sin(ai)**2))/(2./3.- sin(ui)**2)
      go to 888
   36 e = sl
      go to 203
   37 e = 2.0*sl
      go to 203
   38 e = 2.0*(con - v - 2.0*(pm - xi) - 90.0) - (tml+270.0 - sl)
  238 f = ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)*
     1*2))**2
      go to 888
   39 e = 2.0*(con - pm + xi - v) - (tml + 270.0 - sl)
      go to 201
   40 e = con + 2.0*sl - v - pm - pl + 90.0
      go to 258
   41 e = 2.0*(tml + 270.0 - sl) - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   42 e = 2.0*tml - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   43 e = 2.0*tml + pm - pl
      go to 221
   44 e = 4.0*(con + xi - v) - 5.0*pm - 2.0*tml  + pl
      go to 221
   45 e = con - v - 2.0*(pm - xi) + tml - sl + 180.0
      go to 218
   46 e = 4.0*con - 2.0*(pm - xi + v + tml) - vpp
      go to 259
   47 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*(pl - tml)
      go to 221
   48 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + aul
      go to 261
   49 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) - pl + aul
      go to 261
   50 e = 2.0*(tml + pm - xi + v) - vpp
      go to 259
   51 e = 4.0*(xi - pm - v) + 2.0*(tml + vpp)
  251 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   52 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v) - pl - vpp + aul
  252 cra = sqrt(1.0 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 3
     16.0*(sin(.5*ui)**4/cos(.5*ui)**4))
      f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./cra)
      go to 888
   53 e = 6.0*con - 4.0*tml + 2.0*(xi - v - pm - vpp)
  253 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   54 e = 4.0*tml - 2.0*con + pl - pm + vpp
  254 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   55 e = 4.0*con - 2.0*(tml + vpp) + pm - pl
      go to 251
   56 e = 2.0*tml + con - v - 2.0*(pm - xi) - 90.0
      go to 218
   57 e = 2.0*tml + con - vp + 90.0
      go to 205
   58 e = 3.0*(con - v) + 4.0*xi - 5.0*pm + pl - 90.0
  258 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   59 e = 4.0*con - 2.0*(pm - xi + v) - vpp
  259 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   60 e = 2.0*(tml + con + xi - v) - 3.0*pm + pl
      go to 201
   61 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml - pl + aul
  261 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   62 e = 6.0*(con - pm + xi - v) - 2.0*tml
      go to 222
   63 e = 4.0*con - 3.0*pm + 2.0*(xi - v) - pl + aul
  263 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   64 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl
      go to 221
   65 e = 2.0*(con + tml) - pm - pl + aul
      go to 207
   66 e = 5.0*(con - v) - 7.0*pm + 6.0*xi + pl - 90.0
  266 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   67 e = 5.0*(con - v) + 6.0*(xi - pm) - 90.0
      go to 266
   68 e = 5.0*con + 4.0*(xi - pm - v) - vp + 90.0
      go to 229
   69 e = 3.0*con + 2.0*(xi - pm - v + tml) - vp + 90.0
      go to 228
   70 e = 5.0*con + 2.0*(xi - pm - v) - (vp + vpp) + 90.0
  270 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./sqrt(
     2.8965*sin(2.*ui)**2 + .6001*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
   71 e = 4.0*(con - pm + xi - v) + tml - sl + 270.0
      go to 221
   72 e = 6.0*(con - pm + xi - v) - tml + sl - 270.0
      go to 222
   73 e = 5.0*(con - pm) + 4.0*(xi - v) + pl - vp + 90.0
      go to 229
   74 e = 2.0*(con - pm + xi - v) + 4.0*tml
      go to 201
   75 e = 6.0*(con + xi - v) - 7.0*pm + pl
      go to 222
   76 e = 4.0*(con + xi - v) - 5.0*pm + 2.0*tml + pl
      go to 221
   77 e = 4.0*(con - pm + xi - v) + 2.0*tml
      go to 221
   78 e = 8.0*con - 9.0*pm + 6.0*(xi - v) - 2.0*tml + pl + aul
  278 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   79 e = 6.0*(con + xi - v) - 8.0*pm + 2.0*pl
      go to 222
   80 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) - pl + aul
      go to 263
   81 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - pl + aul
      go to 261
   82 e = 4.0*con + 2.0*(xi - pm - v + tml) - vpp
      go to 259
   83 e = 8.0*(con - pm) + 6.0*(xi - v) - 2.0*tml + aul
      go to 278
   84 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml - pl + aul
      go to 278
   85 e = 6.0*con + 4.0*(xi - pm - v) - vpp
      go to 254
   86 e = 7.0*(con - v) - 9.0*pm + 8.0*xi + pl - 90.0
  286 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   87 e = 7.0*con + 6.0*(xi - v) - 8.0*pm + 2.0*pl - vp + 90.0
  287 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   88 e = 5.0*(con - v) + 6.0*(xi - pm) + 2.0*tml - 90.0
      go to 266
   89 e = 5.0*con + 4.0*(xi - pm) - 3.0*v + 2.0*tml - vpp - 90.0
  289 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*( (sin(aw
     2)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   90 e = 6.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl
      go to 222
   91 e = 6.0*(con - pm + xi - v) + 2.0*tml
      go to 222
   92 e = 4.0*(con - pm + xi - v) + 4.0*tml
      go to 221
   93 e = 8.0*(con + xi - v) - 10.0*pm + 2.0*pl
      go to 223
   94 e = 8.0*(con + xi - v) - 9.0*pm + pl
      go to 223
   95 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + 2.0*tml - pl + aul
      go to 261
   96 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml - pl + aul
  296 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   97 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
   98 e = 8.0*con + 6.0*(xi - pm - v) - vpp
  298 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   99 e = 6.0*con + 4.0*(xi - pm - v) + 2.0*tml - vpp
      go to 254
  100 e = 9.0*con - 10.0*pm + 8.0*(xi - v) + 2.0*pl - vp + 90.0
  300 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  101 e = 9.0*(con - pm) + 8.0*(xi - v) + pl - vp + 90.0
      go to 300
  102 e = 9.0*con + 8.0*(xi - pm - v) - vp + 90.0
      go to 300
  103 e = 7.0*con + 6.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 287
  104 e = 10.0*(con + xi - v) - 11.0*pm + pl
  304 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5
      go to 888
  105 e = 10.0*(con - pm + xi - v)
      go to 304
  106 e = 8.0*(con + xi - v) - 9.0*pm + pl + 2.0*tml
      go to 223
  107 e = 8.0*(con - pm + xi - v) + 2.0*tml
      go to 223
  108 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
  109 e = 6.0*(con - pm + xi - v) + 4.0*tml
      go to 222
  110 e = 9.0*con + 8.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 300
  111 e = 10.0*(con + xi - v) - 11.0*pm + pl
      go to 304
  112 e = 10.0*(con - pm + xi - v) + 2.0*tml
      go to 304
  113 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml - pl + aul
      go to 296
  114 e = 8.0*(con - pm + xi - v) + 4.0*tml
      go to 223
  115 e = tml - 2.0*sl + ps - 90.0
      go to 203
  116 e = con + sl - ps + 90.0
      go to 203
  117 e = con + 2.0*sl + 90.0
      go to 203
  118 e = tml + pm - sl + pl - v + 90.0
      go to 258
  119 e = 2.0*con + pm - v - vp - pl + 180.0
  319 f =  ((sin(2.*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui))*(1./sqrt(.8
     1965*sin(2.*ui)**2 + .6011*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
  120 e = 2.0*(con - v) - 5.0*pm + 4.0*xi + pl + 180.0
      go to 238
  121 e = 3.0*(con - v) - 4.0*(pm - xi) - 90.0
      go to 258
  122 e = 2.0*(con + tml) - vpp
      go to 206
  123 e = con + 2.0*(pm - xi - vp) + v + 270.0
  323 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**1*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))**2
      go to 888
  124 e = con - 2.0*v - 4.0*(pm - xi) + vp - 270.0
  324 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**2*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))
      go to 888
  125 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + 2.0*pl - vpp
  325 go to 254
  126 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) + pl - vpp
  326 go to 254
  127 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v - vpp) + pl
  327 go to 253
  128 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  328 go to 254
  129 e = 4.0*con - 3.0*pm + 2.0*(xi - v) + pl - vpp
  329 go to 259
  130 e = 8.0*con - 9.0*pm + 6.0*(xi - v) + 3.0*pl - 2.0*tml - vpp
  330 go to 298
  131 e = 8.0*(con - pm) + 6.0*(xi - v) + 2.0*(pl - tml) - vpp
  331 go to 298
  132 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml + pl - vpp
  332 go to 298
  133 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + pl - vpp
  333 go to 254
  134 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) + pl - vpp
  334 go to 259
  135 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml + pl - vpp
  335 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
  136 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + pl - vpp
  336 go to 298
  137 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  337 go to 254
  138 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl - vpp
  338 go to 298
  139 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml + pl - vpp
  339 go to 335
  140 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl + sl - tml - 270.0
  340 go to 221
  141 e = 3.0*sl - 2.0*vp - tml - 90.0
  341 f = 1.0/(.8965*sin(2.0*ui)**2 + .6001*sin(2.0*ui)*cos(u) + .1006)
      go to 888
  142 e = tml + vp - 3.0*sl + 90.0
  342 go to 205
  143 e = 2.0*tml - vp
  343 go to 205
  144 e = 2.0*(tml - sl) + vpp
  344 go to 206
  145 e = 2.0*(tml - vpp) + 4.0*sl
  345 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**2
      go to 888
  146 e = 2.0*tml - 2.0*(sl - ps)
  346 go to 203
  147 e = 4.0*tml + ps - sl
  347 go to 203
  148 e = 4.0*tml - 2.0*sl + vpp
  348 go to 206
  149 e = 4.0*tml + 6.0*sl - 3.0*vpp
  349 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**3
      go to 888
  150 e = 4.0*tml - 3.0*(sl - ps)
  350 go to 203
  151 e = con - v - 2.0*(pm - xi) + tml + 90.0
  351 go to 218
  152 e = 2.0*(pm - xi) + v + tml - con - 90.0
  352 go to 218
  153 e = pm - xi - 90.0
  353 f = 0.31920/(sin(ui) - sin(ui)**3)
      go to 888
  154 e = pm - sl
  354 go to 235
  155 e = sl - ps
  355 go to 203
  156 e = 3.0*sl
  356 f = 1.0
      go to 888
  157 e = 4.0*sl
  357 go to 356
  158 e = 6.0*sl
  358 go to 356
  159 e = 8.0*sl
  359 go to 356
  160 e = 10.0*sl
  360 go to 356
  161 e = 12.0*sl
  361 go to 356
  162 e = 24.0*sl
  362 go to 356
  163 e = 3.0*tml + 180.0
  363 go to 356
  164 e = 5.0*tml + 180.0
  364 go to 356
  165 e = 2.0*(con - v) - 4.0*(pm - xi) - 180.0
  365 f =((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
     1**2
      go to 888
  166 e = tml + con - vp + 270.0
  366 go to 205
  167 e = 3.0*(con - pm) + 2.0*(xi - v) + pl - vp + 90.0
  367 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/cos(0.5*ui)**4 )*(1.0/sqrt(0.
     18965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  168 e = 3.0*tml + 270.0 - sl
  368 go to 203
  169 e = 3.0*con - vpp - vp + 90.0
  369 f = 1.0/(sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981)*sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     20.1006))
      go to 888
  170 e = 4.0*(con - v) + 6.0*xi - 7.0*pm + pl - 180.0
  370 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*((sin(aw)*c
     1os(0.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))**2
      go to 888
  171 e = 4.0*(con - v) - 6.0*(pm - xi) - 180.0
  371 go to 370
  172 e = 2.0*(tml + con - v) - 4.0*(pm - xi) - 180.0
  372 go to 365
  173 e = 7.0*tml + 180.0
  373 go to 356
  174 e = 8.0*tml
  374 go to 356
  175 e = 10.0*tml
  375 go to 356
  488 print 500, speed
      stop
  888 go to (400,390),itype
  390 e = e + float(ms(ipick))*gonl - spd(ipick)*(tm/15.0)
  400 call twopi0(e)
      f = 1.0/f

      return
      end


      subroutine name(spdd,itag,isub,inum,icode)
c     this subroutine identifies the constituent by its speed,
c     name label or constituent number,
c     and makes it available for labeling.
c     ICODE = 1, by speed
c     ICODE = 2, by label
c     ICODE = 3, by number
c     it also determine the subscript of the constituent
c        order of constituent speeds***  m(2),n(2),s(2),o(1),k(1),k(2)
c      l(2),2n(2)r(2),t(2),lambda(2),mu(2),nu(2),j(1),m(1),oo(1),p(1)
c      q(1),2q(1),rho(1),m(4),m(6),m(8),s(4),s(6),m(3),s(1),mk(3),2mk(3)
c      mn(4),ms(4),2sm(2),mf,msf,mm,sa,ssa
      real*8 spd,spdd
      character*10 lable(180)*10,itag
      dimension spd(180),ip(180)

      common /mmss/ip
      common /speeds/spd
      common /names/lable

    1 format(10x,'Constituent of speed ',f12.7,' not in list.')
    2 format(10x,'Constituent ', a10,' not in the list.')
    3 format(10x,'Constituent no.',i4,' not in list.')

      if (icode.eq.1) then
*       search by speed
        do 100 j = 1,175
        if(spdd.ne.spd(j)) go to 100
        itag = lable(j)
        isub = ip(j)
        inum = j
        return
100     continue
        print 1,spdd
      else if (icode.eq.2) then
*       search by name
        do 200 i = 1,175
        if(itag.ne.lable(i)) go to 200
        spdd = spd(i)
        isub = ip(i)
        inum = i
        return
200     continue
        print 2, itag
      else if (icode.eq.3) then
        write(*,*)'get wrong icode, icode must be 1 or 2'
        stop
        print 3, inum
      end if

      stop '**** Execution terminated in NAME (illegal icode) ****'

      end


      subroutine orbit(xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,nnn)
      implicit real*4(a-h,o-z)
      dimension oex(nnn)
      s = 13.1763968
      p = 0.1114040
      xh = 0.9856473
      p1 = 0.0000471
      xn = -.0529539
      xcan = xyer*0.01 + 0.001
      xcen = aint(xcan)*100.0
      t = -3.0
      yr = 2.5
      gat = 1600.0
      do 10 jk = 1,30
      gp = gat/400.0 + 0.00001
      col = gp - aint(gp)
      if(col.lt.0.010) go to 11
      if(gat.eq.xcen) go to 12
      yr = yr - 1.0
      go to 9
   11 if( gat.eq.xcen) go to 12
    9 gat = gat + 100.0
   10 continue
   12 t = (gat - 1900.0)*0.01
      oex(1) = 270.437422 + 307.892*t + 0.002525*t**2 + .00000189*t**3 +
     1 yr*s
      oex(2) = 334.328019 + 109.032206*t - 0.01034444*t**2 - .0000125*t*
     1*3 + yr*p
      oex(3) = 279.696678 + 0.768925*t + .0003205*t**2 + yr*xh
      oex(4) = 281.220833 + 1.719175*t + 0.0004528*t**2 + .00000333*t**3
     1 + yr*p1
      oex(5) = 259.182533 - 134.142397*t + .00210556*t**2 + .00000222*t*
     1*3 + yr*xn
      call twopi(oex,5)
      do 100 i = 1,5
  100 oex(i) = float(ifix(oex(i)*100.0 + 0.5))*0.01
      xsx = oex(1)
      xpx = oex(2)
      xhx = oex(3)
      xp1x = oex(4)
      xnx = oex(5)

      return
      end


      subroutine table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
      common/boxs/aw,ai,ae,ae1,asp
      v = 0.0
      xi = 0.0
      vp = 0.0
      vpp = 0.0
      an = ang*0.0174533
      ax = ang
      eye = cos(ai)*cos(aw) - sin(ai)*sin(aw)*cos(an)
      c9 = acos(eye)*57.2957795
      vi = float(ifix(c9*100.0 + 0.5))*0.01
      cig = vi*0.0174533
      at = cig
      if(cig.eq.0.0) go to 230
      if(ax.eq.0.0.or.ax.eq.180.0) go to 230
      vxxe = sin(ai)*sin(an)
      vxxn = cos(ai)*sin(aw) + sin(ai)*cos(aw)*cos(an)
      if(vxxe.eq.0.0.or.vxxn.eq.0.0) go to 201
      vxx = vxxe/vxxn
      c10 = atan(vxx)*57.2957795
      v = float(ifix(c10*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.v.gt.0.0) v = -1.0*v
  201 cvx = v*0.0174533
      term = sin(ai)*(cos(aw)/sin(aw))
      exx = term*(sin(an)/cos(an)) + (cos(ai) - 1.0)*sin(an)
      if(exx.eq.0.0) go to 202
      ezz = term + cos(ai)*cos(an) + (sin(an)**2/cos(an))
      if(ezz.eq.0.0) go to 202
      exez = exx/ezz
      if(exez.gt.3450.0) go to 202
      c11 = atan(exez)*57.2957795
      xi = float(ifix(c11*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.xi.gt.0.0) xi = -1.0*xi
  202 cex = xi*0.0174533
      a22 = (0.5 + 0.75*ae**2)*sin(2.0*cig)
      b22 = (0.5 + 0.75*ae1**2)*sin(2.0*aw)*asp
      vpxe = a22*sin(cvx)
      vpxn = a22*cos(cvx) + b22
      if(vpxe.eq.0.0.or.vpxn.eq.0.0) go to 203
      vpx = vpxe/vpxn
      if(vpx.gt.3450.0) go to 203
      vp = atan(vpx)*57.2957795
      if(ax.gt.180.0.and.vp.gt.0.0) vp = -1.0*vp
  203 pvc = vp*0.0174533
      a47 = (0.5 + 0.75*ae**2)*sin(cig)**2
      b47 = (0.5 + 0.75*ae1**2)*asp*sin(aw)**2
      vpye = a47*sin(2.0*cvx)
      vpyn = a47*cos(2.0*cvx) + b47
      if(vpye.eq.0.0.or.vpyn.eq.0.0) go to 204
      vpy = vpye/vpyn
      if(vpy.gt.3450.0) go to 204
      vpp = atan(vpy)*57.2957795
      if(ax.gt.180.0.and.vpp.gt.0.0) vpp = -1.0*vpp
  204 pvcp = vpp*0.0174533

  230 return
      end

