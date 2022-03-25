CC  Write out a table
      SUBROUTINE TABLE_TIDE(FILESHORT,STATIONNAME,KINDAT)
      include 'skills.inc'

      CHARACTER*120 STATIONNAME,FILESHORT
      CHARACTER*20 CMETHOD(2),CGAP(2),CBDATE,CEDATE
      CHARACTER TITLE(6)*36,C1*1,C2*1,TTL*9,TTL2*15
      DIMENSION B(IMX),E(IMX),R(IMX),P(IMX),TZ(IMX),
     1  IFRST(IMX),ILAST(IMX),
     2  HHWP0(NMX),THWP0(NMX),HLWP0(NMX),TLWP0(NMX),
     3  THIGHSR(NMX),HHIGHSR(NMX),IDXR(NMX),
     4  THIGHSP(NMX),HHIGHSP(NMX),IDXP(NMX),
     5  THIGHS0(NMX),HHIGHS0(NMX),IDX0(NMX)
      DIMENSION DIRE(IMX),DIRR(IMX),DIRP(IMX),DIROP(IMX),DIRAP(IMX),
     1  AFCR(NMX),AECR(NMX),TFCR(NMX),TECR(NMX),DFCR(NMX),DECR(NMX),
     2  TSFR(NMX),TEFR(NMX),TSER(NMX),TEER(NMX),   
     3  AFCP(NMX),AECP(NMX),TFCP(NMX),TECP(NMX),DFCP(NMX),DECP(NMX),
     4  TSFP(NMX),TEFP(NMX),TSEP(NMX),TEEP(NMX)   
      DIMENSION OP(IMX),AP(IMX),TP(IMX),P2(IMX),XTMP(IMX),YTMP(IMX)
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE,JULIAN
      REAL*8 YEARB,MONTHB,DAYB,HOURB

      DATA TITLE/'SCENARIO: TIDAL SIMULATION ONLY     ',
     1           'SCENARIO: HINDCAST                  ',
     2           'SCENARIO: SEMI-OPERATIONAL NOWCAST  ',
     3           'SCENARIO: SEMI-OPERATIONAL FORECAST ',
     4           'COMPARISON: PERSISTENCE FORECAST    ',
     5           'COMPARISON: ASTRONOMICAL TIDE ONLY  '/
      DATA CMETHOD/'Cubic Spline','SVD         '/
      DATA CGAP/'Gap is not filled','Gap is filled    '/

      WRITE(*,*) 'PROGRAM TABLE_TIDE'
      II = 0
      DO I = 1, IMAXB
        IF(OF(I) .GT. -900.0) THEN
          II = II+1
          TP(II) = T(I)
        ENDIF
      ENDDO
      NMAX1 = II

      GAP = 1.5*DELT
      CALL CONTINUOUS(TP,NMAX1,GAP,NSEGMENTS,IFRST,ILAST)
      DIFM = -999.0
      WRITE(*,*) 'Nsegments = ', NSEGMENTS
      DO NG = 1, NSEGMENTS
        ISTART = IFRST(NG)
        IEND = ILAST(NG)
        DIF = TP(IEND)-TP(ISTART)
        IF(DIF .GT. DIFM) THEN
          DIFM = DIF
          II = NG
        ENDIF  
      ENDDO
      WRITE(*,*) 'II=',II,' The longest segment=',DIFM,' hours'

      AGAP = 0.0
      IF(NSEGMENTS .GT. 1) THEN
        DO NG = 1, NSEGMENTS-1
          DIF = TP(IFRST(NG+1))-TP(ILAST(NG))
          AGAP = AGAP+DIF
        ENDDO
        AGAP = AGAP/24.0
      ENDIF
      WRITE(*,*) 'Total gap in observations=',AGAP,' days'

      IYR = IYEAR
      ISTART = IFRST(II)
      IEND = ILAST(II)
      YEARB = IYEAR
      MONTHB = 1.0
      DAYB = 1.0
      HOURB = 0.0
      JDAY0 = JULIAN(YEARB,MONTHB,DAYB,HOURB)

      JDAY = TP(1)/24.0+JDAY0-1
      CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
      IYRS = INT(YEARB)
      ICMS = INT(MONTHB+0.001)
      ICDS = INT(DAYB+0.001)
      IHRS = INT(HOURB+0.001)
      IMNS = INT((HOURB-IHRS)*60+0.1)

      JDAY = TP(NMAX1)/24.0+JDAY0-1
      CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
      IYRE = INT(YEARB)
      ICME = INT(MONTHB+0.001)
      ICDE = INT(DAYB+0.001)
      IHRE = INT(HOURB+0.001)
      IMNE = INT((HOURB-IHRE)*60+0.1)
      LU = 10
      LU1 = 12

      OPEN(LU,FILE=TRIM(FILESHORT)//'_table.out',FORM='FORMATTED')
CC  Write column header
      WRITE(*,"(/,' Station: ',a40)") TRIM(STATIONNAME)
      WRITE(*,"(' Variable      X       N     Imax     SM    RMSE    SD',
     1  '     NOF   CF    POF   MDNO  MDPO  WOF  SKILL',
     2   /' Criterion    -     -      -      -  ',
     3  '    -       -     <1%  >90%   <1%    <N    <N  <.5%')")

      WRITE(LU,"(/,' Station: ',a40)") TRIM(STATIONNAME)
      WRITE(LU,"(' Observed data time period from: /',
     1   I2.2,'/',I2.2,'/',I4.4,'  to /',
     2   I2.2,'/',I2.2,'/',I4.4, ' with gaps of ',F7.2, ' days')") ICMS,
     3   ICDS,IYRS,ICME,ICDE,IYRE,AGAP
      IF(KINDAT .EQ. 1) THEN
        OPEN(LU1,FILE=TRIM(FILESHORT)//'phase_table.out')
        WRITE(LU1,"('Station: ',a40)") TRIM(STATIONNAME)
        WRITE(LU1,"(' Observed data time period from: /',
     1     I2.2,'/',I2.2,'/',I4.4,'  to /',
     2     I2.2,'/',I2.2,'/',I4.4, 'with gaps of ',F7.2, ' days')") ICMS,
     3     ICDS,IYRS,ICME,ICDE,IYRE,AGAP
      ENDIF

      IF(IGAPFILL .EQ. 0) THEN
        WRITE(LU,"(' Data gap is not filled')")
        IF(KINDAT .EQ. 1) THEN
          WRITE(LU1,"(' Data gap is not filled')")
        ENDIF
      ELSE
        IF(METHOD .EQ. 0) THEN
          WRITE(LU,"(' Data gap is filled by cubic spline method')")
          IF(KINDAT .EQ. 1) THEN
            WRITE(LU1,"(' Data gap is filled by cubic spline method')")
          ENDIF
        ELSE IF(METHOD .EQ. 1) THEN
          WRITE(LU,"(' Data gap is filled by SVD method')")
          IF(KINDAT .EQ. 1) THEN
            WRITE(LU1,"(' Data gap is filled by SVD method')")
          ENDIF
        ENDIF
      ENDIF  

      IF(TCUT .GT. 0.0) THEN
        WRITE(LU,"(' Data are filtered using ',F5.1,
     1     ' Hour Fourier Filter')") TCUT
        IF(KINDAT .EQ. 1) THEN
          WRITE(LU1,"(' Data are filtered using ',F5.1,
     1       ' Hour Fourier Filter')") TCUT
        ENDIF
      ELSE
        WRITE(LU,"(' Data are not filtered')")
        IF(KINDAT .EQ. 1) THEN
          WRITE(LU1,"(' Data are not filtered')")
        ENDIF
      ENDIF

      WRITE(LU,"('------------------------------------------',
     1   '-------------------------------------------------')")
      WRITE(LU,"(' VARIABLE     X       N    IMAX     SM    RMSE    SD',
     1   '    NOF   CF    POF   MDNO  MDPO  WOF   CORR  SKILL',
     2   /' CRITERION    -       -      -       -  ',
     3   '    -      -    <1%  >90%   <1%    <N    <N  <.5%')")
      WRITE(LU,"('------------------------------------------',
     1   '-------------------------------------------------',/)")
      IF(KINDAT .EQ. 1) THEN
        WRITE(LU1,"('------------------------------------------',
     1     '-------------------------------------------------')")
        WRITE(LU1,"(' VARIABLE     X      N    IMAX    SM    RMSE    SD',
     1     '    NOF   CF    POF   MDNO  MDPO  WOF  CORR  SKILL',
     2     /' CRITERION    -      -      -      -  ',
     3     '    -      -    <1%  >90%   <1%    <N    <N  <.5%')")
        WRITE(LU1,"('------------------------------------------',
     1     '-------------------------------------------------',/)")
      ENDIF

CC  Loop thru scenarios
      DO 220 IS = 1, 6
        IF(ISWITCH(IS) .EQ. 0) GOTO 220
CC  Write header
        WRITE(*,"(/,5x,a36)") TITLE(IS)
        WRITE(lu,"(5x,a36)") TITLE(IS)

        IF(NMAX1 .LT. 2*IPDAY) THEN
          WRITE(LU,"(5x,a36)") 'No observational data'
          GOTO 220
        ENDIF

        IF(KINDAT .EQ. 1) THEN
          WRITE(LU1,"(5x,a36)") TITLE(IS)
        ENDIF

CC  Compute error/difference array/series
        IF(IS .EQ. 4 .OR. IS .EQ. 5) GOTO 130
        IF(IS .EQ. 1) THEN
          NMAX = IMAXB
          II = 0
          DO I = 1, NMAX
            IF((A(I) .GT. -900.0) .AND. (S(I) .GT. -900.0)) THEN
              II = II+1
              R(II) = A(I)
              P(II) = S(I)
              YTMP(II) = OF(I)
              E(II) = P(II)-R(II)
              IF(KINDAT .EQ. 1) THEN
                IF((A(I) .GT. 0.26) .AND. (S(I) .GT. 0.26)) THEN
                  DIRR(II) = DIRA(I)
                  DIRP(II) = DIRS(I)
                  DIRE(II) = DIRP(II)-DIRR(II)
                  IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                  IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
                ELSE 
                  DIRR(II) = DIRA(I)
                  DIRP(II) = DIRS(I)
                  DIRE(II) = 0.0
                ENDIF
              ENDIF                
              TP(II) = T(I)
            ENDIF
            IF(I .LE. 10) THEN
              WRITE(*,"(' is=1 i=',i5,' t=',f11.4,' p=',f9.4,' r=',
     1          f9.4,' e=',f9.4)") ii,tp(ii),p(ii),r(ii),e(ii)
            ENDIF
          ENDDO
          NMAX = II
          WCOF = WOF(P,R,R,IMX,NMAX,2.0*X1)

CC  Test: sum of cf, pof, nof = 1
          WRITE(*,"('Test: nof(x),cf(x),pof(x)=',3f8.3,' sum=',f8.3)")
     1       VOF(E,IMX,NMAX,X1),CF(E,IMX,NMAX,X1),POF(E,IMX,NMAX,X1),
     2       VOF(E,IMX,NMAX,X1)+CF(E,IMX,NMAX,X1)+POF(E,IMX,NMAX,X1)

        ELSE IF(IS .EQ. 2) THEN   !!! model hindcast
          NMAX = IMAXB
          II = 0
          DO I = 1, NMAX
            IF((OF(I) .GT. -900.) .AND. (HI(I) .GT. -900.)) THEN
              II = II+1
              R(II) = OF(I)
              P(II) = HI(I)
              E(II) = P(II)-R(II)
              TP(II) = T(I)
              YTMP(II) = A(I)
              IF(KINDAT .EQ. 1) THEN
                IF((OF(I) .GT. 0.26) .AND. (HI(I) .GT. 0.26)) THEN
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRHI(I)
                  DIRE(II) = DIRP(II)-DIRR(II)
                  IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                  IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
                ELSE 
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRHI(I)
                  DIRE(II) = 0.0
                ENDIF
              ENDIF                
            ENDIF
            IF(I .LE. 10) THEN
              WRITE(*,"(' is=2 i=',i5,' t=',f11.4,' p=',f9.4,' r=',
     1          f9.4,' e=',f9.4)") I,T(I),P(I),R(I),E(I)
            ENDIF
          ENDDO
          NMAX = II
          WCOF = WOF(P,YTMP,R,IMX,NMAX,2.0*X1)

CC  Test: sum of cf, pof, nof = 1
          WRITE(*,"('Test: nof(x),cf(x),pof(x)=',3f8.3,' sum=',f8.3)")
     1       VOF(E,IMX,NMAX,X1),CF(E,IMX,NMAX,X1),POF(E,IMX,NMAX,X1),
     2       VOF(E,IMX,NMAX,X1)+CF(E,IMX,NMAX,X1)+POF(E,IMX,NMAX,X1)

        ELSE IF(IS .EQ. 3) THEN   !!! nowcast
          NMAX = IMAXB
          II = 0
          DO I = 1, nmax
            IF((OF(I) .GT. -900.0) .AND. (C(I) .GT. -900.0)) THEN
              II = II+1
              R(II) = OF(I)
              P(II) = C(I)
              E(II) = P(II)-R(II)
              TP(II) = T(I)
              YTMP(II) = A(I)
              IF(KINDAT .EQ. 1) THEN
                IF((OF(I) .GT. 0.26) .AND. (C(I) .GT. 0.26)) THEN
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRC(I)
                  DIRE(II) = DIRP(II)-DIRR(II)
                  IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                  IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
                ELSE 
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRC(I)
                  DIRE(II) = 0.0
                ENDIF
              ENDIF                
            ENDIF
C            IF(I .LE. 10) THEN
C              WRITE(*,"('is=3 i=',i5,' tp=',f11.4,' p=',f9.4,' r=',
C     1          f9.4,' e=',f9.4)") I,TP(I),P(I),R(I),E(I)
C            ENDIF
          ENDDO
          NMAX = II
          WCOF = WOF(P,YTMP,R,IMX,NMAX,2.0*X1)

        ELSE IF(IS .eq. 6) THEN
          NMAX = IMAXB
          II = 0
          DO I = 1, NMAX
            IF(OF(I) .GT. -900.0) THEN
              II = II+1
              R(II) = OF(I)
              P(II) = A(I)
              E(II) = P(II)-R(II)
              TP(II) = T(I)
              IF(KINDAT .EQ. 1) THEN
                IF((OF(I) .GT. 0.26) .AND. (A(I) .GT. 0.26)) THEN
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRA(I)
                  DIRE(II) = DIRP(II)-DIRR(II)
                  IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                  IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
                ELSE 
                  DIRR(II) = DIRO(I)
                  DIRP(II) = DIRA(I)
                  DIRE(II) = 0.0
                ENDIF
              ENDIF                
            ENDIF
          ENDDO
          NMAX = II
          WCOF = WOF(P,P,R,IMX,NMAX,2.0*X1)
        ENDIF

CC  Create series
        IF(IS .LT. 5) THEN   !!! for is=1:tide only; =2:hindcast; =3:nowcast
          TTL2 = '         '
          GAP = 1.5*DELT

          IF(KINDAT .EQ. 1) THEN
            call prtline1('U        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('u        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
            call prtline1('D        ',ttl2,lu1,dirp,tp,imx,nmax,
     1           X11,1,gap)
            call prtline1('d        ',ttl2,lu1,dirr,tp,imx,nmax,
     1           X11,1,gap)
          ELSE IF(KINDAT .eq. 2) THEN
            call prtline1('H        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('h        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          ELSE IF(KINDAT .eq. 3) THEN
            call prtline1('T        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('t        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          ELSE IF(KINDAT .eq. 4) THEN
            call prtline1('S        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('s        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          ENDIF
        ENDIF

        TTL = 'H-htest'
        IMAX = NMAX
        GAP = 1.5*DELT
        IF(KINDAT .EQ. 1) THEN
          SKILLV = SKILLA(P,R,IMX,NMAX)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          IX1 = INT(X1*100+0.1)
          WRITE(TTL2,'(a3,I3,a9)') '   ',IX1,' cm/s 24h'
          call prtline1('U-u      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)
          WRITE(TTL2,'(a2,F4.1,a9)') '  ',X11,'   dg 24h'
          SKILLV = SKILLA(DIRP,DIRR,IMX,NMAX)
          CORR_C=CORRELATION(DIRP,DIRR,IMX,NMAX)
          call prtline1('D-d      ',ttl2,lu1,dire,tp,imx,nmax,
     1         X11,6,gap)

        ELSE IF(KINDAT .EQ. 2) THEN
          IX1 = INT(X1*100+0.1)
          WRITE(TTL2,'(a3,I3,a9)')'   ',IX1,'   cm 24h'
          SKILLV = SKILLA(P,R,IMX,NMAX)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          call prtline1('H-h      ',ttl2,lu,e,tp,imx,nmax,X1,5,gap)

        ELSE IF(KINDAT .EQ. 3) THEN
          WRITE(TTL2,'(a3,F3.1,a9)')'   ',X1,'    c 24h'
          SKILLV = SKILLA(P,R,IMX,NMAX)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          call prtline1('T-t      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)

        ELSE IF(KINDAT .EQ. 4) THEN
          WRITE(TTL2,'(a3,F3.1,a9)')'   ',X1,'  psu 24h'
          SKILLV = SKILLA(P,R,IMX,NMAX)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          call prtline1('S-s      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)
        ENDIF 
        GOTO 170
 130    CONTINUE
        IF(IS .EQ. 4) INDX1 = 2  ! fcst
        IF(IS .EQ. 5) INDX1 = 3  ! persistence
        NLOOP = INT(NFDURATION/6+0.001) 

CC Use 4 replace jmax to do 00z,06z, 12z, 18z and 24z no matter what jmax is.
CC Extended forecast hours into 30z, 36z,etc. for the skill assessment 
        DO 160 JJ = 0, NLOOP  
          CALL ADDPRED(INDX1,2,YTMP,DIRP,XTMP,IPMAX,JJ)
          WRITE(*,*) 'IPMAX IN ADDPRED=', IPMAX
          CALL COLLECT_OBS(OP,AP,DIROP,DIRAP,XTMP,TZ,IPMAX,NMAX)
          IF(NMAX .NE. IPMAX) THEN
            WRITE(*,*) ' nmax=',nmax,' unequal to ipmax=',ipmax
            WRITE(*,*) 'program stop here'
!          stop
            DO I = 1, NMAX
              DO I1 = 1, IPMAX
                IF(ABS(TZ(I)-XTMP(I1)) .LE. 0.001) THEN
                  YTMP(I) = YTMP(I1)
                  XTMP(I) = XTMP(I1)
                  GOTO 135
                ENDIF
              ENDDO
135           CONTINUE
            ENDDO  
          ENDIF
          II = 0
          DO I = 1, NMAX
            IF(KINDAT .GT. 1) THEN
              IF((OP(I) .GT. -900.0) .AND. (YTMP(I) .GT. -900.0)) THEN
                II = II+1
                R(II) = OP(I)
                P(II) = YTMP(I)
                TP(II) = XTMP(I)
                E(II) = P(II)-R(II)
                AP(II) = AP(I)
C	        WRITE(*,'(4F10.4,2x,I4)') tP(II)/24.0,P(II),R(II),
C     1           E(II),II
	      ENDIF
            ELSEIF(KINDAT .EQ. 1) THEN
              IF((OP(I) .GT. 0.26) .AND. (YTMP(I) .GT. 0.26)) THEN
                II = II+1
                R(II) = OP(I)
                P(II) = YTMP(I)
                TP(II) = XTMP(I)
                AP(II) = AP(I)
                E(II) = P(II)-R(II)
                DIRP(II) = DIRP(I)
                DIROP(II) = DIROP(I)
                DIRAP(II) = DIRAP(I)
                DIRE(II) = DIRP(II)-DIROP(II)
                IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
C               WRITE(*,*) 'II=',II,TP(II)/24.0,P(II),R(II),E(II),
C     1	           DIRP(II),DIROP(II),DIRE(II)
              ENDIF
            ENDIF                
          ENDDO
          NMAX = II

          C1 = CHAR(ICHAR('0')+(JJ*6)/10)
          C2 = CHAR(ICHAR('0')+(JJ*6)-10*((JJ*6)/10))
          IF(KINDAT .EQ. 1) THEN
            IX1 = INT(X1*100+0.1)
            WRITE(TTL2,'(a3,I3,a9)') ' ',IX1,' cm/s 24h'
            TTL = 'U'//c1//c2//'-'//'u'//c1//c2
            WRITE(TTL,3234) 'U',JJ*6,'-u',JJ*6
            WCOF = WOF(P,AP,R,IMX,NMAX,2.0*X1)
            SKILLV = SKILLA(P,R,IMX,NMAX)
            CORR_C = CORRELATION(P,R,IMX,NMAX)           
            GAP = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)
            WRITE(TTL2,'(a2,F4.1,a9)') '  ',X11,'   dg 24h'
            TTL = 'D'//c1//c2//'-'//'d'//c1//c2
            WRITE(TTL,3234) 'D',JJ*6,'-d',JJ*6
            SKILLV = SKILLA(DIRP,DIROP,IMX,NMAX)
            CORR_C = CORRELATION(DIRP,DIROP,IMX,NMAX)            

            call prtline1(ttl,ttl2,lu1,dire,tp,imx,nmax,X11,3,gap)

          ELSE IF(KINDAT .EQ. 2) THEN
            IX1 = INT(X1*100+0.1)
            WRITE(TTL2,'(a3,I3,a9)') '   ',IX1,'   cm 24h'
            TTL = 'H'//c1//c2//'-'//'h'//c1//c2
            WRITE(TTL,3234) 'H',JJ*6,'-h',JJ*6
3234	    FORMAT(A1,I3.3,A2,I3.3)
            WCOF = WOF(P,AP,R,IMX,NMAX,2.0*X1)
            SKILLV = SKILLA(P,R,IMX,NMAX)
            CORR_C = CORRELATION(P,R,IMX,NMAX)  
            GAP = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,2,gap)

          ELSE IF(KINDAT .EQ. 3) THEN
            WRITE(TTL2,'(a3,F3.1,a9)') '   ',X1,'    c 24h'
            TTL = 'T'//c1//c2//'-'//'t'//c1//c2
            WRITE(TTL,3234) 'T',JJ*6,'-t',JJ*6
            SKILLV = SKILLA(P,R,IMX,NMAX)
            CORR_C = CORRELATION(P,R,IMX,NMAX)  

            GAP = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)

          ELSE IF(KINDAT .EQ. 4) THEN
          WRITE(TTL2,'(a3,F3.1,a9)') '   ',X1,'  psu 24h'
            TTL = 'S'//c1//c2//'-'//'s'//c1//c2
            WRITE(TTL,3234) 'S',JJ*6,'-s',JJ*6
            GAP = 25.0
            SKILLV = SKILLA(P,R,IMX,NMAX)
            CORR_C = CORRELATION(P,R,IMX,NMAX)  
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)
          ENDIF
 160    CONTINUE

        IF((KINDAT .EQ. 3) .OR. (KINDAT .EQ. 4)) GOTO 220
        CALL ADDPRED(INDX1,1,YTMP,DIRP,XTMP,IPMAX,0)
        WRITE(*,*) 'ipmax in addpred=',IPMAX
        CALL COLLECT_OBS(OP,AP,DIROP,DIRAP,XTMP,TZ,IPMAX,NMAX)
        IF(NMAX .NE. IPMAX) THEN
          WRITE(*,*)' nmax=',NMAX,' unequal to ipmax=', IPMAX
          WRITE(*,*)'program stop here'
          DO I = 1,nmax
            DO I1 = 1,IPMAX
              IF(ABS(TZ(I)-XTMP(I1)) .LE. 0.001) THEN
                YTMP(I) = YTMP(I1)
                XTMP(I) = XTMP(I1)
                DIRP(I) = DIRP(I1)
                GOTO 165
              ENDIF
            ENDDO
165         CONTINUE
          ENDDO  
        ENDIF

        II = 0
        DO I = 1, NMAX
          IF((OP(I) .GT. -900.0) .AND. (YTMP(I) .GT. -900.0)) THEN
            II = II+1
            R(II) = OP(I)
            P(II) = YTMP(I)
            TP(II) = XTMP(I)
            E(II) = P(II)-R(II)
            IF(KINDAT .EQ. 1) THEN
              IF((OP(I) .GT. 0.26) .AND. (YTMP(I) .GT. 0.26)) THEN
                DIRR(II) = DIROP(I)
                DIRP(II) = DIRP(I)
                DIRE(II) = DIRP(II)-DIRR(II)
                IF(DIRE(II) .GT. 180.0) DIRE(II) = DIRE(II)-360.0
                IF(DIRE(II) .LT. -180.0) DIRE(II) = DIRE(II)+360.0
              ELSE 
                DIRR(II) = DIROP(I)
                DIRP(II) = DIRP(I)
                DIRE(II) = 0.0
              ENDIF
            ENDIF                
          ENDIF
        ENDDO
        NMAX = II

170     CONTINUE
        IF((KINDAT .EQ. 3) .OR. (KINDAT .EQ. 4)) GOTO 220

CC  Use continuous time series to get extrema
        IPR = 0
        IF(IS .EQ. 0) IPR = 1
        WRITE(*,"(/,' REFERENCE SERIES')")
        GAP = 1.5*DELT
        CALL CONTINUOUS(TP,NMAX,GAP,NSEGMENTS,IFRST,ILAST)
        WRITE(*,*) 'NSEGMENTS=', NSEGMENTS

C        RMDO = 0
        NMAXRH9 = 0
        NMAXRL9 = 0
        NMAXPH9 = 0
        NMAXPL9 = 0
        NSMAXR = 0
        NSMAXP = 0
        DO NG = 1, NSEGMENTS
          ISTART = IFRST(NG)
          IEND = ILAST(NG)
          DIF = TP(IEND)-TP(ISTART)
          WRITE(*,*)
          WRITE(*,*) 'START AND END TIME: ',ISTART,IEND
          IF(DIF .GT. 48.0) THEN
            WRITE(*,"(' SEGMENT: ',I3)") NG
            NUMB = IEND-ISTART+1  
            DO I = ISTART, IEND
              I0 = I-ISTART+1 
              XTMP(I0) = TP(I)
              YTMP(I0) = R(I)     ! For observation
            ENDDO
            CALL EXTREMES(XTMP,YTMP,NUMB,IPR1,DELHR,DELAMP,DELPCT,
     1        IOPTA,THIGHS0,HHIGHS0,IDX0,NSMAX0,TCUT,DELT,HHWP0,
     2        THWP0,NMAXRH0,HLWP0,TLWP0,NMAXRL0,NTYPE)

            DO J = 1, NMAXRH0
              NMAXRH9 = NMAXRH9+1
              HHWR(NMAXRH9) = HHWP0(J)
              THWR(NMAXRH9) = THWP0(J)
            ENDDO

            DO J = 1, NMAXRL0
              NMAXRL9 = NMAXRL9+1
              HLWR(NMAXRL9) = HLWP0(J)
              TLWR(NMAXRL9) = TLWP0(J)
            ENDDO

            DO J = 1, NSMAX0
              NSMAXR = NSMAXR+1
              THIGHSR(NSMAXR) = THIGHS0(J)
              HHIGHSR(NSMAXR) = HHIGHS0(J)
              IDXR(NSMAXR) = IDX0(J)
            ENDDO

CC  For prediction
            DO I = ISTART, IEND
              I0 = I-ISTART+1 
              YTMP(I0) = P(I)   ! For nowcast or forecast
            ENDDO
            CALL EXTREMES(XTMP,YTMP,NUMB,IPR1,DELHR,DELAMP,DELPCT,
     1        IOPTA,THIGHS0,HHIGHS0,IDX0,NSMAX0,TCUT,DELT,HHWP0,
     2        THWP0,NMAXPH0,HLWP0,TLWP0,NMAXPL0,NTYPE)

            DO J = 1, NMAXPH0
              NMAXPH9 = NMAXPH9+1
              HHWP(NMAXPH9) = HHWP0(J)
              THWP(NMAXPH9) = THWP0(J)
            ENDDO

            DO J = 1, NMAXPL0
              NMAXPL9 = NMAXPL9+1
              HLWP(NMAXPL9) = HLWP0(J)
              TLWP(NMAXPL9) = TLWP0(J)
            ENDDO

            DO J = 1, NSMAX0
              NSMAXP = NSMAXP+1
              THIGHSP(NSMAXP) = THIGHS0(J)
              HHIGHSP(NSMAXP) = HHIGHS0(J)
              IDXP(NSMAXP) = IDX0(J)
            ENDDO
          ENDIF
        ENDDO
        NMAXRH = NMAXRH9
        NMAXRL = NMAXRL9
        NMAXPH = NMAXPH9
        NMAXPL = NMAXPL9
        WRITE(*,*) 'nsmaxr=',NSMAXR,'nsmaxp=',NSMAXP

        IF(IS .EQ. 1) THEN
          OPEN(50,FILE=TRIM(FILESHORT)//'_timeseries_obsandtides.dat')
          OPEN(51,FILE=TRIM(FILESHORT)//'_extreme_predtides.dat')
          OPEN(52,FILE=TRIM(FILESHORT)//'_extreme_tides.dat')
          DO I = 1, NMAX
            WRITE(50,'(F13.8,2F10.4)') TP(I)/24.0,R(I),P(I)
          ENDDO
          DO I = 1, NSMAXR
            WRITE(51,'(F13.8,F10.4)') THIGHSR(I)/24.0,HHIGHSR(I)
          ENDDO
          DO I = 1, NSMAXP
            WRITE(52,'(F13.8,F10.4)') THIGHSP(I)/24.0,HHIGHSP(I)
          ENDDO
          CLOSE(50)
          CLOSE(51)
          CLOSE(52)
        ELSEIF(IS .EQ. 2) THEN
          OPEN(50,FILE=TRIM(FILESHORT)//'_timeseries_obsandhind.dat')
          OPEN(51,FILE=TRIM(FILESHORT)//'_extreme_obshind.dat')
          OPEN(52,FILE=TRIM(FILESHORT)//'_extreme_hind.dat')
          DO I = 1, NMAX
            WRITE(50,'(F13.8,2F10.4)') TP(I)/24.0,R(I),P(I)
          ENDDO
          DO I = 1, NSMAXR
            WRITE(51,'(F13.8,F10.4)') THIGHSR(I)/24.0,HHIGHSR(I)
          ENDDO
          DO I = 1, NSMAXP
            WRITE(52,'(F13.8,F10.4)') THIGHSP(I)/24.0,HHIGHSP(I)
          ENDDO
          CLOSE(50)
          CLOSE(51)
          CLOSE(52)
        ELSE IF(IS .EQ. 3) THEN
          OPEN(50,FILE=TRIM(FILESHORT)//'_timeseries_obsandnow.dat')
          OPEN(51,FILE=TRIM(FILESHORT)//'_extreme_obsnow.dat')
          OPEN(52,FILE=TRIM(FILESHORT)//'_extreme_now.dat')
          DO I = 1, NMAX
            WRITE(50,'(F13.8,2F10.4)') TP(I)/24.0,R(I),P(I)
          ENDDO
          DO I = 1, NSMAXR
            WRITE(51,'(F13.8,F10.4)') THIGHSR(I)/24.0,HHIGHSR(I)
          ENDDO
          DO I = 1, NSMAXP
            WRITE(52,'(F13.8,F10.4)') THIGHSP(I)/24.0,HHIGHSP(I)
          ENDDO
          CLOSE(50)
          CLOSE(51)
          CLOSE(52)
        ELSEIF(IS .EQ. 4) THEN
          OPEN(50,FILE=TRIM(FILESHORT)//'_timeseries_obsandfr.dat')
          OPEN(51,FILE=TRIM(FILESHORT)//'_extreme_obsfr.dat')
          OPEN(52,FILE=TRIM(FILESHORT)//'_extreme_fr.dat')
          DO I = 1, NMAX
            WRITE(50,'(F13.8,2F10.4)')TP(I)/24.0,R(I),P(I)
          ENDDO
          DO I = 1, NSMAXR
            WRITE(51,'(F13.8,F10.4)') THIGHSR(I)/24.0,HHIGHSR(I)
          ENDDO
          DO I = 1, NSMAXP
            WRITE(52,'(F13.8,F10.4)') THIGHSP(I)/24.0,HHIGHSP(I)
          ENDDO
          CLOSE(50)
          CLOSE(51)
          CLOSE(52)
        ENDIF

        IF(KINDAT .EQ. 1) THEN
CC  Derive current extrema time series of reference (observation)
          WRITE(*,*)
          WRITE(*,*) 'Calculating slack time of observation'
          CALL SLACK_TIME(TP,R,DIRR,THIGHSR,HHIGHSR,IDXR,NSMAXR,NMAX,
     1      DIRFLOOD,REAL(DELT),AFCR,AECR,TFCR,TECR,DFCR,DECR,TSFR,
     2      TEFR,TSER,TEER,NAFCR,NAECR,NTSFR,NTSER,NTSEER,IMX,NMX)

CC  Derive current extrema time series of prediction
          WRITE(*,*) 'Calculating slack time of model'
          CALL SLACK_TIME(TP,P,DIRP,THIGHSP,HHIGHSP,IDXP,NSMAXP,NMAX,
     1      DIRFLOOD,REAL(DELT),AFCP,AECP,TFCP,TECP,DFCP,DECP,TSFP,
     2      TEFP,TSEP,TEEP,NAFCP,NAECP,NTSFP,NTSEP,NTSEEP,IMX,NMX)

          GAP = 25.0
          IX1 = INT(X1*100+0.1)
          WRITE(TTL2,'(a1,I2,a9)')' ',IX1,' cm/s 24h'
          CALL EXTEMS(AFCR,TFCR,NAFCR,AFCP,TFCP,NAFCP,E,XTMP,NMAXB,1)
          SKILLV = SKILLA(AFCP,AFCR,NMX,NMAXB)
          CORR_C = CORRELATION(AFCP,AFCR,NMX,NMAXB)  
          WRITE(*,*) 'NMAXB 1:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('AFC-afc  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1        3,gap)
          END IF

          CALL EXTEMS(AECR,TECR,NAECR,AECP,TECP,NAECP,E,XTMP,NMAXB,1)
          SKILLV = SKILLA(AECP,AECR,NMX,NMAXB)
          CORR_C = CORRELATION(AECP,AECR,NMX,NMAXB)
          WRITE(*,*) 'NMAXB 2:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('AEC-aec  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1        3,gap)
          END IF

          WRITE(TTL2,'(a1,f4.2,a7)')' ',X2,'h   25h'
          CALL EXTEMS(AFCR,TFCR,NAFCR,AFCP,TFCP,NAFCP,E,XTMP,NMAXB,2)
          SKILLV = SKILLA(TFCP,TFCR,NMX,NMAXB)
          CORR_C = CORRELATION(TFCP,TFCR,NMX,NMAXB)

          WRITE(*,*) 'NMAXB 3:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('TFC-tfc  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1        3,gap)
          END IF

          CALL EXTEMS(AECR,TECR,NAECR,AECP,TECP,NAECP,E,XTMP,NMAXB,2)
          SKILLV = SKILLA(TECP,TECR,NMX,NMAXB)
          CORR_C = CORRELATION(TECP,TECR,NMX,NMAXB)

          WRITE(*,*) 'NMAXB 4:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('TEC-tec  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1        3,gap)
          END IF

          X3 = X2
          WRITE(TTL2,'(a1,f4.2,a7)')' ',X3/2,'h   25h'
          CALL EXTEMS(AFCR,TSFR,NTSFR,AFCP,TSFP,NTSFP,E,XTMP,NMAXB,2)
          WRITE(*,*) 'NMAXB 5:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('TSF-tsf  ',ttl2,lu,e,xtmp,imx,nmaxb,X3,
     1        3,gap)
          END IF

          CALL EXTEMS(AFCR,TEFR,NTSFR,AFCP,TEFP,NTSFP,E,XTMP,NMAXB,2)
          WRITE(*,*) 'NMAXB 6:',NMAXB
          IF(NMAXB .GT. 1) THEN
            call prtline1('TEF-tef  ',ttl2,lu,e,xtmp,imx,nmaxb,X3,
     1        3,gap)
          END IF

          WRITE(*,*) 'EXTEMS, SIZE(TSER): ',SIZE(TSER),' NTSER=',NTSEER
          CALL EXTEMS(AFCR,TSER,NTSEER,AFCP,TSEP,NTSEEP,E,XTMP,NMAXB,2)
          IF(NMAXB .GT. 1) THEN
            call prtline1('TSE-tse  ',ttl2,lu,e,xtmp,imx,nmaxb,X3,
     1        3,gap)
          END IF

          CALL EXTEMS(AFCR,TEER,NTSER,AFCP,TEEP,NTSEP,E,XTMP,NMAXB,2)
          IF(NMAXB .GT. 1) THEN
            call prtline1('TEE-tee  ',ttl2,lu,e,xtmp,imx,nmaxb,X3,
     1        3,gap)
          END IF

          WRITE(TTL2,'(a1,F4.1,a7)')' ',X11,' dg 24h'
          CALL EXTEMS(DFCR,TFCR,NAFCR,DFCP,TFCP,NAFCP,E,XTMP,NMAXB,1)
          DO I = 1,NMAXB
            IF(E(I) .GT. 180.0) E(I) = E(I)-360.0
            IF(E(I) .LT. -180.0) E(I) = E(I)+360.0
          ENDDO
          IF(NMAXB .GT. 1) THEN
            call prtline1('DFC-dfc  ',ttl2,lu1,e,xtmp,imx,nmaxb,X11,
     1        3,gap)
          END IF

          CALL EXTEMS(DECR,TECR,NAECR,DECP,TECP,NAECP,E,XTMP,NMAXB,1)
          DO I = 1, NMAXB
            IF(E(I) .GT. 180.0) E(I) = E(I)-360.0
            IF(E(I) .LT. -180.0) E(I) = E(I)+360.0
          ENDDO
          IF(NMAXB .GT. 1) THEN
            call prtline1('DEC-dec  ',ttl2,lu1,e,xtmp,imx,nmaxb,X11,
     1        3,gap)
          END IF
        ENDIF

        IF(KINDAT .EQ. 2) THEN
          IX1 = INT(X1*100+0.1)
          WRITE(TTL2,'(a3,I2,a7)')'   ',IX1,' cm 24h'
          GAP = 25.0
     
CC  Create extrema series for amplitudes
          CALL EXTEMS(HHWR,THWR,NMAXRH,HHWP,THWP,NMAXPH,E,XTMP,NMAXB,1)
          SKILLV = SKILLA(HHWP,HHWR,NMX,NMAXB)
          CORR_C = CORRELATION(HHWP,HHWR,NMX,NMAXB)
          IF((IS .EQ. 4) .OR. (IS .EQ. 5)) THEN
            IF(NMAXB .GT. 1) THEN
              call prtline1('AHW-ahw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          4,gap)
            ENDIF
          ELSE
            IF(NMAXB .GT. 1) THEN
              call prtline1('AHW-ahw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          3,gap)
            ENDIF
          ENDIF

          CALL EXTEMS(HLWR,TLWR,NMAXRL,HLWP,TLWP,NMAXPL,E,XTMP,NMAXB,1)
          SKILLV = SKILLA(HLWP,HLWR,NMX,NMAXB)
          CORR_C = CORRELATION(HLWP,HLWR,NMX,NMAXB)

          IF((IS .EQ. 4) .OR. (IS .EQ. 5)) THEN
            IF(NMAXB .GT. 1) THEN
              call prtline1('ALW-alw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          4,gap)
            ENDIF
          ELSE
            IF(NMAXB .GT. 1) THEN
              call prtline1('ALW-alw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          3,gap)
            ENDIF
          ENDIF

CC  Create extrema series for times 
          TTL2 = '    .5h  25h'
          WRITE(TTL2,'(a1,f4.2,a7)')' ',X2,'  h 25h'
          CALL EXTEMS(HHWR,THWR,NMAXRH,HHWP,THWP,NMAXPH,E,XTMP,NMAXB,2)
          SKILLV = SKILLA(THWP,THWR,NMX,NMAXB)
          CORR_C = CORRELATION(THWP,THWR,NMX,NMAXB)
          IF((IS .EQ. 4) .OR. (IS .EQ. 5)) THEN
            IF(NMAXB .GT. 1) THEN
              call prtline1('THW-thw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          4,gap)
            ENDIF
          ELSE
            IF(NMAXB .GT. 1) THEN
              call prtline1('THW-thw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          3,gap)
            ENDIF
          ENDIF

          CALL EXTEMS(HLWR,TLWR,NMAXRL,HLWP,TLWP,NMAXPL,E,XTMP,NMAXB,2)
          SKILLV = SKILLA(TLWP,TLWR,NMX,NMAXB)
          CORR_C = CORRELATION(TLWP,TLWR,NMX,NMAXB)

          IF((IS .EQ. 4) .OR. (IS .EQ. 5)) THEN
            IF(NMAXB .GT. 1) THEN
              call prtline1('TLW-tlw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          4,gap)
            ENDIF
          ELSE  
            IF(NMAXB .GT. 1) THEN
              call prtline1('TLW-tlw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          3,gap)
            ENDIF
          ENDIF
        ENDIF
220   CONTINUE

      IF(ISWITCH(2) .GT. 0 .AND. ISURGE .GT. 0 .AND.
     1   KINDAT .EQ. 2) THEN
        WRITE(LU,'(5x,a36)') 'SCENARIO: HINDCAST (SURGE ONLY)     '
        NMAX = IMAXB
        II = 0
        DO I = 1,NMAX
          IF((SGO(I) .GT. -900.0) .AND. (SGM(I) .GT. -900.0)) THEN
            II = II+1
            R(II) = SGO(I)
            P(II) = SGM(I)
            E(II) = P(II)-R(II)
            TP(II) = T(I)
            YTMP(II) = 0.0
            DIRR(II) = DIRO(I)
            DIRP(II) = DIRHI(I)
            DIRE(II) = 0.0
          ENDIF

          IF(I .LE. 10) THEN
            WRITE(*,"('  is=2  i=',i5,' t=',f11.4,' p=',f9.4,
     1         ' r=',f9.4,'  e=',f9.4)")i,t(i),p(i),r(i),e(i)
          ENDIF
        ENDDO
        NMAX = II
        WCOF = WOF(P,YTMP,R,IMX,NMAX,2.0*X1)

C  Test: sum of cf, pof, nof = 1
        WRITE(*,"(' test: nof(x),cf(x),pof(x)=',3f8.3,' sum=',f8.3)")
     1    VOF(E,IMX,NMAX,X1),CF(E,IMX,NMAX,X1),POF(E,IMX,NMAX,X1),
     2    VOF(E,IMX,NMAX,X1)+CF(E,IMX,NMAX,X1)+POF(E,IMX,NMAX,X1)
        TTL2 = '         '
        GAP = 1.5*DELT
        CALL PRTLINE1('H      ',TTL2,LU,P,TP,IMX,NMAX,X1,1,GAP)
        CALL PRTLINE1('H      ',TTL2,LU,R,TP,IMX,NMAX,X1,1,GAP)
        IX1 = INT(X1*100+0.1)
        WRITE(TTL2,'(A3,I2,A7)')'   ',IX1,' CM 24H'
        SKILLV = SKILLA(P,R,IMX,NMAX)
        CORR_C = CORRELATION(P,R,IMX,NMAX)
        CALL PRTLINE1('H-H    ',TTL2,LU,E,TP,IMX,NMAX,X1,5,GAP)

        OPEN(50,FILE=TRIM(FILESHORT)//'_timeseries_surgeOM.dat')
        DO I = 1,NMAX
          WRITE(50,'(3F10.4)') TP(I)/24.0,R(I),P(I)
        ENDDO
        CLOSE(50)

C  Pick max/min in time series
        HHIGHP = -900.0
        HLOWP = 900.0
        HHIGHR = -900.0
        HLOWR = 900.0
        DO J = 1,NMAX
          IF(P(J) .GT. HHIGHP) THEN
            HHIGHP = P(J)
            THIGHP = TP(J)
          ENDIF

          IF(P(J) .LT. HLOWP) THEN
            HLOWP = P(J)
            TLOWP = TP(J)
          ENDIF

          IF(R(J) .GT. HHIGHR) THEN
            HHIGHR = R(J)
            THIGHR = TP(J)
          ENDIF

          IF(R(J) .LT. HLOWR) THEN
            HLOWR = R(J)
            TLOWR = TP(J)
          ENDIF
        ENDDO

        WRITE(LU,"('VARIABLE         Obs             Model   ')")
        WRITE(LU,131)'HWM', hhighr, thighr/24., hhighp,thighp/24.
        WRITE(LU,131)'LWM', hlowr, tlowr/24., hlowp,tlowp/24.
      ENDIF
131   FORMAT(A7,4(1X,F8.3))

      CLOSE(LU)
      IF(KINDAT .EQ. 1) THEN
        CLOSE(LU1)
      END IF

      RETURN
      END

