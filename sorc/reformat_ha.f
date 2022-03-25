CC  rewrite the tidal harmonic constants from CO-OPS website as the same format of 
C   lsqha.f and harm29d.f output 'cons.out', and used by pref.f. The 37 tidal constituents are in 
C   the following order: 
C   name         speed
C  1 M2         28.9841042
C  2 S2         30.0000000
C  3 N2         28.4397295
C  4 K1         15.0410686
C  5 M4         57.9682084
C  6 O1         13.9430356
C  7 M6         86.9523127
C  8 MK3        44.0251729
C  9 S4         60.0000000
C 10 MN4        57.4238337
C 11 NU2        28.5125831
C 12 S6         90.0000000
C 13 MU2        27.9682084
C 14 2N2        27.8953548
C 15 OO1        16.1391017
C 16 LAM2       29.4556253
C 17 S1         15.0000000
C 18 M1         14.4966939
C 19 J1         15.5854433
C 20 MM       	0.5443747
C 21 SSA      	0.0821373
C 22 SA       	0.0410686
C 23 MSF      	1.0158958
C 24 MF       	1.0980331
C 25 RHO        13.4715145
C 26 Q1         13.3986609
C 27 T2         29.9589333
C 28 R2         30.0410667
C 29 2Q1        12.8542862
C 30 P1         14.9589314
C 31 2SM2       31.0158958
C 32 M3         43.4761563
C 33 L2         29.5284789
C 34 2MK3       42.9271398
C 35 K2         30.0821373
C 36 M8        115.9364166
C 37 MS4        58.9841042
C  Modified by Lianyuan Zheng on 03/01/2017

      PROGRAM REFORMAT_HA
      CHARACTER*100 FNAME, BUFFER, CTMP1, CTMP2
      CHARACTER*80 FILEI, FILEO, HEADER1, HEADER2, CC, CCC
      CHARACTER*1  C1,C2
      INTEGER ID(37), IW(37)
      REAL*4 AMP(37), PHA(37)
      REAL*4 MHHW, MHW, MTL, MSL, MLW, MLLW
      LOGICAL FEXIST, FOPEN

      CALL GETARG(1, FILEI)

      DATUM = 0.0
      MHHW  = 0.0
      MHW   = 0.0
      MTL   = 0.0
      MSL   = 0.0
      MLW   = 0.0
      MLLW  = 0.0
      FNAME = './'//TRIM(FILEI)//'.datum'
      INQUIRE(FILE = TRIM(FNAME), EXIST = FEXIST)
      IF(.NOT. FEXIST) THEN
        WRITE(*,*) ' File ',TRIM(FNAME), ' does not exist, stop!'
        STOP
      ELSE
        OPEN(24, FILE = TRIM(FNAME), FORM = 'FORMATTED')
60      READ(24,'(A100)',END = 64) BUFFER
        BUFFER = TRIM(ADJUSTL(BUFFER))

        LL = INDEX(BUFFER,'No station was found')
        IF(LL .GT. 0) THEN
          WRITE(*,*) 'No station ',TRIM(FILEI), ' was found'
	  STOP
        ENDIF

        LL = INDEX(BUFFER,'Wrong Station ID')
        IF(LL .GT. 0) THEN
          WRITE(*,*) 'No data at station ',TRIM(FILEI), ' was found'
	  STOP
        ENDIF

        LL = INDEX(BUFFER,'but no datums data is available')
        IF(LL .GT. 0) THEN
          WRITE(*,*) 'No data at station ',TRIM(FILEI), ' was found'
	  STOP
        ENDIF

        LL0 = INDEX(BUFFER,'MHHW')
        LL1 = INDEX(BUFFER,'Mean Higher-High Water')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MHHW, CTMP2
	ENDIF  
        LL0 = INDEX(BUFFER,'MHW')
        LL1 = INDEX(BUFFER,'Mean High Water')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MHW, CTMP2
	ENDIF  
        LL0 = INDEX(BUFFER,'MTL')
        LL1 = INDEX(BUFFER,'Mean Tide Level')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MTL, CTMP2
	ENDIF  
        LL0 = INDEX(BUFFER,'MSL')
        LL1 = INDEX(BUFFER,'Mean Sea Level')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MSL, CTMP2
	ENDIF  
        LL0 = INDEX(BUFFER,'MLW')
        LL1 = INDEX(BUFFER,'Mean Low Water')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MLW, CTMP2
	ENDIF  
        LL0 = INDEX(BUFFER,'MLLW')
        LL1 = INDEX(BUFFER,'Mean Lower-Low Water')
        IF((LL0 .GT. 0) .AND. (LL1 .GT. 0)) THEN
          READ(BUFFER,*) CTMP1, MLLW, CTMP2
	ENDIF  
        DATUM = MSL- MLLW
	GOTO 60 
      ENDIF	
64    CLOSE(24) 
      WRITE(*,*) 'datum = ',DATUM,MHHW,MHW,MTL,MSL,MLW,MLLW

      FNAME = './'//TRIM(FILEI)//'.ha'
      INQUIRE(FILE = TRIM(FNAME), EXIST = FEXIST)
      IF(.NOT. FEXIST) GOTO 99
      OPEN(25,FILE = TRIM(FNAME), FORM = 'FORMATTED')
      DO I = 1, 30
        READ(25,'(A80)',END = 98) CC
      ENDDO  
      NL = 30
5     READ(25,'(A80)',END = 98) CC
      LL = INDEX(CC,'No Data Exists')
      IF(LL .GT. 0) THEN
        WRITE(*,*) 'No Data Exists at station: ',TRIM(FILEI)
	STOP
      ENDIF

      LL = INDEX(CC,'No station was found')
      IF(LL .GT. 0) THEN
        WRITE(*,*) 'No station ',TRIM(FILEI), ' was found'
	STOP
      ENDIF

      LL0 = INDEX(CC,'Amplitude')
      LL1 = INDEX(CC,'Phase')
      LL2 = INDEX(CC,'Speed')
      NL = NL + 1
      IF (LL0 .GT. 0 .AND. LL1 .GT. 0 .AND. LL2 .GT. 0) THEN
        NL1 = NL + 1
      END IF
      LL0 = INDEX(CC,'</pre>')
      IF (LL0 .LE. 0) GOTO 5
      NL2 = NL - 1
      REWIND(25)

      OPEN(10, FILE = TRIM(FILEI)//'.std', status = 'unknown')
      OPEN(11, FILE = TRIM(FILEI)//'.cons', status = 'unknown')

      DO I = 1, NL1
        READ(25,'(A80)',END = 97) CC
        LL = INDEX(CC,'Latitude')
	IF(LL .GT. 0) THEN
	  READ(CC,*) BUFFER, BUFFER, ALT
	ENDIF  

        LL = INDEX(CC,'Longitude')
	IF(LL .GT. 0) THEN
	  READ(CC,*) BUFFER, BUFFER, ALO
	ENDIF
      ENDDO

      WRITE(HEADER2,'(A11,A7,A11,F8.4,A12,F9.4)') 'StationID: ',
     1	TRIM(FILEI), ' Latitude: ', ALT, ' Longitude: ', ALO
      
      DO J = 1, 37
        AMP(J) = 0.0
        PHA(J) = 0.0
      ENDDO
      DO J = 1, NL2 - NL1
        AMP(J) = 0.0
        PHA(J) = 0.0
        READ(25,*,ERR = 1000) ITMP, BUFFER, AMP(J), PHA(J)
      ENDDO
      CLOSE(25)

      HEADER1 = 'Amplitudes are in Meters, Phases are in degrees,'
      HEADER1 = TRIM(HEADER1)//' referenced to UTC (GMT)'
      WRITE(10,'(A80)') HEADER2
      WRITE(10,'(A80)') HEADER1
      WRITE(11,'(A80)') HEADER2
      WRITE(11,'(A80)') HEADER1

      IAVE = (DATUM + 0.0005) * 1000
      WRITE(11,'(I6)') IAVE
      WRITE(10,'(I6)') 0
      DO 15 I = 1, 37                                                               
        IW(I) = NINT(PHA(I)*10.0)                                             
        ID(I) = NINT((AMP(I))*1000.0)              
15    CONTINUE                                                                  

      DO 11 N = 1, 5
        NN = 7*(N-1)
        WRITE(10,6) N, (ID(NN+J),IW(NN+J),J = 1, 7)
        WRITE(11,6) N, (ID(NN+J),IW(NN+J),J = 1, 7)
11    CONTINUE

      WRITE(10,6) 6, ID(36), IW(36), ID(37), IW(37)
      WRITE(11,6) 6, ID(36), IW(36), ID(37), IW(37)
6     FORMAT(7X,I1,7(I5,I4))
      STOP

99    PRINT*, ' File ',TRIM(FNAME), ' does not exist, stop!'
      STOP
98    CONTINUE
      CLOSE(25)
      STOP
97    CONTINUE
      CLOSE(25)
      CLOSE(10)
      CLOSE(11)
      STOP

1000  WRITE(*,*) 'THERE IS ERR AT STATION:',TRIM(FILEI)
      CLOSE(25)
      CLOSE(10)
      CLOSE(11)

      STOP
      END

