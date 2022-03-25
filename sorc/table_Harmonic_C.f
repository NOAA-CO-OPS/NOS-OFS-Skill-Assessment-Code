      PROGRAM Table_Harmonic_C
      PARAMETER(NCON=37)
      DIMENSION PHAM2(99,NCON),AMPM2(99,NCON),FLOODDIR(99)
      DIMENSION AMP_O(NCON),EPOC_O(NCON),AMP_M(NCON),EPOC_M(NCON)
      DIMENSION ISTA(6),NO(6)

      CHARACTER*20 CC,CAZI_O,CAZI_M,BUFFER*100
      CHARACTER*10 ALIST(37),LABEL,LABLE,RATIO_O*9,RATIO_M*9
      CHARACTER*100 HEAD_O(2),HEAD_M(2)       
      CHARACTER*200 FILE_OBS,FILE_MODEL,FNAME,FOUT,LONGNAME*40
      CHARACTER*200 DATADIR
      DATA (ALIST(I),I=1,37) /'M(2)      ','S(2)      ','N(2)      ',
     1                        'K(1)      ','M(4)      ','O(1)      ',
     2                        'M(6)      ','MK(3)     ','S(4)      ',
     3                        'MN(4)     ','NU(2)     ','S(6)      ',
     4                        'MU(2)     ','2N(2)     ','OO(1)     ',
     5                        'LAMBDA(2)  ','S(1)      ','M(1)      ',
     6                        'J(1)      ','MM        ','SSA       ',
     7                        'SA        ','MSF       ','MF        ',
     8                        'RHO(1)    ','Q(1)      ','T(2)      ',
     9                        'R(2)      ','2Q(1)     ','P(1)      ',
     1                        '2SM(2)    ','M(3)      ','L(2)      ',
     2                        '2MK(3)   ','K(2)      ','M(8)      ',
     3                        'MS(4)     '/

      TM = 0.0
      CALL GETARG(1,BUFFER)
      READ(BUFFER,*) KINDAT
      CALL GETARG(2,FNAME)
      CALL GETARG(3,LONGNAME)
      CALL GETARG(4,DATADIR)
      FILE_OBS = '../data/harmonic_con/'//trim(fname)//'.std'
      FILE_OBS = trim(DATADIR)//'/harmonic_con/'//trim(fname)//'.std'
      FILE_MODEL = trim(fname)//'_modeltides.std'
      FOUT = TRIM(FNAME)//'table_HC_comparison.out'
      OPEN(8,FILE=FILE_OBS) 
      OPEN(9,FILE=FILE_MODEL) 
      OPEN(10,FILE=FOUT)
      NOS0 = 0
10    READ(8,550) HEAD_O(1)
      READ(8,550) HEAD_O(2)
      READ(8,532) DATUM,ISTA(1),NO(1),(AMP_O(J),EPOC_O(J),J=1,7),
     1  ISTA(2),NO(2),(AMP_O(J),EPOC_O(J),J=8,14),
     2  ISTA(3),NO(3),(AMP_O(J),EPOC_O(J),J=15,21),
     3  ISTA(4),NO(4),(AMP_O(J),EPOC_O(J),J=22,28),
     4  ISTA(5),NO(5),(AMP_O(J),EPOC_O(J),J=29,35),
     5  ISTA(6),NO(6),(AMP_O(J),EPOC_O(J),J=36,37)
550   FORMAT(A100)
532   FORMAT(F6.3,6(/2I4,7(F5.3,F4.1)))
      IND = INDEX(HEAD_O(2),'along')+5
      CAZI_O=HEAD_O(2)(IND:IND+4)
      IND = INDEX(HEAD_O(2),'R=')
      RATIO_O = HEAD_O(1)(IND:IND+8)
      Ind = INDEX(HEAD_O(2),'Hour')+9
      HEAD_O(2) = 'Observation: '//HEAD_O(2)(1:IND)
      IF(KINDAT .EQ. 2) THEN
        HEAD_O(2) = 'Observation: CO-OPS Accepted Harmonic Constants'
      END IF

      READ(9,550)  HEAD_M(1),HEAD_M(2)
      READ(9,532)DATUM,ISTA(1),NO(1),(AMP_M(J),EPOC_M(J),J=1,7),
     1  ISTA(2),NO(2),(AMP_M(J),EPOC_M(J),J=8,14),
     2  ISTA(3),NO(3),(AMP_M(J),EPOC_M(J),J=15,21),
     3  ISTA(4),NO(4),(AMP_M(J),EPOC_M(J),J=22,28),
     4  ISTA(5),NO(5),(AMP_M(J),EPOC_M(J),J=29,35),
     5  ISTA(6),NO(6),(AMP_M(J),EPOC_M(J),J=36,37)
      IND = INDEX(HEAD_M(2),'along')+5
      CAZI_M = HEAD_M(2)(IND:IND+4)
      IND = INDEX(HEAD_O(2),'R=')
      RATIO_M = HEAD_M(1)(IND:IND+8)
      IND = INDEX(HEAD_M(2),'Hour')+9
      HEAD_M(2) = 'Model: '//HEAD_M(2)(1:IND)
      IF (NOS0 .EQ. 0) THEN
      WRITE(10,"(/,'Station: ',A40)") LONGNAME
      WRITE(10,"(A100)") HEAD_O(2)
      WRITE(10,"(A100)") HEAD_M(2)

      IF(KINDAT .EQ. 1) THEN
        IF(TM .EQ. 0.0) THEN
          WRITE(10,*) 
     1      'Amplitudes are in m/s, and Phase is in degrees (GMT)'
        ELSE IF(TM .EQ. 75.0) THEN
          WRITE(10,*) 
     1      'Amplitudes are in m/s, and Phase is in degrees (EST)'
        END IF

      ELSE IF(KINDAT .EQ. 2) THEN
        IF(TM .EQ. 0.0) THEN
          WRITE(10,*)
     1      'Amplitudes are in meters, and Phase is in degrees (GMT)'
        ELSE IF(TM .EQ. 75.0) THEN
          WRITE(10,*)
     1      'Amplitudes are in meters, and Phase is in degrees (EST)'
        END IF
      ENDIF
      WRITE(10,"('------------------------------------------',
     1 '-----------------------------------------')")

      IF(KINDAT .EQ. 1) THEN
        WRITE(10,171) RATIO_O, RATIO_M
      ELSE IF(KINDAT .EQ. 2) THEN
        WRITE(10,"('                            Observed            ',
     1 '   Modeled             Difference')")
      ENDIF
171   FORMAT(22x,9HObserved(,a9,1H),3x,8HModeled(,a9,1H),7x
     1   ,10HDifference)

      WRITE(10,"('    N   Constituent   Amplitude    Epoch    ',
     1 'Amplitude    Epoch    Amplitude    Epoch')")
      write(10,"('------------------------------------------',
     1 '-----------------------------------------',/)")
      ENDIF
      IF((KINDAT .EQ. 1) .AND. (NOS0 .EQ. 0)) THEN
        WRITE(10,"('CURRENT ALONG PCD',10x,4HDIR=,A4,14x,4HDIR=,A4)")
     1          CAZI_O,CAZI_M 
      ELSE IF( (KINDAT .EQ. 1) .and. (NOS0 .eq. 1) )THEN
       WRITE(10,"(/'CURRENT ACROSS PCD ',10x,4HDIR=,A4,14x,4HDIR=,A4)")
     1          CAZI_O,CAZI_M 
      ENDIF

      DO J = 1,NCON
        IF(AMP_O(J) .le. 0.00001) THEN
          WRITE(10,667) J,ALIST(J),AMP_O(J),EPOC_O(J),AMP_M(J)
     1    ,EPOC_M(J),0.0,0.0
        ELSE 
          DIFF = EPOC_M(J)-EPOC_O(J)
          IF(DIFF .LT. -180.0) DIFF = DIFF+360.0
          IF(DIFF .GT.  180.0) DIFF = 360.0-DIFF
          WRITE(10,667)J,ALIST(J),AMP_O(J),EPOC_O(J),AMP_M(J),
     1      EPOC_M(J),AMP_M(J)-AMP_O(J),DIFF
        endif
      enddo
      IF((KINDAT .EQ. 1) .AND. (NOS0 .eq. 0)) THEN
        NOS0 = 1
        GOTO 10
      ENDIF
        
667   FORMAT(I5,5x,a10,2x,F7.3,4x,F6.1,5x,F7.3,4x,F6.1,5x,F7.3,4x,F6.1)
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)

      STOP
      END

