C      calculate maximum continuous time of a time series (gap < 1.5 delt)
C      lf95 gaps_zaj.f utility.f -o gaps_zaj.x
      parameter (NT=200000)
      character BUFFER*100,FNAME*100 
      real*8 jdayb,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      DIMENSION TIME(NT),VIN(NT),DIN(NT),SIYR(NT),SMON(NT)
     1          ,SDD(NT),SHH(NT),SMIN(NT),ifrst(NT),ilast(NT)
     2          ,XTMP(NT),YTMP(NT)
      CALL GETARG(1,FNAME)
      CALL GETARG(2,BUFFER)
      READ(BUFFER,*)KINDAT
      CALL GETARG(3,BUFFER)
      READ(BUFFER,*)DELT 
      DELT=DELT/60.0  !! 11/30/2006  convert from minutes into hours

      OPEN(2,file=trim(FNAME))
      OPEN(11,file='tmp.dat')
      OPEN(12,file='tmp1.dat')
      OPEN(13,file='tmp2.dat')
      N=0
      IF(KINDAT .EQ. 1)THEN
         I=1
5        READ(2,*,END=10)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I),DIN(I)       
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
        II=0
        DO I=1,NMAX
          IF ((VIN(i) .GT. -900.) .and. (DIN(i) .GT. -900.) )then
            II=II+1
            TIME(II)=TIME(I)
            SIYR(II)=SIYR(I)
            SMON(II)=SMON(I)
             SDD(II)=SDD(I)
             SHH(II)=SHH(I)
            SMIN(II)=SMIN(I)
             VIN(II)=VIN(I)
             DIN(II)=DIN(I)
          ENDIF
        ENDDO
        NMAX=II
      ELSE IF(KINDAT .EQ. 2)THEN
         I=1
15       READ(2,*,END=20)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I)
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
        II=0
        DO I=1,NMAX
          IF (VIN(i) .GT. -900. )then
            II=II+1
            TIME(II)=TIME(I)
            SIYR(II)=SIYR(I)
            SMON(II)=SMON(I)
             SDD(II)=SDD(I)
             SHH(II)=SHH(I)
            SMIN(II)=SMIN(I)
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
      ENDDO
      NMAX=ilast(INDEX)-ifrst(INDEX)+1
      max_time=int(NMAX*DELT/24)
      print *,'longest continuous time=',max_time
      write(11,'(I6)')max_time
      ISTART=ifrst(INDEX)
      IYRS=SIYR(Istart)
      IMMS=SMON(Istart)
      IDDS=SDD(Istart)
      IHHS=SHH(Istart)
      IMNS=SMIN(Istart)
      IEND=ilast(INDEX)
      IYRE=SIYR(IEND)
      IMME=SMON(IEND)
      IDDE=SDD(IEND)
      IHHE=SHH(IEND)
      IMNE=SMIN(IEND)
      write(12,100)IYRS,IMMS,IDDS,IHHS,IMNS
      write(13,100)IYRE,IMME,IDDE,IHHE,IMNE
      close(2)
      close(11)
      close(12)
      close(13)
100   FORMAT(I4,4(1x,I2))
      stop
      end
