      character*120 fileshort,station_file,staid,stationname
      character*120 filename,CTMP,CTMP1,FMT,filein
      character*120 FVCOM, ROMS,SELFE,POM,FCTL, BUFFER
      REAL RMSE(20),SM(20),NOF(20),CF(20),POF(20),SD(20)
      LOGICAL FEXIST
      print *,'input KINDAT.. '
      READ(5,*)KINDAT
      IF(KINDAT .EQ. 1)FCTL='hsofs_cu_station.ctl'
      IF(KINDAT .EQ. 2)FCTL='hsofs_wl_station.ctl'
      IF(KINDAT .EQ. 3)FCTL='hsofs_temp_station.ctl'
      IF(KINDAT .EQ. 4)FCTL='hsofs_salt_station.ctl'
      IF(KINDAT .EQ. 5)FCTL='hsofs_cu_station.ctl'

      IF(KINDAT .EQ. 1)OPEN(16,file='CURRENT_COMPARISON.out')
      IF(KINDAT .EQ. 2)OPEN(16,file='WL_COMPARISON.out')
      IF(KINDAT .EQ. 3)OPEN(16,file='TEMP_COMPARISON.out')
      IF(KINDAT .EQ. 4)OPEN(16,file='SALT_COMPARISON.out')
      IF(KINDAT .EQ. 5)OPEN(16,file='CURRENTPhase_COMPARISON.out')
      BUFFER='STATION_ID  STA_NAME               NOWCAST  FH00   FH06'
      BUFFER=trim(BUFFER)//'   FH12   FH18   FH24   FH30   FH36'
      BUFFER=trim(BUFFER)//'   FH42   FH48'
      WRITE(16,*)trim(BUFFER)
      open(1,file='../control_files/'//trim(FCTL))
CC  generate table for each station, loop through all stations      
10    read(1,*,err=999,end=999)staid,fileshort,stationname
  !   write(6,*)'fileshort=',fileshort
      write(6,*)'stationname=',trim(stationname)
      read(1,*)alat,alon,dirflood,sdepth
      ROMS='../work/'//trim(fileshort)//'_table.out'
      write(*,*)trim(ROMS)
      IF(KINDAT .EQ. 5)THEN
        ROMS='../work/'//trim(fileshort)//'phase_table.out'
      ENDIF
      INQUIRE(FILE=trim(ROMS),EXIST=FEXIST)
      IF(.NOT.FEXIST)THEN
           write(6,*)'filein=',trim(ROMS), ' does not exist'
           goto 10
      ENDIF
      close(10)
      open(10,file=trim(ROMS),form='formatted',status='old')
      DO I=1,13
         read(10,100,err=10,end=10)CTMP
         if(I .EQ. 12)THEN
          IF(LEN(TRIM(CTMP(1:10))) .EQ. 0)GOTO 10
         ENDIF
      ENDDO
      read(10,100,err=10)ROMS
      read(ROMS,130)SM(1),RMSE(1),SD(1),NOF(1),CF(1),POF(1)
      write(*,*)trim(ROMS)

      DO I=1,5
         read(10,100,err=10,end=10)CTMP
      ENDDO
      DO N=1,9 !NLOOP
        N1=N+1
        read(10,100,err=10)ROMS
        read(ROMS,130)SM(N1),RMSE(N1),SD(N1),NOF(N1),CF(N1),POF(N1)
      ENDDO
      IF(KINDAT .EQ. 2)THEN
         DO N=1,NLOOP+1
           SM(N)=SM(N)*100.0
	   RMSE(N)=RMSE(N)*100.
	   SD(N)=SD(N)*100.0
         ENDDO
      ENDIF	 
      LL=len(trim(stationname))
      if(LL .GT. 25)LL=25
      IF(KINDAT .EQ. 2)THEN
        write(FMT,'(a6,a1,I2.2,a1,I2.2,a1,a12)')'A8,2x,',
     1'A',LL,',',25-LL,'X',',10(1x,f6.2)'
        FMT="("//trim(FMT)//")"
        print *,trim(FMT)
        write(16,FMT)trim(staid),trim(stationname),(RMSE(I),I=1,10)
!      STOP	
      ELSE        
        write(FMT,'(a1,I2.2,a1,I2.2,a1,a40)')
     1'A',LL,',',12-LL,'X',',f5.1,3x,f5.1,5x,4f7.3,2x,4f6.1,2x,4f6.2'
        FMT="("//trim(FMT)//")"
        write(16,FMT)trim(stationname),tdepth,sdepth,RMS1,RMS2,RMS3,
     1    RMS4,CF1,CF2,CF3,CF4,SKILL1,SKILL2,SKILL3,SKILL4
      ENDIF
        goto 10
100     format(a120)
 130  format(27x,3f7.3,3f6.1,18x,f6.2) 
! 130  format(27x,3f7.3,3f6.1,1x,f6.1,f5.1,f6.2) 
! 130  format(a7,1x,a12,i6,1x,3f7.3,3f6.1,1x,f6.1,f5.1,f6.2) 
999     stop
        end

        
