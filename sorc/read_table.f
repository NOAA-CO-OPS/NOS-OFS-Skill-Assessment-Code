      character*120 fileshort,station_file,staid,stationname
      character*120 filename,CTMP,CTMP1,FMT,filein
      character*120 FVCOM, ROMS,SELFE,POM,FCTL
      LOGICAL FEXIST
      print *,'input KINDAT.. '
      READ(5,*)KINDAT
      IF(KINDAT .EQ. 1)FCTL='cbofs_cu_station.ctl'
      IF(KINDAT .EQ. 2)FCTL='tbofs_wl_station.ctl'
      IF(KINDAT .EQ. 3)FCTL='DEL_temp_stations.ctl'
      IF(KINDAT .EQ. 4)FCTL='DEL_salt_stations.ctl'

      IF(KINDAT .EQ. 1)OPEN(16,file='cbofs_CURRENT_COMPARISON.out')
      IF(KINDAT .EQ. 2)OPEN(16,file='tbofs_WL_COMARISON.out')
      IF(KINDAT .EQ. 3)OPEN(16,file='cbofs_TEMP_COMPARISON.out')
      IF(KINDAT .EQ. 4)OPEN(16,file='cbofs_SALT_COMPARISON.out')
!      WRITE(16,200)'RMSE','CF','SKILL'
      WRITE(16,210)'STATION','SM','RMSE','NOF','CF','POF','SKILL'
200   FORMAT(44x,a4,23x,a2,23x,a5)
210   FORMAT(a7,25x,a4,6x,a4,3x,a5,3x,a2,3x,a4,3x,a5,2x,a5,1x,a5
     1  ,4x,a4,2x,a3,2x,a5,1x,a5)
      open(1,file='../control_files/'//trim(FCTL))
CC  generate table for each station, loop through all stations      
10      read(1,*,err=999,end=999)staid,fileshort,stationname
  !      write(6,*)'fileshort=',fileshort
        write(6,*)'stationname=',stationname
        read(1,*)alat,alon,dirflood,sdepth
        ROMS='../work/'//trim(fileshort)//'_table.out'
        FVCOM='../work/'//trim(fileshort)//'_table.out'
        POM='../work_POM/'//trim(fileshort)//'_table.out'
        SELFE='../work_SELFE/'//trim(fileshort)//'_table.out'
        INQUIRE(FILE=trim(ROMS),EXIST=FEXIST)
        IF(.NOT.FEXIST)THEN
           write(6,*)'filein=',trim(fileshort), ' does not exist'
           goto 10
        ENDIF
        close(10)
        close(11)
        close(12)
        close(13)
        open(10,file=trim(ROMS),form='formatted',status='old')
 !       open(11,file=trim(POM),form='formatted',status='old')
 !       open(12,file=trim(FVCOM),form='formatted',status='old')
 !       open(13,file=trim(SELFE),form='formatted',status='old')
        DO I=1,13
         read(10,100,err=10,end=10)CTMP
  !       read(11,100,err=10,end=10)CTMP
   !      read(12,100,err=10,end=10)CTMP
   !      read(13,100,err=10,end=10)CTMP
         if(I .EQ. 12)THEN
          IF(LEN(TRIM(CTMP(1:10))) .EQ. 0)GOTO 10
         ENDIF
        ENDDO

        read(10,100,err=10)ROMS
  !      read(11,100,err=10)POM
  !      read(12,100,err=10)FVCOM
  !      read(13,100,err=10)SELFE
        read(ROMS,130)SM1,RMS1,SD1,ANOF1,CF1,POF1,SKILL1
  !      read(POM,130)SM2,RMS2,SD2,ANOF2,CF2,POF2,SKILL2
  !      read(FVCOM,130)SM3,RMS3,SD3,ANOF3,CF3,POF3,SKILL3
  !      read(SELFE,130)SM4,RMS4,SD4,ANOF4,CF4,POF4,SKILL4
        LL=len(trim(stationname))
	IF(LL .GT. 30)LL=29
        write(FMT,'(a1,I2.2,a1,I2.2,a1,a31)')
     1'A',LL,',',30-LL,'X',',f7.3,2x,f7.3,2x,3f6.1,2x,4f6.2'
        FMT="("//trim(FMT)//")"
        print *,trim(FMT)
        write(16,FMT)trim(stationname),SM1,RMS1,ANOF1,CF1,POF1,SKILL1

!        write(16,FMT)trim(stationname),RMS1,RMS2,SD1,SD2
!     1    ,ANOF1,ANOF2,CF1,CF2,POF1,POF2
        
        goto 10
100     format(a120)
 130  format(27x,3f7.3,3f6.1,18x,f6.2) 
! 130  format(27x,3f7.3,3f6.1,1x,f6.1,f5.1,f6.2) 
! 130  format(a7,1x,a12,i6,1x,3f7.3,3f6.1,1x,f6.1,f5.1,f6.2) 
800     FORMAT(32x,f6.3,1x,F6.3,3(1x,f5.1),1x,f5.2)
999     stop
        end

        
