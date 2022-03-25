CC   lf95 generate_skillscore_table.f -o generate_skillscore_table.x
CC ifort -extend-source -O generate_skillscore_table.f -o generate_skillscore_table.x
      PARAMETER(NSTA=100,NMAX=30)
      character*120 OFS,fileshort,station_file,staid,stationname
      character*120 DIR_CTL,DIR_WORK,sid(NSTA),sname(NSTA)
      character*120 FOUT,CTMP,CTMP1,FMT,filein,criteria*5
      character*120 FVCOM, ROMS,SELFE,POM,FCTL, BUFFER
      REAL RMSE(NSTA,NMAX),SM(NSTA,NMAX),NOF(NSTA,NMAX),CF(NSTA,NMAX),
     1POF(NSTA,NMAX),SD(NSTA,NMAX),ASKILL(NSTA,NMAX)
      LOGICAL FEXIST
      NFDURATION=48
      NLOOP=INT(NFDURATION/6+0.01)+1
      write(*,*)'NLOOP=',NLOOP
      READ(5,"(a120)")OFS
      READ(5,"(a120)")DIR_CTL
      READ(5,"(a120)")DIR_WORK
      READ(5,*)KINDAT
      READ(5,*)NFDURATION     
      NLOOP=INT(NFDURATION/6+0.01)+2 !including nowcast 00, 06 12 ... 
      IF(KINDAT .EQ. 1)FCTL=TRIM(OFS)//'_cu_station.ctl'
      IF(KINDAT .EQ. 2)FCTL=TRIM(OFS)//'_wl_station.ctl'
      IF(KINDAT .EQ. 3)FCTL=TRIM(OFS)//'_temp_station.ctl'
      IF(KINDAT .EQ. 4)FCTL=TRIM(OFS)//'_salt_station.ctl'
      IF(KINDAT .EQ. 5)FCTL=TRIM(OFS)//'_cu_station.ctl'
      CLOSE(16)
      CLOSE(17)
      IF(KINDAT .EQ. 1)THEN
        FOUT=TRIM(OFS)//'_CU_skillscores.out'
        OPEN(16,file=TRIM(FOUT))
	BUFFER=TRIM(OFS)//' Skill Assessment Summary Table - Currents'
        WRITE(16,*)trim(BUFFER)
        WRITE(16,*)'BIAS (Model - Observation) in cm/s'
        WRITE(16,*)'RMSE - Root Mean Square Error in cm/s'
        WRITE(16,*)'SD - Standard Deviation in cm/s' 
        WRITE(16,*)'NOF - Nagative Outlier Frequency %'
        WRITE(16,*)'CF - Central Frequency %'
        WRITE(16,*)'NOF - Positive Outlier Frequency %'
        WRITE(16,*)'CORR - Correlation Coefficient between model and data'
        WRITE(16,*)'SKILL - non-dimensional skill' 
      ELSE IF(KINDAT .EQ. 2)THEN
         FOUT=TRIM(OFS)//'_WL_skillscores.out'
         OPEN(16,file=TRIM(FOUT) )
	 BUFFER=TRIM(OFS)
     1	 //' Skill Assessment Summary Table - Water Levels'
         WRITE(16,*)trim(BUFFER)
        WRITE(16,*)'BIAS (Model - Observation), in cm '
        WRITE(16,*)'RMSE - Root Mean Square Error in cm'         
        WRITE(16,*)'SD - Standard Deviation in cm' 
        WRITE(16,*)'NOF - Nagative Outlier Frequency %'
        WRITE(16,*)'CF - Central Frequency %'
        WRITE(16,*)'NOF - Positive Outlier Frequency %'
        WRITE(16,*)'CORR - Correlation Coefficient between model and data'
        WRITE(16,*)'SKILL - non-dimensional skill'
      ELSEIF(KINDAT .EQ. 3)THEN
         FOUT=TRIM(OFS)//'_TEMP_skillscores.out'
         OPEN(16,file=TRIM(FOUT))
	 BUFFER=TRIM(OFS)
     1	 //' Skill Assessment Summary Table - Temperature'
         WRITE(16,*)trim(BUFFER)
        WRITE(16,*)'BIAS (Model - Observation) Degree Celsius'
        WRITE(16,*)'RMSE - Root Mean Square Error in Degree Celsius'
        WRITE(16,*)'SD - Standard Deviation in Degree Celsius'
        WRITE(16,*)'NOF - Nagative Outlier Frequency %'
        WRITE(16,*)'CF - Central Frequency %'
        WRITE(16,*)'NOF - Positive Outlier Frequency %'
        WRITE(16,*)'CORR - Correlation Coefficient between model and data'
        WRITE(16,*)'SKILL - non-dimensional skill'
      ELSEIF(KINDAT .EQ. 4)THEN
         FOUT=TRIM(OFS)//'_SALT_skillscores.out'
         OPEN(16,file=TRIM(FOUT))
	 BUFFER=TRIM(OFS)
     1	 //' Skill Assessment Summary Table - Salinity'
         WRITE(16,*)trim(BUFFER)
        WRITE(16,*)'BIAS (Model - Observation) in PSU'
        WRITE(16,*)'RMSE - Root Mean Square Error in PSU'
        WRITE(16,*)'SD - Standard Deviation in PSU'
        WRITE(16,*)'NOF - Nagative Outlier Frequency %'
        WRITE(16,*)'CF - Central Frequency %'
        WRITE(16,*)'NOF - Positive Outlier Frequency %'
        WRITE(16,*)'CORR - Correlation Coefficient between model and data'
        WRITE(16,*)'SKILL - non-dimensional skill'
      ELSEIF(KINDAT .EQ. 5)THEN
         FOUT=TRIM(OFS)//'CU_Phase_skillscores.out'
         OPEN(16,file=TRIM(FOUT))
	 BUFFER=TRIM(OFS)
     1	 //' Skill Assessment Summary Table - Current Phase'
         WRITE(16,*)trim(BUFFER)
      ENDIF
!      BUFFER='STATION_ID  STA_NAME               NOWCAST  FH00   FH06'
!      BUFFER=trim(BUFFER)//'   FH12   FH18   FH24   FH30   FH36'
!      BUFFER=trim(BUFFER)//'   FH42   FH48'
!      WRITE(16,*)trim(BUFFER)
      open(1,file=TRIM(DIR_CTL)//'/'//trim(FCTL))
CC  generate table for each station, loop through all stations 
      NS=0     
10    read(1,*,err=999,end=999)staid,fileshort,stationname
      write(6,*)'fileshort=',fileshort
      write(6,*)'stationname=',trim(stationname)
      read(1,*)alat,alon,dirflood,sdepth
      ROMS=TRIM(DIR_WORK)//'/'//trim(fileshort)//'_table.out'
!      write(*,*)trim(ROMS)
      IF(KINDAT .EQ. 5)THEN
        ROMS=TRIM(DIR_WORK)//'/'//trim(fileshort)//'phase_table.out'
      ENDIF
      INQUIRE(FILE=trim(ROMS),EXIST=FEXIST)
      IF(.NOT.FEXIST)THEN
           write(6,*)'filein=',trim(fileshort), ' does not exist'
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
      write(*,100)'string=',trim(ROMS)
      IF(KINDAT .EQ. 4 )THEN
        read(ROMS,140)criteria,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR,SKILL1
      ELSEIF((KINDAT .EQ. 2) .OR. (KINDAT .EQ. 3) )THEN
        read(ROMS,130)criteria,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR,SKILL1
      ELSEIF((KINDAT .EQ. 1) .OR. (KINDAT .EQ. 5) )THEN
        read(ROMS,130)criteria,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR,SKILL1
      ENDIF
      WRITE(6,*)criteria,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR,SKILL1
130   format(11x,a5,14x,3f7.3,3f6.1,18x,2f6.2) 
140   format(14x,a5,11x,3f7.3,3f6.1,18x,2f6.2)
!S-s            3.5 24h  6961  -0.892  1.031  0.517   0.0 100.0   0.0    0.0  0.0        0.46
      NS=NS+1
      SM(NS,1)=SM1
      RMSE(NS,1)=RMSE1
      SD(NS,1)=SD1
      NOF(NS,1)=ANOF1
      CF(NS,1)=CF1
      POF(NS,1)=POF1
      CORR_C(NS,1)=CORR
      ASKILL(NS,1)=SKILL1
      IF(KINDAT .LE. 2)THEN
        read(10,100,err=10,end=10)CTMP
        ll0=index(ctmp,'FORECAST')
        do while (ll0 <= 0)
          read(10,100,err=10,end=10)CTMP
          ll0=index(ctmp,'FORECAST') 
        enddo        
        DO N=1,NLOOP  !read Extreme stats as well
         N1=N+1
         read(10,100,err=10)ROMS
S000-s000      3.5 24h   124   6.087  6.401  1.987   0.0   9.7  30.6 0.0198.0
T000-t000    3.0 c 24h   124  -0.461  0.548  0.299   0.0 100.0   0.0
0.0  0.0

!        read(ROMS,130)SM1,RMSE1,SD1,ANOF1,CF1,POF1
         read(ROMS,130)buffer,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR
         SM(NS,N1)=SM1
         RMSE(NS,N1)=RMSE1
         SD(NS,N1)=SD1
         NOF(NS,N1)=ANOF1
         CF(NS,N1)=CF1
         POF(NS,N1)=POF1
         CORR_C(NS,N1)=CORR
        ENDDO
        IF(KINDAT .EQ. 2)THEN
          DO N=1,NLOOP
            SM(NS,N)=SM(NS,N)*100.0
	    RMSE(NS,N)=RMSE(NS,N)*100.
	    SD(NS,N)=SD(NS,N)*100.0
          ENDDO
        ENDIF
      ELSEIF(KINDAT .GT. 2)THEN
        read(10,100,err=10,end=10)CTMP
        ll0=index(ctmp,'FORECAST')
        do while (ll0 <= 0)
          read(10,100,err=10,end=10)CTMP
          ll0=index(ctmp,'FORECAST') 
        enddo 
        DO N=1,NLOOP
          N1=N+1
          read(10,100,err=10,end=50)ROMS
          read(ROMS,130)buffer,SM1,RMSE1,SD1,ANOF1,CF1,POF1,CORR
          SM(NS,N1)=SM1
          RMSE(NS,N1)=RMSE1
          SD(NS,N1)=SD1
          NOF(NS,N1)=ANOF1
          CF(NS,N1)=CF1
          POF(NS,N1)=POF1
          CORR_C(NS,N1)=CORR
        ENDDO
      ENDIF	
50    CONTINUE
      stationname=trim(adjustL(stationname))	 
      LL=len(trim(stationname))
!      DO L=1,LL
        sname(NS)=trim(adjustL(stationname))
 !     ENDDO
      sid(NS)=trim(adjustL(staid))	
      if(LL .GT. 25)LL=25
      IF(KINDAT .EQ. 2)THEN
        write(FMT,'(a6,a1,I2.2,a1,I2.2,a1,a12)')'A8,2x,',
     1'A',LL,',',25-LL,'X',',10(1x,f6.2)'
        FMT="("//trim(FMT)//")"
!        print *,trim(FMT)
!        write(17,FMT)trim(staid),trim(stationname),(RMSE(NS,I),I=1,10)
!      STOP	
      ELSE        
        continue
      ENDIF
      write(*,*)'sname=',sname(NS)
      goto 10
999   CONTINUE
! Write output 
      WRITE(16,*)
      IF(KINDAT .EQ. 1)THEN
	BUFFER=TRIM(OFS)//' Nowcast Skill Assessment Summary Table - Currents'
        WRITE(16,*)trim(BUFFER)
      ELSE IF(KINDAT .EQ. 2)THEN
	BUFFER=TRIM(OFS)//' Nowcast Skill Assessment Summary Table - Water Levels'
        WRITE(16,*)trim(BUFFER)
      ELSEIF(KINDAT .EQ. 3)THEN
	BUFFER=TRIM(OFS)//' Nowcast Skill Assessment Summary Table - Temperature'
        WRITE(16,*)trim(BUFFER)
      ELSEIF(KINDAT .EQ. 4)THEN
	BUFFER=TRIM(OFS)//' Nowcast Skill Assessment Summary Table - Salinity'
        WRITE(16,*)trim(BUFFER)
      ELSEIF(KINDAT .EQ. 5)THEN
	BUFFER=TRIM(OFS)//' Nowcast Skill Assessment Summary Table - Current Phase'
        WRITE(16,*)trim(BUFFER)
      ENDIF
      WRITE(16,*)
      BUFFER='STATION_ID  STA_NAME                  BIAS   RMSE    SD    NOF    CF     POF   CORR SKILL'
      LL=len(trim(adjustL(BUFFER)) )
      write(FMT,'(a1,I3.3)')'A',LL
      FMT="("//trim(FMT)//")"
      WRITE(16,FMT)trim(BUFFER)
      DO N=1,NS      
        write(16,131)trim(sid(N)),sname(N),SM(N,1),RMSE(N,1),SD(N,1),
     1  NOF(N,1),CF(N,1),POF(N,1),CORR_C(N,1),ASKILL(N,1)
      ENDDO
      WRITE(16,*)
      WRITE(16,*)
      IF(KINDAT .EQ. 1)
     1WRITE(16,*)TRIM(OFS)//' - Root Mean Square Error (RMSE) in cm/s' 
      IF(KINDAT .EQ. 2)
     1WRITE(16,*)TRIM(OFS)//' - Root Mean Square Error (RMSE) in cm' 
      IF(KINDAT .EQ. 3)
     1WRITE(16,*)TRIM(OFS)//' - Root Mean Square Error (RMSE) in Deg.' 
      IF(KINDAT .EQ. 4)
     1WRITE(16,*)TRIM(OFS)//' - Root Mean Square Error (RMSE) in PSU' 
      BUFFER='STATION_ID  STA_NAME               NOWCAST  FH00   FH06'
      BUFFER=trim(BUFFER)//'   FH12   FH18   FH24   FH30   FH36'
      BUFFER=trim(BUFFER)//'   FH42   FH48'
      LL=len(trim(adjustL(BUFFER)) )
      write(FMT,'(a1,I3.3)')'A',LL
      FMT="("//trim(FMT)//")"
      WRITE(16,FMT)trim(adjustL(BUFFER))
      DO N=1,NS
        buffer=sname(N)
	LL=len(trim(buffer))
        if(LL .GT. 25)LL=25
        IF(KINDAT .EQ. 2)THEN
          write(FMT,'(a6,a1,I2.2,a1,I2.2,a1,a12)')'A8,2x,',
     1    'A',LL,',',25-LL,'X',',10(1x,f6.2)'
          FMT="("//trim(FMT)//")"
        ENDIF
        write(*,*)'FMT=',trim(FMT),'LL=',LL,'sname=',sname(N)	
        write(16,131)trim(sid(N)),sname(N),(RMSE(N,I),I=1,10)
      ENDDO
      WRITE(16,*)
      WRITE(16,*)
      IF(KINDAT .EQ. 1)
     1 WRITE(16,'(a90)')TRIM(OFS)//' - Percentage of abs(obs - model) <  in cm/s'
     1 //trim(criteria)//'  -- Central Frequency (CF)'
      IF(KINDAT .EQ. 2)
     1 WRITE(16,'(a90)')TRIM(OFS)//' - Percentage of abs(obs - model) <  in cm'
     1 //trim(criteria)//'  -- Central Frequency (CF)'
      IF(KINDAT .EQ. 3)
     1 WRITE(16,'(a90)')TRIM(OFS)//' - Percentage of abs(obs - model) <  in Deg.'
     1 //trim(criteria)//'  -- Central Frequency (CF)'
      IF(KINDAT .EQ. 4)
     1 WRITE(16,'(a90)')TRIM(OFS)//' - Percentage of abs(obs - model) <  in PSU'
     1 //trim(criteria)//'  -- Central Frequency (CF)'

      BUFFER='STATION_ID  STA_NAME               NOWCAST  FH00   FH06'
      BUFFER=trim(BUFFER)//'   FH12   FH18   FH24   FH30   FH36'
      BUFFER=trim(BUFFER)//'   FH42   FH48'
      LL=len(trim(adjustL(BUFFER)) )
      write(FMT,'(a1,I3.3)')'A',LL
      FMT="("//trim(FMT)//")"
      WRITE(16,FMT)trim(BUFFER)
      DO N=1,NS      
        write(16,131)trim(sid(N)),sname(N),(CF(N,I),I=1,10)
      ENDDO

      CLOSE(16)      
100   format(a120)
131   FORMAT(A8,2x,A25,10(1x,f6.2))      	
132   FORMAT(A8,2x,A25,10(1x,f6.2))      	
      stop
      end

