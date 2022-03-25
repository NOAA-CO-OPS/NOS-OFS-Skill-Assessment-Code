C   compiling with:  lf95 filter.f foufil.f utility.f -o filter.x

!#   DELT:       desired time interval (in hours) of observation, tide prediction and model outputs.
!#   DELT_O:     actual time interval (in hours) of observation, tide prediction
!#   DELT_M:     actual time interval (in hours) of model outputs.
!#   CUTOFF:     CUTOFF period (in hours) for Fourier filtering, =0 no filtering
!#   METHOD:     index of interpolation method
!#               0: cubic spline  1:Singular Value Decomposition(SVD); 
!#   CRITERIA1:  (in hours)means using linear and cubic spline interpolation
!#                when gap is less than criteria1    
!#   CRITERIA2:  (in hours) means using cubic spline or SVD interpolation method
!#                when criteria1 < gap < criteria2.
!#               fill gaps using missing value -999.0 while gap > criteria2
         parameter(num=190000)
         dimension time(num),wl(num),xnew(num),ynew(num),wlnew(num)
     1     ,wltmp(99),sp(num),diro(num),thour(num),u(num),v(num)
         character*80 filein,fileout,BUFFER
      real*8 jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,dayb,hourb
         CALL GETARG(1,BUFFER)
         READ(BUFFER,*)KINDAT
         CALL GETARG(2,BUFFER)
         READ(BUFFER,*)DELT
         CALL GETARG(3,BUFFER)
         READ(BUFFER,*)CUTOFF
         CALL GETARG(4,filein)
         CALL GETARG(5,fileout)
         DELT=DELT/60.0  !! 11/30/2006  convert from minutes into hours
         open(10,file=filein,form='formatted',status='old')
         open(11,file=fileout,form='formatted')
        i=0
        time0=0.0
        IF (KINDAT .EQ. 1)THEN
10       READ(10,*,end=20)ttmp,IYR0,ICM,ICD,IHR,IMN,(wltmp(n),n=1,4)
         if (ttmp .le. time0 .or. abs(wltmp(1)) .gt. 900.)goto 10
         if(i .eq. 0)then
           IYRS=IYR0
           yearb=IYRS
           monthb=1.
           dayb=1.
           hourb=0.
           jbase_date=JULIAN(yearb,monthb,dayb,hourb)
         endif
         yearb=IYR0
         monthb=ICM
         dayb=ICD
         hourb=IHR+IMN/60.0
         jday=JULIAN(yearb,monthb,dayb,hourb)
!         time1=jday-jbase_date+1
         time1=ttmp
         if(time1 .le. time0)goto 10
         i=i+1
         time(i)=time1
         time0=time1
         sp(i)=wltmp(1)
         diro(i)=wltmp(2)
         u(I)=wltmp(3)
         v(I)=wltmp(4)
!         u(I)=sp(i)*sin(diro(i)*3.1415926/180.)
!         v(I)=sp(i)*cos(diro(i)*3.1415926/180.)
         goto 10
20       continue
         close(10)
         NT=I
        ELSE IF (KINDAT .GE. 2)THEN
15       READ(10,*,end=25)ttmp,IYR0,ICM,ICD,IHR,IMN,wltmp(1)
!         if(ttmp .lt. 300.)goto 15
         if (abs(wltmp(1)) .gt. 90.)goto 15
         if(i .eq. 0)then
           IYRS=IYR0
           yearb=IYRS
           monthb=1.
           dayb=1.
           hourb=0.
           jbase_date=JULIAN(yearb,monthb,dayb,hourb)
         endif
         yearb=IYR0
         monthb=ICM
         dayb=ICD
         hourb=IHR+IMN/60.0
         jday=JULIAN(yearb,monthb,dayb,hourb)
         time1=jday-jbase_date+1
         if(time1 .le. time0)goto 15
         i=i+1
         time(i)=time1
         time0=time1
         wl(i)=wltmp(1)
         goto 15
25       continue
         close(10)
         NT=I
        ENDIF
        print *,'total number of observations = ',nt
        tstart=INT(time(1)*24.0/DELT)*DELT/24.
        tend=INT(time(nt)*24.0/DELT)*DELT/24.
        if(tstart.lt. time(1))tstart=tstart+DELT/24.
        if(tend.gt. time(nt))tend=tend-DELT/24.
        tstart=time(1)+cutoff/24./2.
        tend=time(NT)-cutoff/24./2.
        IF (KINDAT .EQ. 1)THEN
 !!          low-pass Fourior filtering
          call foufil(NT,DELT*60.,cutoff,U,Ynew)
          call foufil(NT,DELT*60.,cutoff,V,wlnew)
          DO i=1,NT
              IF( (ABS(YNEW(i)) .Gt. 900.) .or. 
     1                 (ABS(wlnew(i)) .GT. 900.) )then
                     sp(i)=-999.0
                     diro(i)=-999.0
              ELSE   
                 AAVN=wlnew(i)
                 AAVE=ynew(I)
                 CALL VELDIR (AAVN,AAVE,AANGLE,AVEL)
                 IF (AANGLE.GT.360.)   AANGLE = AANGLE - 360.
                 sp(i)=avel
                 diro(i)=AANGLE
              ENDIF
              u(i)=ynew(i)
              v(i)=wlnew(i)
          ENDDO
        ELSE IF (KINDAT .GE. 2)THEN
          call foufil(NT,DELT*60.,cutoff,WL,wlnew)
        ENDIF
        yearb=IYRS
        monthb=1.
        dayb=1.
        hourb=0.
        jbase_date=JULIAN(yearb,monthb,dayb,hourb)
        do i=1,NT
            IF(time(i) .GE.TSTART .and. time(i).LE.TEND)THEN
              jday=time(i)+jbase_date-1
              call GREGORIAN(jday,yearb,monthb,dayb,hourb)
              IYEAR=INT(yearb)
              ICM=int(monthb+0.001)
              ICD=INT(dayb+0.001)
              IHR=INT(hourb+0.001)
              IMN=INT((hourb-IHR)*60+0.1)
              IF (KINDAT .EQ. 1)THEN
                write(11,571)time(i),IYEAR,ICM,ICD,IHR,IMN,
     1          sp(i),diro(i),u(i),v(i)
              ELSE IF (KINDAT .GE. 2)THEN
                write(11,571)time(i),IYEAR,ICM,ICD,IHR,IMN,wlnew(i)
              ENDIF
            ENDIF
        enddo
  571   format(f10.5,I5,4I3,4f10.4)
        STOP 
        END
!         include 'foufil.f'
!         include './equal_interval.f'
!         include './svd.f'
!         include './utility.f'
