CC     KINDAT:  =1 for vector data (current speed and direction)
C               =2 dor scalar (water level)
C      DELT:    equally-spaced data time interval in minutes, = 6 minutes time interval
C      gonl:    longitude of the station in decimal degrees (positive for west longitude)
C      TM:      time meridian of the starting time 
C               =0 for Greenwich time, =75 for standard eastern time
C                
      PARAMETER(NT=100000)
      character*100 filestation,outputf0,FNAME,AQ,BUFFER*100 
      DIMENSION TIME(NT),VIN(NT),DIN(NT),SIYR(NT),SMON(NT)
     1          ,SDD(NT),SHH(NT),SMIN(NT),ifrst(NT),ilast(NT)
     2          ,XTMP(NT),YTMP(NT),U(NT),V(NT)
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(NT,2),RATIO 
      real*8 jday,jdayb,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      REAL*8 TIME8(NT),TSTART8,TFINISH8,DELT8,XTMP8(NT)
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
!      READ(BUFFER,*)XMAJOR  
      CALL GETARG(5,FNAME)
      CALL GETARG(6,BUFFER)
      READ(BUFFER,*)IYRS,MMS,IDDS,IHHS,IMINS 
      yearb=IYRS
      monthb=1.
      dayb=1.
      hourb=0.
      jbase_date=JULIAN(yearb,monthb,dayb,hourb)
      yearb=IYRS
      monthb=MMS
      dayb=IDDS
      hourb=IHHS+IMINS/60.0
      jdayb=JULIAN(yearb,monthb,dayb,hourb)
      stime=jdayb-jbase_date
      CALL GETARG(7,BUFFER)
      READ(BUFFER,*)IYRE,MME,IDDE,IHHE,IMINE
      yearb=IYRE
      monthb=MME
      dayb=IDDE
      hourb=IHHE+IMINE/60.0
      jdayb=JULIAN(yearb,monthb,dayb,hourb)
      etime=jdayb-jbase_date
      AQ='free format ascii file'
100     format(a7,1x,a4,2x,a40,2x,f8.3)
110     format(I2,1x,F5.2,1x,a1,1x,I3,1x,F5.2,1x,a1) 
!      call ncrght(FNAME,ncut)
      OPEN(2,file=trim(FNAME))
      N=0
      IF(KINDAT .EQ. 1)THEN
         I=1
5        READ(2,*,END=10)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I),DIN(I)       
         yearb=SIYR(I)
         monthb=SMON(I)
         dayb=SDD(I)
         hourb=SHH(I)+SMIN(I)/60.
         jdayb=JULIAN(yearb,monthb,dayb,hourb)
         TTMP=jdayb-jbase_date
         IF( (TTMP .GE. stime) .and.(TTMP .LE. etime))THEN
           TIME(I)=TTMP*24.
           I=I+1
         ENDIF
         GOTO 5
10      NMAX=I-1
!! fill up the small gaps in the original time series, i.e.  criteria1 < gap < criteria2 
         tstart=TIME(1)/24.
	 tfinish=TIME(NMAX)/24.
	 IMAXB=INT((tfinish-tstart)*24/delt)
	 method=1     ! using Singular Value Decomposition (SVD)
	 criteria1=2  ! in hours
	 criteria2=6  ! in hours
	 
         DO I=1,NMAX
            sp=VIN(I)
            dr=DIN(I)
            u(I)=sp*sin(dr*3.1415926/180.)
            v(I)=sp*cos(dr*3.1415926/180.)
         ENDDO 
! modified by AJ 06/16/2021
        TSTART8=TSTART
        TFINISH8=TFINISH
        DELT8=DELT
        DO I = 1,NMAX
          time8(I) = TIME(I)
        ENDDO
         CALL equal_interval(tstart8,tfinish8,delt8,delt8,method,
     1        criteria1,criteria2,TIME8,u,XTMP8,YTMP,NMAX,IMAXB)
         DO I=1,IMAXB
            U(I)=YTMP(I)
	 ENDDO   
         CALL equal_interval(tstart8,tfinish8,delt8,delt8,method,
     1        criteria1,criteria2,TIME8,v,XTMP8,YTMP,NMAX,IMAXB)
         DO I=1,IMAXB
            V(I)=YTMP(I)
	    TIME(I)=XTMP8(I)
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
      ELSE IF(KINDAT .EQ. 2)THEN
         I=1
15       READ(2,*,END=20)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I)
         yearb=SIYR(I)
         monthb=SMON(I)
         dayb=SDD(I)
         hourb=SHH(I)+SMIN(I)/60.
         jdayb=JULIAN(yearb,monthb,dayb,hourb)
         TTMP=jdayb-jbase_date
         IF( (TTMP .GE. stime) .and.(TTMP .LE. etime))THEN
           TIME(I)=TTMP*24.
           I=I+1
         ENDIF
         GOTO 15
20      NMAX=I-1
!! fill up the small gaps in the original time series, i.e.  criteria1 < gap < criteria2 
         tstart=TIME(1)/24.
	 tfinish=TIME(NMAX)/24.
	 IMAXB=INT((tfinish-tstart)*24/delt)
	 method=1     ! using Singular Value Decomposition (SVD)
	 criteria1=2  ! in hours
	 criteria2=6  ! in hours
        TSTART8=TSTART
        TFINISH8=TFINISH
        DELT8=DELT
        DO I = 1,NMAX
          time8(I) = TIME(I)
        ENDDO
        CALL equal_interval(tstart8,tfinish8,delt8,delt8,method,
     1        criteria1,criteria2,TIME8,VIN,XTMP8,YTMP,NMAX,IMAXB)
        DO I=1,IMAXB
            VIN(I)=YTMP(I)
	    TIME(I)=XTMP8(I)
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
      gap=1.5*delt
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
      ELSE IF(KINDAT .EQ. 2)THEN
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
!      AZI=XMAJOR
      N=NMAX
      NSPH=int(1./delt+0.0001)
      CVAR=0.0
      umean=0.0
      vmean=0.0
      INDATA=1
      VFAC=1.0
      ITYPE=3
C*****************************
      NJ=1
      IAND=0
      IREF=0
      IIT=0
      ISKIP9=0
      IEL=KINDAT-1
C****************************
        OPEN(20,file='ha_analysis.ctl')
        WRITE(20,'(a80)')FNAME
        WRITE(20,600)NJ,IAND,NSPH,IREF,IIT
        WRITE(20,610)AZI,TM,GONL,CVAR,IEL
        WRITE(20,'(a40)')AQ
        WRITE(20,620)IYRS,MONS,IDDS,IHHS,MINS,ISKIP9,N
        close(20)
        OPEN(25,file='ha_analysis.ctl')
        CALL HARM15D
        CLOSE(25)
!      ENDIF
500   FORMAT(f10.6,2I6)
510   FORMAT(f9.3,I8,I4,1x,3f9.3)
520   FORMAT(I5,4I4,f5.1,f9.3)
530   FORMAT(2I4,f9.3,2I4)
600   FORMAT(5I6)
610   FORMAT(f9.3,1x,3f9.3,I5)
620   FORMAT(I5,4I4,I4,I8)
      STOP
      END
C*****************************************************************************
      SUBROUTINE prcmp(ndata,x,y,pcd,rxy,RATIO)
C principle component computation
      parameter (nmax=120000)
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
                                                                     

CC***************************************************************************************
       SUBROUTINE HARM15D
C     PROGRAM CURAN (INPUT,OUTPUT,TAPE9,PUNCH)
C     ENVIRONMENTAL SCIENCE SERVICES ADMINISTRATION - C.+ G.S.
C     FOURIER - HARMONIC ANALYSIS    PROJECT 131403   15 DAYS
      COMMON/PARAM3/KINDAT,KAZI,AZI,AMP(37),pha(37),RDATA(100000,2)
     1      ,RATIO                                             
      COMMON/COSTX/CXX,OEX
      COMMON XYER,MONTH,DAY,STT
      CHARACTER*10 IDENS(16),IHD1,IHD2,IHD3
      CHARACTER*6 KAPR
      CHARACTER*8 AQ(4)
      CHARACTER*30 IDENST
       character*40 tapein
       character*40 altout
c      logical extst
      REAL*4 CR(20)
      DIMENSION ISP(25),SAVK(25),SAVH(25),SAVV(25)
      DIMENSION VIN(60000),DIN(60000),CHSIP(10),CRP(10),COVU(25),CF(25),
     1TXXIX1( 5,20), TXXIX2(5,20),CKAPAP(20), CKAP(20), CHSI(20),
     2BELL(20),BINGO(20),G(5),PAK(5),D(25),W(24),ZETA(5),CRC(5),PX(25)
      DIMENSION RES(4),ZEP(4),XTEMP(60000),VN(60000)
      DIMENSION CXX(30),DEG(40), TAB1(40), TAB2(40), TAB3(40), TAB4(40)
      DIMENSION TAB5(40), TAB6(40), OEX(5)
      DIMENSION TXXIXA(20),TXXIXB(20),TXXIXC(20),TXXIXD(20),TXXIXE(20)
      DIMENSION TXXIXF(20),TXXIXG(20),TXXIXH(20),TXXIXI(20),TXXIXJ(20)
      DATA SXM2,SXS2,SXN2,SXO1,SXK1/0.1610228,0.1666667,0.1579985,    
     1 0.0774613,0.0835615/                                                        
      DATA(TAB1(JZ),JZ=1,18)/6.5,13.9,17.9,19.0,17.6,14.7,10.8,6.4,1.5,   
     1 -3.5,-8.2,-12.5,-16.0,-18.4,-18.9,-16.8,-11.3,-2.8/                    
      DATA(TAB2(JZ),JZ=1,18)/-0.31,-0.25,-0.15,-0.04,0.07,0.17,0.25,    
     1 0.30,0.33,0.32,0.28,0.22,0.13,0.03,-0.08,-0.19,-0.28,-0.32/               
      DATA(TAB3(JZ),JZ=1,18)/3.2,7.2,10.8,13.7,2*15.4,13.2,8.6,1.9,-5.5,    
     1 -11.2,-14.6,-15.6,-14.7,-12.4,-9.1,-5.3,-1.1/                         
      DATA(TAB4(JZ),JZ=1,18)/0.26,0.23,0.18,0.10,0.01,-0.08,-0.17,-0.23,    
     1 -0.27,-0.25,-0.20,-0.13,-0.03,0.06,0.14,0.21,0.25,0.27/               
      DATA(TAB5(JZ),JZ=1,37)/-0.4,-1.0,-1.5,-2.0,-2.4,-2.7,-3.0,-3.2,    
     1 2*-3.4,-3.3,-3.1,-2.8,-2.4,-1.9,-1.4,-0.8,-0.2,0.5,1.1,1.6,2.1,    
     2 2.6,2.9,3.2,3.3,3.4,3.3,3.2,2.9,2.6,2.2,1.7,1.2,0.7,0.1,-0.4/               
      DATA(TAB6(JZ),JZ=1,37)/2*1.06,2*1.05,1.04,1.03,1.02,1.01,2*1.00,    
     1 0.99,0.98,0.97,0.96,2*0.95,3*0.94,2*0.95,2*0.96,0.97,0.98,0.99,    
     2 1.00,1.01,1.02,1.03,2*1.04,1.05,4*1.06/                                    
      DATA(TXXIXA( J ),J=1,20)/0.000,0.579,0.016,2*0.0,0.088,0.579,   
     1 0.080,0.054,0.026,0.672,0.016,0.672,7*0.0/                                
      DATA(TXXIXB( J ),J=1,20)/20*0.0/                                      
      DATA(TXXIXC( J ),J=1,20)/0.016,0.200,3*0.0,0.989,0.672,0.049,   
     1 2*0.997,0.579,0.016,0.214,7*0.0/                                         
      DATA(TXXIXD( J ),J=1,20)/4*0.0,0.024,8*0.0,0.207,0.615,0.024,  
     1 0.054,0.626,0.016,0.711/                                                 
      DATA(TXXIXE( J ),J=1,20)/3*0.0,0.024,9*0.0,2*0.626,0.024,0.990, 
     1 0.207,0.020,0.216/                                                     
      DATA(TXXIXF( J ),J=1,20)/0.0,262.,3.,2*0.0,18.,98.,344.,10.,175., 
     1 85.,357.,275.,7*0.0/                                                 
      DATA(TXXIXG( J ),J=1,20)/20*0.0/                                    
      DATA(TXXIXH( J ),J=1,20)/357.,259.,3*0.0,15.,275.,341.,7.,353.,  
     1 262.,354.,272.,7*0.0/                                                  
      DATA(TXXIXI( J ),J=1,20)/4*0.0,4.,8*0.0,96.,93.,9.,171.,269.,   
     1 357.,281./                                                                
      DATA(TXXIXJ( J ),J=1,20)/3*0.0,356.,9*0.0,91.,269.,4.,346.,264.,   
     1 353.,276./                                                            
      DATA(ISP(I),I = 1,24)/155854433,150410686,300821373,295284789,
     1 144966939,289841042,579682084,869523127,1159364169,284397295,
     2 278953548,139430356,161391017,149589314,133986609,128542862,
     3 300410667,300000000,600000000,900000000,299589333,294556253,
     4 285125831,134715145/
      DATA(PX(I),I=1,24)/1.,1.,2.,2.,1.,2.,4.,6.,8.,2.,2.,1.,1.,1.,1.,
     *                   1.,2.,2.,4.,6.,2.,2.,2.,1./
    5 FORMAT( 28X,7H1/F(K2), 6X,2HRA, 8X, 2HQA, 6X, 6HH~O(1))
    6 FORMAT(1H ,25X, F9.4, 2F10.5, F10.4/)
    7 FORMAT(21X,49HEquilibrium (V + U) and Elimination (F) Arguments /)
    9 FORMAT( 42X,22HAnalysis and Inference /)
   11 FORMAT(3X,33H Elimination of Component Effects //15X,4HR(A),5X
     1,4HZeta,5X, 5HKappa, 5X, 4HH(A) )
   14 FORMAT(1H1,20X,26HHARMONIC ANALYSIS, 15-DAYS, 2X,13HN(2) INFERRED)
   15 FORMAT(/ 14X,18HHarmonic Constants,1X,13H(H) and Kappa,A6/)
   16 FORMAT(16X,4HJ(1),6X,4HK(1),6X,4HK(2),6X,4HL(2),6X,4HM(1))
   17 FORMAT(16X,4HM(2),6X,4HM(4),6X,4HM(6),6X,4HM(8),6X,4HN(2))
   18 FORMAT(16X,4H(2N),6X,4HO(1),6X,4H(OO),6X,4HP(1),6X,4HQ(1))
   19 FORMAT(16X,4H(2Q),6X,4HR(2),6X,4HS(2),6X,4HS(4),6X,4HS(6))
  220 FORMAT(16X,4HT(2),4X,6HLAMBDA,5X,5HNU(2),4X,6HRHO(1))
   71 FORMAT(5X, 5HKappa, 5F10.2)
   72 FORMAT( 6X, 4HH(A), 5F10.4 /)
   73 FORMAT( 9X,11HZeta(Prime), 5F12.3)
   74 FORMAT(24F3.1)
   76 FORMAT( F5.0,7I5, F5.0, 2I5, F10.8, F5.0)
   77 FORMAT(60X,F3.0,24X,F4.2)
   78 FORMAT(24F3.2)
   79 FORMAT(24F3.0)
  221 FORMAT(I5)
  222 FORMAT(16I5)
  555 FORMAT(I5,F10.0,F10.2)
  802 FORMAT(30X,'Mean Current'/31X,'Speed    Dir'/31X,F7.2,F5.0)
  804 FORMAT (30X,'VALUES AFTER SERIES HAS BEEN DEMEANED')
  805 FORMAT(105X,'Mean Current'/106X,'Speed    Dir'/106X,F7.2,F5.0)
  808 FORMAT(8A10 / 8A10)
  810 FORMAT(/8A10/8A10/)
  811 FORMAT(27X,4HM(2),8X,4HN(2),8X,4HS(2),8X,4HO(1),8X,4HK(1))
  812 FORMAT(10X,8HR(Prime), 5F12.3/)
  813 FORMAT(27X,4HM(4),8X,4HM(6),8X,4HM(8),8X,4HS(4),8X,4HS(6))
  814 FORMAT( 12X,7H(V + U),1X, 5F12.3)
  815 FORMAT(13X, 4HF(A),2X, 5F12.5 /)
  816 FORMAT(27X,4HK(2),8X,4HL(2),8X,4H(2N),8X,4HR(2),8X,4HT(2))
  817 FORMAT(26X,6HLAMBDA,6X,5HMU(2),7X,5HNU(2),8X,4HJ(1),8X,4HM(1))
  818 FORMAT(26X,5HOO(1),8X,4HP(1),8X,4HQ(1),8X,4H(2Q),7X,6HRHO(1))
  819 FORMAT(15X,5HKappa, 5F12.3)
  820 FORMAT(15X, 4HZeta, 5F12.1)
  821 FORMAT(15X, 4HR(A), 5F12.4 /)
  822 FORMAT( 5X, 4HM(2), F10.5, 2F10.2, F10.5)
  823 FORMAT( 5X, 4HN(2), F10.5, 2F10.2, F10.5)
  824 FORMAT( 5X, 4HS(2), F10.5, 2F10.2, F10.5)
  825 FORMAT( 5X, 4HO(1), F10.5, 2F10.2, F10.5)
  826 FORMAT( 5X, 4HK(1), F10.5, 2F10.2, F10.5)
  223 FORMAT(1X/////1X,6HNSPH =,I5/////)
  827 FORMAT (1H1,////1X,33HHARMONIC ANALYSIS FOR TEMPERATURE,////)
  828 FORMAT (1H1,////1X,34HHARMONIC ANALYSIS FOR CONDUCTIVITY,////)
  829 FORMAT (1H1,////1X,30HHARMONIC ANALYSIS FOR PRESSURE,////)
  830 FORMAT (1H0,////1X,6HNSPH =,I10)
 2037 FORMAT(5F7.3,F5.0)
 2038 FORMAT(F5.0,I5,F5.0,F5.2,2F6.2,3I5,F5.0,2I5)
 2039 FORMAT(4A8)
 4040 FORMAT(1H0,10X, F6.0, 6F10.2)
 4041 FORMAT(1H0, 4F10.5)
 7358 FORMAT (/1H0,'JOB TERMINATED ..... NSPH = 0'/)
 7359 FORMAT(5X, F10.1,4X, F10.1,4X, 2F10.1)
 7360 FORMAT(1H0)
 7361 FORMAT(10X,5H PAGE,I3,7H  OF  3 )
 8109 FORMAT(4X,A6)
 8229 FORMAT(2F10.3,2F10.2,2I5,3X,A4,2A10)
      DO 8008 MAXJ = 1,20
      TXXIX1(1,MAXJ) = TXXIXA(MAXJ)
      TXXIX1(2,MAXJ) = TXXIXB(MAXJ)
      TXXIX1(3,MAXJ) = TXXIXC(MAXJ)
      TXXIX1(4,MAXJ) = TXXIXD(MAXJ)
      TXXIX1(5,MAXJ) = TXXIXE(MAXJ)
      TXXIX2(1,MAXJ) = TXXIXF(MAXJ)
      TXXIX2(2,MAXJ) = TXXIXG(MAXJ)
      TXXIX2(3,MAXJ) = TXXIXH(MAXJ)
      TXXIX2(4,MAXJ) = TXXIXI(MAXJ)
 8008 TXXIX2(5,MAXJ) = TXXIXJ(MAXJ)
       lout=6
       PRINT*, '       ********** 15 DAY HARMONIC ANALYSIS**********'
       print*, ' '
       read(25,4701) tapein
 4701 format(a40)
       call ncrght(tapein,ncut)
       tapein=tapein(1:ncut)
c      print*, '   analyze data in ',tapein
c      print*, ' '
c      print*, ' '
      DO 6005 I=1,16
      IDENS(I)='          '
 6005 CONTINUE
       IDENS(1) = 'Harmonic A'
       IDENS(2) = 'nalysis of'
       IDENS(3) = ' Data in  '
       IDENS(4) = TAPEIN(1:10)
       IDENS(5) = TAPEIN(11:20)
       IDENS(6) = TAPEIN(21:30)
       IDENS(7) = TAPEIN(31:40)
c      read(5,4711) altout
 4711  format(a40)
c      call ncrght(altout,ncut)
c      altout=altout(1:ncut)
c       inquire(file=altout,exist=extst)
c      if(extst) then
c       print*, '  sorry the output file specified in your control  '
c       print*,'   stream all ready exists'
c       print*,'   try again with new output file'
c       stop
c      end if
       altout = 'cons.out'
       open(unit=10,file=altout,form=
     1'formatted',access='sequential')
c
c
c      nj= number of harmonic analysis to do
c      iand is the default short input for reading from cdf tape
c      usually use iand =1 for the defaults
c
c       nsph is self explanatory
c
c          iref=0 use equil tide theory to infer the small const.)
c          iref=1 use a reference station results to infer
c          iit= number of iterations in eliminations and inference
c
c
      READ (25,*) NJ,IAND,NSPH,IREF,IIT
      print *, '  NJ IAND NSPH IREF  IIT'
      write (6,222) NJ,IAND,NSPH,IREF,IIT
C     *****WARNING*****
C     IREF AND IIT ARE NOT IMPLEMENTED!!!!!!!!!!
C
      IF(NSPH.EQ.0) STOP
      DO 1111 JN = 1,NJ
      JOBX=15
      N = 60000
      LCHOP = 24*NSPH
      N15 = 24*NSPH*15
      IF(JN.GT.1) GO TO 6003
      DO 6002 I=1,60000
      DIN(I) = 0.0
 6002 CONTINUE
      GO TO 6004
 6003 N = NSAVE
      K = KSAVE
      GO TO 8022
 6004 READ (25,*) AZI,TM,GONL,CVAR,IEL
      stm = tm
      print *, '  AZI   TM    GONL   CVAR  IEL'
      write (6,7134)AZI,TM,GONL,CVAR,IEL
      READ (25,2039)AQ
      PRINT 880,GONL
  880 FORMAT(/' GONL =',F10.5,' degrees west')
 7134 format(2x,4f6.2,i5)
 6006 FORMAT(2F5.0,3I5,F5.0,2I5,F5.0,I5,F5.0)
c       azi= major axis direction (if equals 0 do N-E analysis)
c
c       stm= time meridian for calculating kappa primes from kappa
c
c          cvar = magnetic compass variation
c         iel= xdata element to perform analysis on
c            =3 temperature
c            =0 do speed and direction
c            =4 do conductivity
c            =6 do pressure
c            =5 do residuals(major axis)
c            =8 do residuals (minor axis
c
c          tm = time meridian of data
c             =0.0 greenwich
c
c
      IF(IAND.LE.0) GO TO 6200
       open(11,file=tapein,status='old',form='unformatted',
     #access='sequential')
       print*, ' '
      IF (IEL.EQ.3)   PRINT 827
      IF (IEL.EQ.4)   PRINT 828
      IF (IEL.EQ.6)   PRINT 829
      GO TO 8022
 6200 CONTINUE
 !     open(9,file=tapein,status='old',form='formatted')
 8022 CONTINUE
c     READ (5,808) (IDENS(JK),JK = 1,16)
      READ (25,*) XYER,MONTH,DAY,STT,STTM,ISKIP9,N
      STT=STT+STTM/60.
      IDENS(9)  = '15-Day H.A'
      IDENS(10) = '.  Beginni'
      IDAY = DAY
      IYR = XYER
      write(IDENST,3011) month,iday,iyr,stt
 3011 FORMAT(3Hng ,I2,1H-,I2,1H-,I4,10H  at Hour ,F5.2,2H     )
      idens(11)(1:10)=idenst(1:10)
      idens(12)(1:10)=idenst(11:20)
      idens(13)(1:10)=idenst(21:30)
      NYNN = N15
      KAPR = ' '
      KP = 0
c     PRINT 14
      PRINT 810,(IDENS(KJ),KJ = 1,16)
       if(iskip9.ne.0) then
       print*,'Will skip ',iskip9,' data points as requested'
       endif
      PRINT 3045, XYER,MONTH,DAY,STT
      IF(JN.GT.1.AND.IEL.LE.0)  GO TO 3009
      IF(JN.GT.1.AND.IEL.GT.0)  GO TO 123
C
C     DATA READ IN (VARIOUS FORMATS)
C
      IF(IAND.EQ.1) GO TO 667
      IF(IEL) 668,102,101
  101 DO I=1,N
        VIN(I)=RDATA(I,1)
!      READ (9,AQ,END=1010)VIN(I)
      END DO
 1010 NE = I-1
      PRINT 3042,VIN(1),0.0,VIN(NE),0.0
  123 NECOM = 1
  208 PRINT 3043,ISKIP9+1,VIN(ISKIP9+1)
      PRINT 3043,N15+ISKIP9,VIN(N15+ISKIP9)
      IF(N15+ISKIP9.gt.NE) then
       print*, '    ATTENTION '
       print*, '    ATTEMPTED TO SKIP TOO MANY POINTS'
       print*, ' '
       print*, '   last point is # ',ne,'   needed  ',n15+iskip9
       go to 1111
       endif
      N = N15
      NSAVE = NE
      CALL TEMP(VIN,XTEMP,NE,ISKIP9)
C
C**** CALCULATE MEAN AND DEMEAN THE DATA
C
      CALL XMEAN(XTEMP,SUTZ,N15,AVE,2)
C
      GO TO 108
  102 DO I=1,N
        VIN(I)=RDATA(I,1)
        DIN(I)=RDATA(I,2)
!      READ(9,AQ,END = 112) VIN(I),DIN(I)
      END DO
  112 NE=I-1
      PRINT 3042,VIN(1),DIN(1),VIN(NE),DIN(NE)
      GO TO 3009
  667 CALL LOCATE (11,IHD1,IHD2,IHD3)
      CALL DBLOCK(VIN,DIN,N,IEL)
      NE = N
      IF (IEL.GT.0)   GO TO 123
      GO TO 3009
C
C   DATA ON CARDS (24F3.2/24F3.0), TILL 1 IN COL.80, CHECKS NOS.
C
  668 ISD1=0
!      ISS=1HS
!      IDD=1HD
      ION=0
      NB=1
      NE=24
 3000 READ (9,3001)(VIN(I),I=NB,NE),IS,IN1
      READ (9,3002) (DIN(I),I=NB,NE),ID,IN2,I80
      IF (ION.NE.0) GO TO 3005
      ION=97
      GO TO 3007
 3005 IF(IN22.LT.IN2) GO TO 3006
      PRINT 3040,IS,IN1,ID,IN2
      GO TO 3007
 3006 IF(IN2.EQ.(IN22+1)) GO TO 3007
      PRINT 3039,IN22
      STOP
 3007 IF(ISD1.NE.0) GO TO 3008
!      IF (IS.NE.ISS) PRINT 3037,IN1
!      IF (ID.NE.IDD) PRINT 3038,IN2
      ISD1=97
 3008 IN11=IN1
      IN22=IN2
      IF (I80.EQ.1) GO TO 3009
      NB=NE+1
      NE=NE+24
      GO TO 3000
 3009 PRINT 3043,ISKIP9+1,VIN(ISKIP9+1),DIN(ISKIP9+1)
      PRINT 3043,ISKIP9+N15,VIN(ISKIP9+N15),DIN(ISKIP9+N15)
       if(iskip9+n15.gt.ne) then
       print*, '    ATTENTION '
       print*, '    ATTEMPTED TO SKIP TOO MANY POINTS'
       print*, ' '
       print*, '  last point is ',ne,'   needed  ',iskip9+n15
       go to 1111
       endif
C
 3001 FORMAT(24F3.2,1X,A1,1X,I4)
 3002 FORMAT(24F3.0,1X,A1,1X,I4,I1)
 3040 FORMAT(' REPEAT CARD  NUMBER.   CARD NO. ',A1,I5/
     *24X,'CARD NO. ',A1,I5)
 3039 FORMAT (/' CARDS(S) MISSING BEFORE CARD NO.',I5)
 3037 FORMAT (' S NOT FOUND ON CARD NO.',I5)
 3038 FORMAT (' D NOT FOUND ON CARD NO.',I5)
 3041 FORMAT (/ I8,' DATA POINTS READ IN BUT WANTED TO USE',I8/)
 3042 FORMAT ( ' First data point read in',F7.2,F5.0/
     *' Last data point read in ',F7.2,F5.0)
 3043 FORMAT (/' Data point No.',I5,' is',F8.2,F5.0)
 3045 FORMAT(27X,'Begin at    Year  Mo Day Time'/39X,F5.0,I3,F4.0,F6.2/)
C
      NECOM = 4
      KSAVE = K
      NSAVE = NE
      IF(K.NE.2) K = 1
      N = N15
      PRINT 7370,AZI
 7370 FORMAT (/' Azimuth used =',F5.1)
 8355 CALL AZIM(DIN,VIN,AZI,VN,K,NE,CVAR)
      AMINOR = AZI + 90.0
      IF(AMINOR.GE.360.0) AMINOR = AMINOR - 360.0
  107 CALL TEMP(VN,XTEMP, NE,ISKIP9)
C
C**** CALCULATE AND DEMEAN THE DATA
C
      CALL XMEAN(XTEMP,SUTZ,N15,AVE,2)
      IF(NECOM.LT.3.OR.NECOM.GT.4) GO TO 83
      IF(K.EQ.1) GO TO 83
!      CALL VELDIR(ASAV,AVE,VL,DR,AZI)
      CALL VELDIR(ASAV,AVE,VL,DR)
      DR = DR + AZI
      IF(DR.GE.360.0) DR = DR - 360.0
      PRINT 802,VL,DR
      GO TO 108
   83 ASAV = AVE
C
  108 XM2 = SXM2/FLOAT(NSPH)
      XS2 = SXS2/FLOAT(NSPH)
      XN2 = SXN2/FLOAT(NSPH)
      XO = SXO1/FLOAT(NSPH)
      XK = SXK1/FLOAT(NSPH)
 8357 CALL FORAN(XTEMP,XM2,RES,ZEP,N,4)
      CHSIP(1) = ZEP(1)
      CHSIP(6) = ZEP(2)
      CHSIP(7) = ZEP(3)
      CHSIP(8) = ZEP(4)
      CRP(1) = RES(1)
      CRP(6) = RES(2)
      CRP(7) = RES(3)
      CRP(8) = RES(4)
      IF(JOBX -15) 109,109,110
  109 CHSIP(2) = 0.0
      CRP(2)  =  0.0
      GO TO 111
  110 CALL FORAN(XTEMP,XN2, RES,ZEP, N, 1)
      CHSIP(2)= ZEP(1)
      CRP(2) = RES(1)
  111 CALL FORAN(XTEMP,XS2, RES,ZEP, N, 3)
      CHSIP(3) = ZEP(1)
      CHSIP(9) = ZEP(2)
      CHSIP(10) =ZEP(3)
      CRP(3) = RES(1)
      CRP(9) = RES(2)
      CRP(10) = RES(3)
      N = N - LCHOP
      CALL FORAN( XTEMP,XO, RES, ZEP, N, 1)
      CHSIP(4) = ZEP(1)
      CRP(4) = RES(1)
      CALL  FORAN(XTEMP,XK, RES, ZEP, N, 1)
      CHSIP(5) = ZEP(1)
      CRP(5) =  RES(1)
      CALL TWOPI(CHSIP,10)
      N = N + LCHOP
      IF(NECOM.EQ.3.OR.NECOM.EQ.4) GO TO 8058
      GO TO 8057
 8058 IF(K.EQ.2) GO TO 8059
C     DETERMINATION OF VO+ U FOLLOWS WITH LOGF AND ARGUMENTS
 8057 CALL DAYXX(MONTH,DAY,STT,TM,DAYB,DAYM,GRBS,GRMS,CP)
      DO 7135 JAN = 1,19
      LXN = JAN + 18
      TAB1(LXN) = TAB1(JAN)
      TAB2(LXN) = TAB2(JAN)
      TAB3(LXN) = TAB3(JAN)
 7135 TAB4(LXN) = TAB4(JAN)
      CALL ASTRO (DAYB,DAYM,GRBS,GRMS,JOBX)
      PM = CXX(1)
      PL = CXX(2)
      SL = CXX(3)
      PS = CXX(4)
      PLM = CXX(5)
      SKYN = CXX(6)
      VI = CXX(9)
      V = CXX(10)
      XI = CXX(11)
      VP = CXX(12)
      VPP = CXX(13)
      P = CXX(14)
      AUL = CXX(15)
      AUM = CXX(16)
      PAPA = 0.0
      DO 3073 NAP = 1,37
      DEG(NAP) = PAPA
 3073 PAPA = PAPA + 10.0
C     DETERMINATION  OF EQUILIBRIUM  ARGUMENTS (V0 + U)
 7777 TML = CP - GONL
      CON = SL + TML
      COVU(1) = 2.0*(CON - PM + XI - V)
      COVU(2) = COVU(1) - (PM - PL)
      COVU(3) = 2.0*TML
      COVU(6) = 2.0*CON - VPP
      COVU(7) = 2.0*CON - PM - PL + AUL
      COVU(8) = COVU(2) - (PM - PL)
      COVU(9) = SL - PS + 180.0 + 2.0*(TML)
      COVU(10) = COVU(3) - (SL - PS)
      COVU(11) = COVU(1) - (2.0*(SL - PM)) - (PM - PL) + 180.0
      COVU(12) = COVU(1) + 2.0*(SL - PM)
      COVU(13) = COVU(12) + (PM - PL)
      COVU(21) = 2.0*(COVU(1))
      COVU(22) = 3.0*(COVU(1))
      COVU(23) = 4.0*(COVU(1))
      COVU(24) = 2.0*(COVU(3))
      COVU(25) = 3.0*(COVU(3))
C     DETERMINATION  OF NODE FACTOR RECIPROCALS (F)
      V = V*.0174533
      P = P*.0174533
      VI =VI*.0174533
      CRA = 1.0/SQRT(1.0 - 12.0*(SIN(.5*VI)**2/COS(.5*VI)**2)*COS(2.0*P)
     1+ 36.0*(SIN(.5*VI)**4/COS(.5*VI)**4))
      CF(1) = 1.0/(COS(.5*VI)**4/.91544)
      CF(2) = CF(1)
      CF(3) = 1.0
      CF(6) = 1.0/SQRT(19.0444*SIN(VI)**4 + 2.7702*SIN(VI)**2*COS(2.0*V)
     1+ 0.0981) - 0.000999
      CF(7) = CF(1)*CRA
      CF(8) = CF(1)
      CF(9) = 1.0
      CF(10) = 1.0
      CF(11) = CF(1)
      CF(12) = CF(1)
      CF(13) = CF(1)
      CF(21) = CF(1)**2
      CF(22) = CF(1)**3
      CF(23) = CF(1)**4
      CF(24) = 1.0
      CF(25) = 1.0
      CKKF = 1.0/CF(6)
      PLM = CXX(7)
      SKYN = CXX(8)
      V = CXX(10)
      P = CXX(17)
      AUM = CXX(19)
      VI = CXX(20)
      COVU(4) = CON - V - 2.0*(PM - XI) - 90.0
      COVU(5) = CON - VP + 90.0
      COVU(14) = CON + PM - PL - V + 90.0
      COVU(15) = CON - PM + AUM
      COVU(16) = CON - V + 2.0*(PM - XI) + 90.0
      COVU(17) = TML + 270.0 - SL
      COVU(18) = COVU(4) - (PM - PL)
      COVU(19) = COVU(18) - (PM - PL)
      COVU(20) = COVU(18) + 2.0*(SL - PL)
      CALL TWOPI(COVU,25)
      SAVV(01) = COVU(14)
      SAVV(02) = COVU(05)
      SAVV(03) = COVU(06)
      SAVV(04) = COVU(07)
      SAVV(05) = COVU(15)
      SAVV(06) = COVU(01)
      SAVV(07) = COVU(21)
      SAVV(08) = COVU(22)
      SAVV(09) = COVU(23)
      SAVV(10) = COVU(02)
      SAVV(11) = COVU(08)
      SAVV(12) = COVU(04)
      SAVV(13) = COVU(16)
      SAVV(14) = COVU(17)
      SAVV(15) = COVU(18)
      SAVV(16) = COVU(19)
      SAVV(17) = COVU(09)
      SAVV(18) = COVU(03)
      SAVV(19) = COVU(24)
      SAVV(20) = COVU(25)
      SAVV(21) = COVU(10)
      SAVV(22) = COVU(11)
      SAVV(23) = COVU(13)
      SAVV(24) = COVU(20)
      V = V*0.0174533
      P = P*0.0174533
      VI =VI*0.0174533
      CQA = 1.0/SQRT(2.310 + 1.435*COS(2.0*P))
      CF(4) = 1.0/(SIN(VI)*COS(.5*VI)**2/.37988)
      CF(5) = 1.0/SQRT(0.8965*SIN(2.0*VI)**2 + 0.6001*SIN(2.0*VI)*COS(V)
     1+ 0.1006)
      CF(14) = 1.0/(SIN(2.0*VI)/0.72137)
      CF(15) = CF(4)*CQA
      CF(16) = 1.0/(SIN(VI)*SIN(.5*VI)**2/0.016358)
      CF(17) = 1.0
      CF(18) = CF(4)
      CF(19) = CF(18)
      CF(20) = CF(19)
      P = CXX(17)
      VI = CXX(20)
      V = CXX(10)
      C5F = CF(5)
C     ENTER  SECONDARY EFFECTS AND  DETERMINE KAPPA
C     (WITH INFERENCE)
      CALL SASTR(DEG,C5F,CKKF,SMAC,PROD,ACCP,RESAM,37,TAB1,TAB2,TAB3,TAB
     14,TAB5,TAB6,SL,VPP,PS,VP,CF,25)
 8059 DO 50 I = 1,5
   50 CKAPAP(I) = CHSIP(I) + COVU(I)
      CALL TWOPI(CKAPAP, 5)
      AHNAT = CF(4)*CRP(4)
      CKAPAP(2) = 0.0
      CKAP(1) = CKAPAP(1)
      CKAP(3) = CKAPAP(3) + SMAC
      SMD = CKAP(3) - CKAP(1)
      CALL ONEPI( SMD)
      CKAP(2) = CKAP(1) - 0.536*(SMD)
      CKAP(4) = CKAPAP(4)
      CKAP(5) = CKAPAP(5) + ACCP
      CKAP(6) = CKAP(3)
      DMN = CKAP(1) - CKAP(2)
      CALL ONEPI(DMN)
      CKAP(7) = CKAP(1) + DMN
      CKAP(8) = CKAP(2) - DMN
      CKAP(9) = CKAP(3)
      CKAP(10) = CKAP(3)
      CKAP(11) = CKAP(1) + 0.464*(SMD)
      CKAP(12) = CKAP(1) - SMD
      CKAP(13) = CKAP(1) - 0.866*(DMN)
      OKD = CKAP(5) - CKAP(4)
      CALL ONEPI(OKD)
      CKAP(14) = CKAP(5) + 0.5*(OKD)
      CKAP(15) = CKAP(5) - 0.5*(OKD)
      CKAP(16) = CKAP(5) + OKD
      CKAP(17) = CKAP(5)
      CKAP(18) = CKAP(5) - 1.5*(OKD)
      CKAP(19) = CKAP(5) - 2.0*(OKD)
      CKAP(20) = CKAP(5) - 1.43*(OKD)
      CALL TWOPI(CKAP,20)
C     DETERMINATION OF ZETA(A) AND R(A) FOLLOWS
C     (WITH INFERENCE)
      DO 51 I = 1,20
   51 CHSI(I) = CKAP(I) - COVU(I)
      CALL TWOPI(CHSI, 20)
      CR(1) = CRP(1)
      CR(2) = 0.194*CR(1)
      CR(3) = CRP(3)/PROD
      CR(4) = CRP(4)
      CR(5) = CRP(5)/RESAM
      CR(6) = 0.272*CR(3)/CF(6)
      CR(7) = 0.145*CR(2)/CRA
      CR(8) = 0.133*CR(2)
      CR(9) = 0.008*CR(3)
      CR(10) = 0.059*CR(3)
      CR(11) = 0.007*CR(1)
      CR(12) = 0.024*CR(1)
      CR(13) = 0.194*CR(2)
      CR(14) = 0.079*AHNAT/CF(14)
      CR(15) = 0.071*CR(4)/CQA
      CR(16) = 0.043*AHNAT/CF(16)
      CR(17) = 0.331*CR(5)*CF(5)
      CR(18) = 0.194*CR(4)
      CR(19) = 0.026*CR(4)
      CR(20) = 0.038*CR(4)
C     ELIMINATION OF COMPONENT EFFECTS AND DETERMINE KAPPA AND H(A)
      DO 67  I = 1,5
      SUM = 0.0
      SOQ = 0.0
      IF(I - 2) 53,52,53
   52 PAK(I)= 0.0
      CRC(I)= 0.0
      ZETA(I) = 0.0
      G(I) = 0.0
      GO TO  67
   53 DO 61 J = 1, 20
  441 BELL(J) = CR(J)*TXXIX1(I,J)
      IF( BELL(J) - 0.0005)54,55,55
   54 BELL(J) = 0.0
   55 BINGO(J)= TXXIX2(I,J) - CHSI(J) + CHSIP(I)
      IF( BINGO(J)) 56,59,57
   56 BINGO(J) = (BINGO(J) + 360.0)*0.0174533
      GO TO 60
   57 IF( BINGO(J) - 360.0) 59,59,58
   58 BINGO(J) = ( BINGO(J) - 360.0)*0.0174533
      GO TO 60
   59 BINGO(J) = BINGO(J)*0.0174533
   60 SUM = SUM + BELL(J)*SIN(BINGO(J))
   61 SOQ = SOQ + BELL(J)*COS(BINGO(J))
      QOS = CRP(I) - SOQ
      CALL FITAN(SUM,QOS,DCCH,2)
      CDCH = DCCH*0.0174533
      CRC(I) = 0.0
      IF(DCCH.EQ.90.0.OR.DCCH.EQ.270.0) GO TO 66
      CRC(I) = (CRP(I) - SOQ)/ COS(CDCH)
   66 G(I) = CRC(I)*CF(I)
      ZETA(I) = CHSIP(I) + CDCH*57.29578
      PAK(I)  = ZETA(I) + COVU(I)
   67 CONTINUE
      CALL TWOPI(ZETA,5)
      CALL TWOPI(PAK, 5)
      D(1) = 0.079*G(4)
      D(2) = G(5)
      D(3) = 0.272*G(3)
      D(4) = 0.028*G(1)
      D(5) = 0.071*G(4)
      D(6) = G(1)
C     INFERENCE  OF     H(A)    FROM MAJOR CONSTITUENTS
      D(7) = CRP(6)*CF(21)
      D(8) = CRP(7)*CF(22)
      D(9) = CRP(8)*CF(23)
      D(10) = 0.194*D(6)
      D(11) = 0.133*D(10)
      D(12) = G(4)
      D(13) = 0.043*G(4)
      D(14) = 0.331*G(5)
      D(15) = 0.194*G(4)
      D(16) = 0.026*G(4)
      D(17) = 0.008*G(3)
      D(18) = G(3)
      D(19) = CRP(9)*CF(24)
      D(20) = CRP(10)*CF(25)
      D(21) = 0.059*G(3)
      D(22) = 0.007*G(1)
      D(23) = 0.194*D(10)
      D(24) = 0.038*G(4)
C     INFERENCE  OF   KAPPA(A)  FROM MAJOR CONSTITUENTS
      PAFX = PAK(5) - PAK(4)
      GAFX = PAK(5) - PAK(4)
      FPAX = PAK(5) + PAK(4)
      PATTO = PAK(3) - PAK(1)
      CALL ONEPI(PAFX)
      IF(ABS(GAFX).LT.180.0) GO TO 6781
      IF(GAFX)6779,6781,6780
 6779 FPAX = FPAX + 360.0
      GO TO 6781
 6780 FPAX = FPAX - 360.0
 6781 CALL ONEPI(PATTO)
      PAMNN = PAK(1) - 0.536*PATTO
      CALL TWOPI(PAMNN,1)
      PAMN = PAK(1) - PAMNN
      CALL ONEPI(PAMN)
      W(1)  = PAK(5) + 0.5*PAFX
      W(2) = PAK(5)
      W(3) = PAK(3)
      W(4) = 2.0*PAK(1) - PAMNN
      W(5) = 0.5*FPAX
      W(6) = PAK(1)
      W(7) = CHSIP(6) + COVU(21)
      W(8) = CHSIP(7) + COVU(22)
      W(9) = CHSIP(8) + COVU(23)
      W(10) = PAMNN
      W(11) = 2.0*W(10) - PAK(1)
      W(12) = PAK(4)
      W(13) = 2.0*PAK(5) - PAK(4)
      W(14) = PAK(5)
      W(15) = PAK(5) -1.5*PAFX
      W(16) = PAK(5) - 2.0*PAFX
      W(17) = PAK(3)
      W(18) = PAK(3)
      W(19) = CHSIP(9) + COVU(24)
      W(20) = CHSIP(10) + COVU(25)
      W(21) = PAK(3)
      W(22) = PAK(3) - 0.536*PATTO
      W(23) = PAK(1) - 0.866*PAMN
      W(24) = PAK(5) - 1.43*PAFX
      CALL TWOPI(W,24)
C
C
      IPAGE = 1
      if(iel.le.0) then
      if(k.eq.1) write(idens(14),701) nint(azi)
      if(k.eq.2) write(idens(14),701) nint(aminor)
  701 format(6Halong ,I3,X)
      idens(15) = 'degrees   '
      end if
c     PRINT 14
      PRINT 810,(IDENS(KJ),KJ = 1,16)
c     PRINT 7361, IPAGE
      PRINT 811
      PRINT 73, (CHSIP(I),I = 1, 5)
      PRINT 812,(CRP(I),I = 1, 5)
      PRINT 813
      PRINT 73, (CHSIP(I),I = 6,10)
      PRINT 812,(CRP(I),I = 6,10)
      PRINT 5
C     PRINT EQUILIBRIUM ARGUMENTS
      PRINT 6, CKKF,CRA,CQA,AHNAT
      PRINT 7
      PRINT 811
      PRINT 814,(COVU(I),I = 1, 5)
      PRINT 815,(CF(I),I = 1, 5)
      PRINT 816
      PRINT 814,(COVU(I),I = 6,10)
      PRINT 815,(CF(I),I = 6,10)
      PRINT 817
      PRINT 814,(COVU(I),I =11,15)
      PRINT 815,(CF(I),I =11,15)
      PRINT 818
      PRINT 814,(COVU(I),I =16,20)
      PRINT 815,(CF(I),I =16,20)
      PRINT 813
      PRINT 814,(COVU(I),I =21,25)
      PRINT 815,(CF(I),I =21,25)
      IPAGE = 2
c     PRINT 14
c     PRINT 810,(IDENS(KJ),KJ = 1,16)
c     PRINT 7361, IPAGE
C     PRINT ANALYSIS AND INFERENCE
      PRINT 9
      PRINT 811
      PRINT 819,(CKAP(I),I = 1, 5)
      PRINT 820,(CHSI(I),I = 1, 5)
      PRINT 821,(CR(I),I = 1, 5)
      PRINT 816
      PRINT 819,(CKAP(I),I = 6,10)
      PRINT 820,(CHSI(I),I = 6,10)
      PRINT 821,(CR(I),I = 6,10)
      PRINT 817
      PRINT 819,(CKAP(I),I =11,15)
      PRINT 820,(CHSI(I),I =11,15)
      PRINT 821,(CR(I),I =11,15)
      PRINT 818
      PRINT 819,(CKAP(I),I =16,20)
      PRINT 820,(CHSI(I),I =16,20)
      PRINT 821,(CR(I),I =16,20)
C     PRINT ELIMINATION
      PRINT 11
      PRINT 822, CRC(1),ZETA(1),PAK(1),G(1)
      PRINT 823, CRC(2),ZETA(2),PAK(2),G(2)
      PRINT 824, CRC(3),ZETA(3),PAK(3),G(3)
      PRINT 825, CRC(4),ZETA(4),PAK(4),G(4)
      PRINT 826, CRC(5),ZETA(5),PAK(5),G(5)
      IPAGE = 3
c     PRINT 14
c     PRINT 810,(IDENS(KJ),KJ = 1,16)
c     PRINT 7361, IPAGE
 8106 PRINT 15,KAPR
C     PRINT H(A) AND KAPPA
      PRINT 16
      PRINT 71,(W(I),I = 1,5)
      PRINT 8109,KAPR
      PRINT 72,(D(I),I = 1,5)
      PRINT 17
      PRINT 71,(W(I),I = 6,10)
      PRINT 8109,KAPR
      PRINT 72,(D(I),I = 6,10)
      PRINT 18
      PRINT 71,(W(I),I = 11,15)
      PRINT 8109,KAPR
      PRINT 72,(D(I),I = 11,15)
      PRINT 19
      PRINT 71,(W(I),I = 16,20)
      PRINT 8109,KAPR
      PRINT 72,(D(I),I = 16,20)
      PRINT 220
      PRINT 71,(W(I),I = 21,24)
      PRINT 8109,KAPR
      PRINT 72,(D(I),I = 21,24)
      IF(KP.EQ.1) GO TO 8120
      IF(IEL.LE.0) GO TO 8107
c     PRINT 14
      PRINT 810,(IDENS(KJ),KJ = 1,16)
      EQUIB = PM - XI - 90.0
      CALL TWOPI(EQUIB,1)
      CALL ONE80(W,D,GONL,EQUIB,DUM)
 8107 continue
c     PRINT 14
c     PRINT 810,(IDENS(KJ),KJ=1,16)
      IF(K.EQ.1) PRINT 806,AZI
      IF(K.EQ.2) PRINT 807,AMINOR
      PRINT 803,STM
  803 FORMAT(30X,'Kappa Primes are for Time Meridian',F5.0)
  806 FORMAT(/' Major Axis (',F4.0,')')
  807 FORMAT(/' Minor Axis (',F4.0,')')
C
C**** CALCULATE THE MEAN EVENTHOUGH THE SERIES HAS BEEN DEMEANED
C
c     CALL XMEAN(XTEMP,SUTZ,N15,AVEDM,1)
C
      IF(NECOM.LT.3.OR.NECOM.GT.4) GO TO 809
      IF(K.EQ.1) GO TO 809
!      CALL VELDIR (ASAVDM,AVEDM,VLDM,DRDM,AZI)
      CALL VELDIR (ASAVDM,AVEDM,VLDM,DRDM)
      DRDM = DRDM + AZI
      IF(DRDM.GE.360.0) DRDM = DRDM - 360.0
c     PRINT 804
c     PRINT 802,VLDM,DRDM
      GO TO 8809
  809 ASAVDM = AVEDM
c     PRINT 7360
c     PRINT 7360
C
C     CHANGING KAPPA TO KAPPA PRIME
C
 8809 DO 8110 I=1,24
      SPD = ISP(I)/10000000.
      W(I) = W(I) + PX(I)*GONL - SPD*STM/15.0
      CALL TWOPI(W(I),1)
 8110 CONTINUE
      CALL CARDS(W,D,AVE,NSPH,TFAC,IDENS,4,AZI)
      KP = 1
      KAPR = ' Prime'
C
      IF(NECOM.EQ.4.AND.K.EQ.2) GO TO 8106
      IF(NECOM.EQ.3.AND.K.EQ.2) GO TO 8106
      DO 8118 IQ=1,24
      SAVK(IQ) = W(IQ)
      SAVH(IQ) = D(IQ)
 8118 CONTINUE
      GO TO 8106
 8120 KP = 0
      N = NYNN
      KAPR = '      '
      IF(NECOM.EQ.3.AND.K.EQ.1) GO TO 5998
      IF(NECOM.EQ.4.AND.K.EQ.1) GO TO 5998
      IF(NECOM.NE.4.AND.NECOM.NE.3) GO TO 5899
C
C    DETERMINE ELLIPSE PARAMETERS
C
      idens(14) = '          '
      idens(15) = '          '
c     WRITE (LOUT,14)
      WRITE (LOUT,810)(IDENS(KJ),KJ=1,16)
      CALL ELIPSE (D,SAVH,W,SAVK,AZI,STM,GONL)
      PRINT 805, VL,DR
      GO TO 5899
C
 5998 K = 2
      if(necom.eq.3.or.necom.eq.4) write(idens(14),701) nint(aminor)
c     PRINT 14
c     PRINT 810,(IDENS(KJ),KJ = 1,16)
      GO TO (5899,5899,8355,8355),NECOM
 5899 CONTINUE
 1111 CONTINUE
  !    STOP
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C     FOURIER ANALYSIS -  DETERMINE ZETA(PRIME)AND R(PRIME)
      SUBROUTINE FORAN(SER, SP, RES, ZEP, N, M )
      DIMENSION  SER(N), RES(M), ZEP(M)
  200 FORMAT(1H0,10X, I5, 5X, 2F10.4, 5X, 2F12.3)
      DO 55 JR = 1, M
      RES(JR) = 0.0
      ZEP(JR) = 0.0
      C = 0.0
      S = 0.0
      DO 100 L = 1,N
      C = C + COS(FLOAT(L-1)*(3.141593*(FLOAT(JR) *SP)))*SER(L)
  100 S = S + SIN(FLOAT(L-1)*(3.141593*(FLOAT(JR) *SP)))*SER(L)
      C = 2.0/ FLOAT(N)* C
      S = 2.0/ FLOAT(N)* S
      IF(C.EQ.0.0) GO TO 14
      X = S/C
      GO TO 15
   14 IF(S)16,55,17
   16 ZEP(JR)= 270.0
      GO TO 18
   17 ZEP(JR) = 90.0
      GO TO 18
   15 ZEP(JR) = ATAN(X)*57.29578
      IF(S) 32,32,34
   32 IF(C) 33,55,37
   33 ZEP(JR) = 180.0 + ZEP(JR)
      GO TO 18
   34 IF(C) 35,17,36
   35 ZEP(JR) = ZEP(JR) + 180.0
      GO TO 18
   37 ZEP(JR) = ZEP(JR) + 360.0
      GO TO 18
   36 ZEP(JR) = ZEP(JR)*1.0
   18 RES(JR) = SQRT( C**2 + S**2)
c     PRINT 200, JR, C, S, ZEP(JR), RES(JR)
   55 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C     TEST FOR 180 DEGREE DIFFERENCE BETWEEN TWO ANGLES
      SUBROUTINE ONEPI( AUX)
      IF(ABS(AUX) -180.0)170,170,163
  163 IF(AUX)164,164,165
  164 AUX = AUX + 360.0
      GO TO 170
  165 AUX = AUX - 360.0
  170 RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE AZIM(DX,VR,A,VN,K,J,COMPV)
      DIMENSION DX(J), VR(J), VN(J)
      DO 100 I = 1,J
  100 VN(I) = (DX(I) + COMPV - A )*0.0174533
      GO TO (1,2),K
    1 DO 66 I = 1,J
   66 VN(I) = VR(I)*COS(VN(I))
      GO TO 109
    2 DO 77 I = 1,J
   77 VN(I) = VR(I)*SIN(VN(I))
  109 RETURN
      END
C
C ----------------------------------------------------------------------
C
C     STORE SERIES FOR FOURIER ANALYSIS
      SUBROUTINE TEMP(XSER,XTEMP,I,ISKIP9)
      DIMENSION   XSER(I) , XTEMP(I)
      DO 100 J = 1, I
      IF(J.LE.I-ISKIP9) THEN
      XTEMP(J) = XSER(J+ISKIP9)
      ELSE
      XTEMP(J) = 0.
      END IF
  100 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C     INTERPOLATE TABLES 21,22,23,24,25,AND 26
      SUBROUTINE TERPO(ANG,TXX,LID,ANSWER,AUGM)
      DIMENSION ANG(LID), TXX(LID)
    1 FORMAT(///21H JOB TERMINATED......,26HFAILED TO INTERPOLATE ANG.)
      LIN = LID - 1
      DO 200  ION = 1,LIN
      IF(AUGM.GT.ANG(ION).AND.AUGM.LT.ANG(ION+1)) GO TO 28
      IF(AUGM.EQ.ANG(ION)) GO TO 26
      IF(AUGM.EQ.ANG(ION+1)) GO TO 27
      GO TO 200
   26 ANSWER = TXX(ION)
      GO TO 201
   27 ANSWER = TXX(ION+1)
      GO TO 201
   28 REW = AUGM*0.1
      WER = REW - AINT(REW)
      FRA = WER*10.0
      RET = ABS(TXX(ION))
      TEX = ABS(TXX(ION+1))
      IF(RET.GT.TEX) GO TO 29
      IF(RET.LT.TEX) GO TO 30
      IF(RET.EQ.TEX) GO TO 31
   29 DIFR = (RET - TEX)*0.1
      IF(TXX(ION).GE.0.0.AND.TXX(ION+1).LT.0.0) DIFR = (RET + TEX)*0.1
      IF(TXX(ION).LT.0.0.AND.TXX(ION+1).GE.0.0) DIFR = (RET + TEX)*0.1
      DIFX = DIFR*FRA
      CHS = SIGN(DIFX,TXX(ION))
      ANSWER = TXX(ION) - CHS
      GO TO 201
   30 DIFR = (TEX - RET)*0.1
      IF(TXX(ION).GE.0.0.AND.TXX(ION+1).LT.0.0) DIFR = (RET + TEX)*0.1
      IF(TXX(ION).LT.0.0.AND.TXX(ION+1).GE.0.0) DIFR = (RET + TEX)*0.1
      DIFX = DIFR*FRA
      CHS = SIGN(DIFX,TXX(ION+1))
      ANSWER = TXX(ION) + CHS
      GO TO 201
   31 IF(TXX(ION).GE.0.0.AND.TXX(ION+1).GE.0.0) GO TO 26
      IF(TXX(ION).LT.0.0.AND.TXX(ION+1).LT.0.0) GO TO 26
      IF(TXX(ION).GE.0.0.AND.TXX(ION+1).LT.0.0) GO TO 32
      IF(TXX(ION).LT.0.0.AND.TXX(ION+1).GE.0.0) GO TO 32
      PRINT 1
      STOP
   32 DIFR = (TXX(ION+1) - TXX(ION))*0.1
      DIFX = DIFR*FRA
      ANSWER = TXX(ION) + DIFX
      GO TO 201
  200 CONTINUE
  201 RETURN
      END
C
C ----------------------------------------------------------------------
C
C     DETERMINE SECONDARY CONSTITUENT EFFECTS
      SUBROUTINE SASTR(DEG,C5F,CKKF,ACS,APS,ACK,APK,IOT,TAB1,TAB2, TAB3,
     1TAB4,TAB5,TAB6,SL,VPP,PS,VP,CF,JOP)
      DIMENSION DEG(IOT),TAB1(IOT),TAB2(IOT),TAB3(IOT),TAB4(IOT),TAB5(IO
     1T),TAB6(IOT),CF(JOP)
 2043 FORMAT(3x,4F10.5)
 2044 FORMAT(/6X,4HSMAC,6X,4HPROD,6X,4HACCP,4X,5HRESAM)
 2045 FORMAT( 2X,6F10.4/)
      AMT = SL - 0.5*VPP
      CALL TWOPI(AMT, 1)
      CALL TERPO(DEG,TAB3,37,RAK,AMT)
      TAZ = RAK
      AMR = SL - PS
      CALL TWOPI(AMR,1)
      CALL TERPO(DEG,TAB5,37,RAL,AMR)
      TAZ1 = RAL
      ACS = RAK*CKKF + RAL
      CALL TERPO(DEG,TAB4,37,SKO,AMT )
      TAZ2 = SKO
      CALL TERPO(DEG,TAB6,37,STOX,AMR)
      TAZ3 = STOX
      APS = STOX*(1.0 + (SKO*CKKF))
      EALT = SL - 0.5*VP
      CALL TWOPI(EALT,1)
      CALL TERPO(DEG,TAB1,37,PKX,EALT)
      TAZ4 = PKX
      ACK = PKX*C5F
      CALL TERPO(DEG,TAB2,37,PKXA,EALT)
      TAZ5 = PKXA
      APK = 1.0 + (C5F*PKXA)
      PRINT 2044
      PRINT 2043, ACS,APS,ACK,APK
      PRINT 2045,TAZ, TAZ1, TAZ2, TAZ3, TAZ4, TAZ5
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C     MEAN SERIES FOR NON-TIDAL CURRENT OR MEAN SEA LEVEL
      SUBROUTINE XMEAN(VIN,EH,N,BH,ISW)
C
C     IF ISW = 1 --- CALCULATE MEAN AND RETURN
C     IF ISW = 2 --- CALCULATE AND DEMEAN THE DATA
C
      DIMENSION VIN(N)
   15 FORMAT(/ 1X,14HSum of series ,3X,F12.3 / 1X, 7HDivisor,12X,I10
     1 / 1X,14HMean of series,5X,F10.5 )
      MA = N
      EH = 0.0
      DO 10 NA = 1,N
   10 EH = EH + VIN(NA)
      BH = EH/FLOAT(MA)
      PRINT 15, EH, MA, BH
C
      IF (ISW.EQ.1)   RETURN
C
C
      CALL DEMEAN (VIN,BH,N)
C
      RETURN
      END
      SUBROUTINE DEMEAN (VIN,BH,N)
C
C**** THIS SUROUTINE WILL DEMEAN DATA BEFORE PERFORMING H.A.
C
      DIMENSION VIN(N)
C
 1000 FORMAT (/' Data demeaned before performing H. A. '/)
 1001 FORMAT (1H0)
C
      WRITE (6,1000)
c     WRITE (6,1001)
C
      DO 1   I=1,N
C
      VIN(I) = VIN(I) - BH
C
    1 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C  DETERMINATION OF ORBITAL ELEMENTS FOR OBSERVED SERIES
      SUBROUTINE ASTRO(DAYB,DAYM,GRBS,GRMS,JOBX)
      COMMON/COSTX/CXX(30),OEX(5)
      COMMON XYER,MONTH,DAY,STT
 2110 FORMAT(1H0, 38HERROR.....PARAMETER (I) EXCEEDS LIMITS )
 2111 FORMAT(/15X,1HR,10X,1HQ,15X,9HU OF M(2))
 2112 FORMAT( 8X, F10.3, 3X, F10.3, 9X, F10.4)
 2040 FORMAT(/10F10.3/10F10.3/10F10.3)
      CALL ORBIT(XCEN,XSX,XPX,XHX,XP1X,XNX,OEX,T,XYER,5)
      DO 30 NOE = 1,30
   30 CXX(NOE) = 0.0
      CPLEX = XCEN/400.0 + 0.0001
      DICF = CPLEX - AINT(CPLEX)
      XCET = XCEN + 1.0
      DOBY = 0.25*(XYER - XCET)
      AMI = XYER - XCEN
      AMIT = AINT(DOBY)
      FARM = 0.25*(AMI)
      FARX = FARM - AINT(FARM)
      IF(DICF.LT.0.001.AND.AMI.NE.0) AMIT = AMIT + 1.0
      IF(DICF.GT.0.001.AND.AMI.EQ.0) FARX = 1
      IF(FARX.EQ.0.0.AND.MONTH.EQ.3.AND.DAY.EQ.1) GO TO 32
      IF(FARX.EQ.0.0.AND.DAYB.LE.60.0) DAYB = DAYB - 1.0
      IF(FARX.GT.0.0) DAYB = DAYB - 1.0
      IF(FARX.EQ.0.0.AND.DAYM.LE.60.0) DAYM = DAYM - 1.0
      IF(FARX.GT.0.0) DAYM = DAYM - 1.0
      IF(DAYB.EQ.61.0.AND.DAY.EQ.29.0) DAYM = DAYM - 1.0
      IF(DAYB.EQ.61.0.AND.DAY.EQ.29.0) DAYB = DAYB - 1.0
   32 CXX(1) = XSX + 129.38482032*AMI + 13.176396768*(DAYB+AMIT) +
     * 0.549016532*GRBS
      CXX(2) = XPX + 40.66246584*AMI + 0.111404016*(DAYB + AMIT) + 0.004
     1641834*GRBS
      CXX(3) = XHX - 0.238724988*AMI + 0.9856473288*(DAYB + AMIT) + 0.04
     110686387*GRBS
      CXX(4) = XP1X + 0.01717836*AMI + 0.000047064*(DAYB + AMIT) + 0.000
     1001961*GRBS
      CXX(5) = XPX + 40.66246584*AMI + 0.111404016*(DAYM + AMIT) + 0.004
     1641834*GRMS
      CXX(6) = XNX - 19.328185764*AMI - 0.0529539336*(DAYM + AMIT) - 0.0
     1022064139*GRMS
      IF(JOBX.LE.15.AND.JOBX.GT. 0 ) GO TO 40
      CALL TWOPI(CXX,6)
      GO TO 41
   40 CXX(7) = XPX + 40.66246584*AMI + 0.111404016*(DAYM + AMIT) + 0.004
     1641834*GRBS
      CXX(8) = XNX - 19.328185764*AMI - 0.0529539336*(DAYM + AMIT) - 0.0
     1022064139*GRBS
      CALL TWOPI(CXX, 8)
   41 AN =  CXX(6)*0.0174533
      AX =  CXX(6)
      EYE = 0.9136949 - 0.0356926*COS(AN)
      CXX(9) = ACOS(EYE)*57.2957795
      IF(CXX(9).LT.17.0.OR.CXX(9).GT.30.0) PRINT 2110
      CIG =  CXX(9)*0.0174533
      IF(CIG.EQ.0.0) GO TO 230
      IF(AX.EQ.0.0.OR.AX.EQ.180.0) GO TO 230
      VXX = 0.0896831*SIN(AN)/SIN(CIG)
      CXX(10) = ASIN(VXX)*57.2957795
      IF(AX.GT.180.0.AND. CXX(10).GT.0.0) CXX(10) = -1.0*CXX(10)
      CVX = CXX(10)*0.0174533
      EXX = 0.2067274*SIN(AN)*(1.0 - 0.0194926*COS(AN))
      IF(EXX.EQ.0.0) GO TO 202
      EZZ = 0.9979852 + 0.2067274*COS(AN) - 0.0020148*COS(2.0*AN)
      IF(EZZ.EQ.0.0) GO TO 202
      EXEZ  = EXX/EZZ
      IF(EXEZ.GT.3450.0) GO TO 202
      CXX(11) = ATAN(EXEZ)*57.2957795
      IF(AX.GT.180.0.AND.CXX(11).GT.0.0) CXX(11) = -1.0*CXX(11)
      CEX = CXX(11)*0.0174533
  202 VPX = SIN(CVX)/(COS(CVX) + 0.334766/SIN(2.0*CIG))
      CXX(12) = ATAN(VPX)*57.2957795
      IF(AX.GT.180.0.AND.CXX(12).GT.0.0) CXX(12) = -1.0*CXX(12)
      PVC = CXX(12)*0.0174533
      VPY = SIN(2.0*CVX)/(COS(2.0*CVX) + 0.0726184/SIN(CIG)**2)
      CXX(13) = ATAN(VPY)*57.2957795
      IF(AX.GT.180.0.AND.CXX(13).GT.0.0) CXX(13) = -1.0*CXX(13)
      PVCP = CXX(13)*0.0174533
  230 PGX = CXX(5) - CXX(11)
      CALL TWOPI(PGX, 1)
      XPG = PGX*0.0174533
      CXX(14) = PGX
      RAX = SIN(2.0*XPG)/(COS(0.5*CIG)**2/(6.0*SIN(0.5*CIG)**2) - COS(2.
     10*XPG))
      RXX = ATAN(RAX)*57.2957795
      UM2 = 2.0*(CXX(11) - CXX(10))
      UL2 =  UM2 - RXX
      IF(UL2.GT.180.0) UL2 = UL2 - 180.0
      IF(UL2.LT.180.0) UL2 = UL2 + 180.0
      CXX(15) = UL2
      ZES = ( 5.0*COS(CIG) - 1.0)*SIN(XPG)
      ZEC = (7.0*COS(CIG) + 1.0)*COS(XPG)
      CALL FITAN(ZES,ZEC,QXX,2)
      CRAV   =  0.5*( UM2)  + QXX  + 90.0
      CALL TWOPI(CRAV,1)
      CXX(16) = CRAV
      IF(JOBX.GT.15.AND.JOBX.GT.0 ) GO TO 88
      PGXX =  CXX(7) - CXX(11)
      CALL TWOPI(PGXX,1)
      XXPG = PGXX*0.0174533
      CXX(17) = PGXX
      BATX = CXX(8)*0.0174533
      EYEX = 0.9136949 - 0.0356926*COS(BATX)
      CXX(20) = ACOS(EYEX)*57.2957795
      UM2X = 2.0*(CXX(11) - CXX(10))
      EYIT = CXX(20)*0.0174533
      ZEXS = (5.0*COS(EYIT) - 1.0)*SIN(XXPG)
      ZEXC = ( 7.0*COS(EYIT) + 1.0)*COS(XXPG)
      CALL FITAN( ZEXS,ZEXC,QXXX,2)
      CXX(19) =  0.5*(UM2X) + QXXX + 90.0
      CRAVX = CXX(19)
      CALL TWOPI(CRAVX,1)
      CXX(19) = CRAVX
   88 DO 33 NT = 1,30
      DROP = CXX(NT)*10000. + 5.0
   33 CXX(NT) = AINT(DROP)*0.0001
 2141 FORMAT(10X  /  5X, 6HDAYB =, F6.0, 5X, 6HDAYM =, F6.0 ) 
 2142 FORMAT(        5X, 6HGRBS =, F7.2, 5X, 6HGRMS =, F7.2 )
      PRINT 2141,DAYB + 1,DAYM + 1
      PRINT 2142,GRBS,GRMS
      PRINT 2040,(CXX(IAX),IAX = 1,30)
      PRINT 2111
      PRINT 2112, RXX,QXX,UM2
      RETURN
      END
C
C ----------------------------------------------------------------------
C
C  DETERMINATION OF ORBITAL ELEMENTS FOR BEGINNING OF CENTURY
      SUBROUTINE ORBIT(XCEN,XSX,XPX,XHX,XP1X,XNX,OEX,T,XYER,NNN)
      DIMENSION OEX(NNN)
    7 FORMAT(/29H ***  Astronomical  Data  ***  )
    8 FORMAT(5F12.4, 5X, F5.0, 2F6.1)
      S = 13.1763968
      P = 0.1114040
      XH = 0.9856473
      P1 = 0.0000471
      XN = -0.0529539
      XCAN = XYER*0.01 + 0.001
      XCEN = AINT(XCAN)*100.0
      T = -3.0
      YR = 2.5
      GAT = 1600.0
      DO 10 JK = 1,30
      GP = (GAT)*0.01/4.0 + 0.00001
      FPX = ABS(GP)
      FP = AINT(FPX)
      COL = FPX - FP
      IF(COL.LT.0.010) GO TO 11
      IF(GAT.EQ.XCEN) GO TO 12
      YR = YR - 1.0
      GO TO 9
   11 IF( GAT.EQ.XCEN) GO TO 12
    9 GAT = GAT + 100.0
   10 CONTINUE
   12 T = (GAT - 1900.0)*0.01
      OEX(1) = 270.43742222 + 307.892*T + 0.002525*T**2 + 0.00000189*T**
     13 + YR*S
      OEX(2) = 334.32801944 + 109.0322055*T - 0.01034444*T**2 - 0.000012
     15*T**3 + YR*P
      OEX(3) = 279.69667778 + 0.768925*T + 0.0003025*T**2 + YR*XH
      OEX(4) = 281.22083333 + 1.719175*T + 0.00045278*T**2 + 0.00000333*
     1T**3 + YR*P1
      OEX(5) = 259.18253333 - 134.1423972*T + 0.00210556*T**2 + 0.000002
     122*T**3 + YR*XN
      CALL TWOPI(OEX,5)
      XSX = OEX(1)
      XPX = OEX(2)
      XHX = OEX(3)
      XP1X = OEX(4)
      XNX = OEX(5)
      PRINT 7
      PRINT 8, (OEX(JX),JX = 1,5) , GAT, T, YR
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE DAYXX(MONTH,DAY,STT,TM,DAYB,DAYM,GRBS,GRMS,CP)
 8800 FORMAT (/5X,8H DAYB = ,F5.0,10X,8H DAYM = ,F5.0/ 5X,8H GHBS = ,
     1 F7.2,10X,8H GHMS = ,F7.2/5X,17H ORIGINAL T.M. = ,F7.2,
     210X,18H CORRECTED T.M. = ,F7.2)
      GO TO (207,209,211,213,215,217,219,251,223,225,227,229),MONTH
  207 DAYB = DAY 
      GO TO 230 
  209 DAYB = DAY + 31.                                             
      GO TO 230
  211 DAYB = DAY + 59.
      GO TO 230  
  213 DAYB = DAY + 90. 
      GO TO 230   
  215 DAYB = DAY + 120.
      GO TO 230 
  217 DAYB = DAY + 151. 
      GO TO 230 
  219 DAYB = DAY + 181. 
      GO TO 230
  251 DAYB = DAY + 212.
      GO TO 230 
  223 DAYB = DAY + 243. 
      GO TO 230  
  225 DAYB = DAY + 273.
      GO TO 230 
  227 DAYB = DAY + 304.
      GO TO 230  
  229 DAYB = DAY + 334.
  230 IF(TM.LT.0.0.AND.STT.GT.12.0) GO TO 231
      IF(STT.GT.12.0) GO TO 233
  231 CP = TM + 15.0*STT
      IF(CP.LE.180.0) GO TO 235
  233 CP = TM - (24.00 - STT)*15.0 
      DAYB = DAYB + 1. 
  235 GRBS = CP/15.0  
      GRMS = GRBS + 12.00 
      IF(GRMS.GT.24.0) GO TO 240 
      DAYM = DAYB + 7.0 
      GO TO 242 
  240 GRMS = GRMS - 24.00
      DAYM = DAYB + 8.0
  242 CONTINUE
c 242 PRINT 8800,DAYB,DAYM,GRBS,GRMS,TM,CP
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE DBLOCK(VIN,DIN,N,IEL)
      REAL*4 DATAB(14,50)
      DIMENSION VIN(N),DIN(N),JUL(12),IAN(12)
      DATA(IAN(I),I=1,12)/31,59,90,120,151,181,212,243,273,304,334,365/
      DO 30 I=1,12
      JUL(I) = IAN(I)
   30 CONTINUE
      NP = N
      I = 0
      NP = NP + 50
      DO 102 L = 1,NP,50
      MPT = L + 49
      CALL NTRAN (11,2,700,DATAB,IS,22)
      IF (IS .NE. -2) GO TO 50
      GO TO 103
   50 LYNE = 0
      DO 100 MX = L,MPT
      LYNE = LYNE + 1
      I = I + 1
      IF(I.NE.1) GO TO 90
      XYER = DATAB(14,LYNE)
C
C**** TEST FOR LEAP YEAR
C
      IXYER = XYER/4.
      ITEST = 4*IXYER
      TEST = ITEST
      IF (TEST.EQ.XYER)   GO TO 70
      GO TO 71
   70 DO 60 IX=2,12
      JUL(IX) = JUL(IX) + 1
   60 CONTINUE
   71 STT = DATAB(12,LYNE) + DATAB(11,LYNE)/60.
      DO 91 IX=1,12
      IF(DATAB(13,LYNE).LE.JUL(IX)) GO TO 89
   91 CONTINUE
   89 MONTH = IX
      IF(IX.EQ.1) GO TO 92
      DAY = DATAB(13,LYNE) - JUL(IX-1)
      GO TO 90
   92 DAY = DATAB(13,LYNE)
   90 IF (IEL.EQ.0)   GO TO 93
      VIN(I) = DATAB(IEL,LYNE)
      GO TO 94
   93 VIN(I) = DATAB(9,LYNE)
      DIN(I) = DATAB(8,LYNE)
   94 IF (I.EQ.N)   GO TO 103
  100 CONTINUE
  102 CONTINUE
  103 PRINT 200, XYER,MONTH,DAY,STT
  200 FORMAT(' First point in CDF file  '/' Year  Month   Day   Time'/
     #F5.0,2X,I3,4X,F4.0,2X,F6.2)
      N = I
  104 RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE LOCATE(N,IHD1,IHD2,IHD3) 
      CHARACTER*10 IHD1,IHD2,IHD3
      DIMENSION IHEAD(34)
  400 CALL NTRAN (N,2,34,IHEAD,IS,22)
      PRINT 21,(IHEAD(NZ),NZ=1,33)
   21 FORMAT(' SUMMARY OF HEADING INFORMATION FOR THIS FILE',//,
     1' DATA GENERATION DATE ',A4,A4,A2,/,' S T A T I O N --     '
     *  ,A4,A4,A2,//,
     2' METER SERIAL NUMBER ',I7,/,' METER REFERENCE NUMBER ',I4,//,
     3' LATITUDE,  DIRECTION = ',A1,', DEGREES = ',I3,', MINUTES = ',I2,
     4', SECONDS = ',I2,/,' LONGITUDE, DIRECTION = ',A1,', DEGREES = ', 
     5I3,', MINUTES = ',I2,', SECONDS = ',I2,//,' START TIME, ', 
     6'MINUTES = ',I2,', HOURS = ',I2,', DAYS = ',I3,', YEAR  = ',
     7I2,/,' STOP  TIME, MINUTES = ',I2,', HOURS = ',I2,', DAYS = ',
     8I3,', YEAR  = ',I2,//,  
     9' TIME ANCHOR DROPPED, MINUTES =',I3,', HOURS =',I3,', DAYS =',I4,
     1', YEARS =',I3,/,' TIME ANCHOR RELEASE, MINUTES =',I3,', HOURS =',
     2I3,', DAYS =',I4,', YEARS =',I3,//,' DEPTH OF METER =',I5,' FEET,
     3IF NEGATIVE, FEET FROM THE SURFACE',/,29X,
     4'IF POSITIVE, FEET FROM THE BOTTOM'//)

      write(IHD2,22) IHEAD(33)   
      write(ihd1,1923) ihead(4),ihead(5),ihead(6)
 1923 format(a4,a4,a2)
   22 FORMAT(' AT ',I5,' ') 
      IF(IHEAD(33).GT.0) IHD3='ABOVE BOT.'
      IF(IHEAD(33).LT.0) IHD3='BELOW MLLW' 
   23 CONTINUE  
      RETURN  
      END  
C
C ----------------------------------------------------------------------
C
C     DETERMINE ARCTANGENT AND PLACE IN CORRECT QUADRANT
      SUBROUTINE FITAN( AUS,AUC,RTA,JMAP)
      IF(AUC.EQ.0.0) GO TO 14
      BXG =  AUS/AUC
      GO TO 15
   14 IF(AUS)16,55,17
   16 RTA = 270.0
      GO TO 88
   17 RTA = 90.0
      GO  TO 88
   15 RTA = ATAN(BXG)*57.2957795
      GO TO (88,38),JMAP
   38 IF(AUS)32,32,34
   32 IF(AUC)33,55,37
   33 RTA = 180.0 + RTA
      GO TO 88
   34 IF(AUC)35,17,88
   35 RTA = RTA + 180.0
      GO TO 88
   37 RTA = RTA + 360.0
      GO TO 88
   55 RTA = 0.0
   88 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE TWOPI( AUG, IO)
      DIMENSION  AUG(IO)
      DO 114 MO = 1,IO
  120 IF(AUG(MO))112,114,113
  112 AUG(MO) = AUG(MO)+ 360.0
      GO TO 120
  113 IF(AUG(MO) - 360.0)114,114,115
  115 AUG(MO) = AUG(MO) - 360.0
      GO TO 120
  114 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE ONE80(W,D,GONL,EQUIB,COVU)
      DIMENSION W(24),D(24),COVU(27)
    5 FORMAT(1H0,13X,'AMPLITUDE',37X,'EPOCH' / 14X,'RELATIONS',35X,
     1 'RELATIONS' ///  7X,'M(4)/M(2)', 7X, F7.3,20X,'2M(2)-M(4)', 5X, 
     2 F7.2  // 7X,'M(6)/M(2)', 7X, F7.3,20X,'3M(2)-M(6)', 5X,F7.2 //  
     3 7X,'S(2)/M(2)', 7X, F7.3,20X,'S(2)-M(2)', 6X,F7.2 //  7X,F
     4 'N(2)/M(2)', 7X, F7.3,20X,'M(2)-N(2)', 6X, F7.2 //  7X,
     5 'O(1)/K(1)', 7X, F7.3,20X,'K(1)-O(1)', 6X, F7.2 //  7X, 
     6 'K(1)+O(1)', 7X, F7.3,20X,'K(1)+O(1)', 6X, F7.2 //  6X,
     7 '(K(1)+O(1))/M(2) ', F7.3,20X,'M(2)-K(1)-O(1) ', F7.2 // 50X,
     8 'MKO',12X,F7.2)
    4 FORMAT(1H0,23X,'HARMONIC CONSTANT REDUCTIONS' / )
    6 FORMAT( //  18X,'AGE OF INEQUALITIES IN HOURS' //   3X,'PHASE AGE
     1 ', F6.1, 5X,'PARALLAX AGE ', F6.1, 5X,'DIURNAL AGE ', F6.1 //  )
    7 FORMAT(22X,'MEAN TIDE LEVEL = ', F6.3 // )
    8 FORMAT(15X,'ACC. IN HW DUE TO M(4) AND M(6)  =', F8.3 //
     1       15X,'ACC. IN LW DUE TO M(4) AND M(6)  =', F8.3 /)
    9 FORMAT(18X, 4F8.3)
   10 FORMAT(//4X,'GREENWICH HWI  ', F7.3,15X,'HW AMPLITUDE ', F8.3
     1       //4X,'GREENWICH LWI  ', F7.3,15X,'LW AMPLITUDE ', F8.3)
   11 FORMAT(//5X,'SPRING RANGE  ', F8.3,11X,'PERIGEAN RANGE  ', F8.3
     1       //6X,'NEAP  RANGE  ',  F8.3,11X,'APOGEAN  RANGE  ', F8.3)
   12 FORMAT(//27X,'TROPIC INTERVALS' // 11X,'TCHHWI ', F7.3,' (A)',
     1 14X, 'TCLLWI  ', F7.3, ' (A)'  // 11X,'TCLHWI  ', F7.3, ' (A)',
     2 14X, 'TCHLWI  ', F7.3, ' (A)')
   23 FORMAT(/// 12X,'SEQUENCE OF TIDE')
   14 FORMAT(1H+,29X,'    ( HHW    LW   LHW    LW )')
   15 FORMAT(1H+,29X,'    ( HHW   LLW   LHW   HLW )')
   16 FORMAT(1H+,29X,'    (  HW   LLW    HW   HLW )')
   17 FORMAT(1H+,29X,'    ( LHW   LLW   HHW   HLW )')
   18 FORMAT(1H+,29X,'    ( LHW    LW   HHW    LW )')
   19 FORMAT(1H+,29X,'    ( LHW   HLW   HHW   LLW )')
   20 FORMAT(1H+,29X,'    (  HW   HLW    HW   LLW )')
   21 FORMAT(1H+,29X,'    ( HHW   HLW   LHW   LLW )')
   22 FORMAT(///26X,'TROPIC    HEIGHTS' //12X, 'TCHHW  ', F8.3,18X,
     1 'TCLLW  ', F8.3 //12X, 'TCLHW  ', F8.3, 18X, 'TCHLW  ', F8.3 )
   24 FORMAT(///28X,'MEAN   HEIGHTS' //13X, 'MHHW  ', F8.3,20X,
     1 'MLLW  ', F8.3 //13X, 'MLHW  ', F8.3,20X, 'MHLW  ', F8.3 )
   25 FORMAT(//14X, 'DHQ  ', F7.3,21X, 'DLQ  ', F7.3 )
   26 FORMAT(///28X,'DIURNAL   TIDE'  //11X, 'TCHHWI  ', F8.3,18X,
     1 'TCLLWI  ', F8.3  //12X, 'TCHHW  ', F8.3,19X, 'TCLLW  ', F8.3//
     2 13X, 'MHHW  ', F8.3,19X, 'MLLW  ', F8.3 )
   29 FORMAT(/// ' THE COMPOUND WAVE IS SIX-DIURNAL FOR ALL PHASE REL.')
   30 FORMAT(/// ' DUE TO M(6), DOUBLE AMPLITUDES HAVE OCCURED ' // 
     1 ' SEE MANUAL OF HARMONIC CONSTANT REDUCTION  ( PAGES 38-39 ) ' )
   31 FORMAT( /// ' COMPOUND WAVE IS QUARTER-DIURNAL FOR ALL PHASE ',
     1 'RELATIONS ' )
   32 FORMAT(/// ' A DOUBLE LOW AMPLITUDE IS PRESENT ' )
   33 FORMAT(/// ' A DOUBLE HIGH AMPLITUDE IS PRESENT ' )
   34 FORMAT(/// ' DUE TO M(4), THE WAVE IS APPROACHING CRITICAL ',
     1 'LIMITS FOR DOUBLE AMPLITUDES TO OCCUR ' )
   37 FORMAT(/// ' DOUBLE HIGH AMPLITUDES OF EQUAL HEIGHTS AND ',
     1 'DOUBLE LOW AMPLITUDES OF EQUAL HEIGHTS ARE PRESENT ' )
   38 FORMAT(/// ' DUE TO M(6), THE COMPOUND WAVE IS SIX-DIURNAL ',
     1 'FOR ALLPHASE RELATIONS ' )
   39 FORMAT(// 10F10.4)
   40 FORMAT( 12X, 4F10.4)
   41 FORMAT(12X, 2F10.4)
   42 FORMAT(/// 22X,12HMEAN RANGE   , F8.4)
   43 FORMAT(////// 1H0,'**** COMPUTE MEAN HEIGHTS MANUALLY ***',
     1 'MEAN HEIGHT SUBROUTINE HAS TO BE DEBUGGED ****')

      HK1 = D(2)
      HO1 = D(12)
      HP1 = D(14)
      HS2 = D(18)
      HM2 = D(6)
      HM4 = D(7)
      HM6 = D(8)
      HN2 = D(10)
      GK1 = W(2)
      GO1 = W(12)
      GP1 = W(14)
      GS2 = W(18)
      GM2 = W(6)
      GM4 = W(7)
      GM6 = W(8)
      GN2 = W(10)
      IF( HM2.EQ.0.0) GO TO 888
      IF(HK1.EQ.0.0) GO TO 888
      PT = GM2 - GM4
      IF(ABS(PT).LT.180.0) GO TO 102
      IF(PT)101,102,103
  101 RM2 = GM2 + 360.0
      RM4 = GM4
      GO TO 104
  103 RM4 = GM4 + 360.0
      RM2 = GM2
      GO TO 104
  102 RM2 = GM2
      RM4 = GM4
  104 P4  = 2.0*RM2 - RM4
      CALL TWOPI(P4,1)
      E4 = P4
      P4 = P4*0.0174533
      PXT = GM2 - GM6
      IF(ABS(PXT).LT.180.0) GO TO 106
      IF(PXT)105,106,107
  105 RM2 = GM2 + 360.0
      RM6 = GM6
      GO TO 108
  107 RM6 = GM6 + 360.0
      RM2 = GM2
      GO TO 108
  106 RM2 = GM2
      RM6 = GM6
  108 P6 = 3.0*RM2 - RM6
      CALL TWOPI(P6,1)
      E6 = P6
      P6 = P6*0.0174533
      R4 = HM4/HM2
      R6 = HM6/HM2
      ASM = HS2/HM2
      ANM = HN2/HM2
      AOK = HO1/HK1
      AOKK = HK1 + HO1
      AOKM = AOKK/HM2
      ESM = GS2 - GM2
      CALL ONEPI(ESM)
      EMN = GM2 - GN2
      CALL  ONEPI(EMN)
      RGKO = GK1 + GO1
      EKO  = GK1 - GO1
      IF(ABS(EKO).GT.180.0) RGKO = RGKO + 360.0
      CALL ONEPI(EKO)
      IF(RGKO.GT.720.0) RGKO = RGKO - 720.0
      EMKO = GM2 - RGKO
      IF(EMKO.LT.0.0) EMKO = EMKO + 720.0
      RMKO  =  EMKO/2.0
      TRMKO = RMKO
      CALL TWOPI(RMKO,1)
      CALL TWOPI(RMKO,1)
      TEMKO = EMKO*0.0174533
      XMTL = HM4*COS(P4) - 0.03*AOKK*COS(TEMKO)*AOKM
      PHAGE = 0.984*ESM
      PARAGE = 1.837*EMN
      DIAGE  = 0.911*EKO
      PRINT 4
      PRINT 5,  R4, E4, R6, E6, ASM, ESM, ANM, EMN, AOK, EKO, AOKK, RGKO
     1,AOKM,EMKO,RMKO
      PRINT  6, PHAGE, PARAGE, DIAGE
      PRINT  7, XMTL
      IF( AOKM.GT.4.00) GO TO 777
      CALL RELATE(E4,180.0,PP4)
      CALL RELATE(E4,180.0,ET)
      IF(R4.LE.0.250) GO TO 121
      PRINT 34
      IF(R4.GT.0.50) PRINT 31
      IF(E4.EQ.000.0.OR.E4.EQ.360.0) PRINT 32
      IF(E4.EQ.180.0) PRINT 33
  121 CALL EFFECT(VP,R4,2.0,E4,1.0,0.25,6.2103006,1,1.0,1.0)
      CALL EFFECT(WP,R4,2.0,ET,1.0,0.25,6.2103006,1,1.0,1.0)
      XVP = VP*0.0174533
      XWP = WP*0.0174533
      BETA = ( E6 - 3.0*VP)*0.0174533
      ALFA = ( E6 - 3.0*WP)*0.0174533
      CALL RELATE(E6,180.0,ER)
      IF(R6.LT.0.111) GO TO 129
      IF(E6.GT.090.0.AND.E6.LT.270.0) PRINT 30
      IF(R6.GE.0.3333) PRINT 38
      IF(R6.GT.0.1111.AND.E6.EQ.180.0) PRINT 37
      IF(R6.GT.0.1111.AND.E6.EQ.180.0) IDUB2 = 1
  129 CONTINUE
      IF(R6.LE.0.3333) GO TO 125
      CALL EFFECT(VZ,R6,3.0,E6,1.0,0.11,6.2103006,1,1.0,1.0)
      CALL EFFECT(VY,R6,3.0,ER,1.0,0.11,6.2103006,1,1.0,1.0)
      VPP = VZ
      WPP = VY
      GO TO 130
  125 CONTINUE
      GIV = 3.0*R6*SIN(BETA)
      GIW = 3.0*R6*SIN(ALFA)
      VPP = ASIN(GIV)*57.29578
      WPP = ASIN(GIW)*57.29578
  130 CONTINUE
      V = VP + VPP
      XW = WP + WPP
      PRINT 8, V,XW
      PRINT 9, VP,VPP,WP,WPP
      GHWI = 0.0345*(GM2 - V) + 0.069*GONL
      IF(GHWI.GE.12.42) GHWI = GHWI - 12.42
      SI = GM2 -XW
      CALL RELATE(SI,180.0,SI)
  131 GLWI = 0.0345*SI + 0.069*GONL
      IF(GLWI.GE.12.42) GLWI = GLWI - 12.42
      RV = V*0.0174533
      RW =XW*0.0174533
      ZX =  P4 - 2.0*RV
      ZY =  P4 - 2.0*RW
      ZA =  P6 - 3.0*RV
      ZB =  P6 - 3.0*RW
      RANGE = 1.02*HM2*(COS(RV) + COS(RW) + (0.020 + 0.577*ASM**2) + 0.0
     172*AOKM**2) + 1.02*HM4*(COS(ZX) - COS(ZY)) + 1.02*HM6*(COS(ZA) + C
     2OS(ZB))
      DIAR = 0.5*(0.020 + 0.577*ASM**2 + 0.072*AOKM**2)
      XHW = 1.02*HM2*( COS(RV) + DIAR ) + 1.02*HM4*COS(ZX) + 1.02*HM6*CO
     1S(ZA)
      XLW = -1.02*HM2*(COS(RW) + DIAR ) + 1.02*HM4*COS(ZY) - 1.02*HM6*CO
     1S(ZB)
      PRINT 10,GHWI, XHW, GLWI, XLW
      PRINT 42, RANGE
      IF(AOKM.GT.1.50) GO TO 666
      SPRG = RANGE - 0.536*(HS2**2/HM2) + HS2*(1.96 - 0.08*AOKM**2)
      PAEN = RANGE - 0.536*(HS2**2/HM2) - HS2*(1.96 - 0.08*AOKM**2)
      PERIR = ( 1.00 + ANM)*RANGE
      APOR = ( 1.00 - 0.75*ANM)*RANGE
      PRINT 11, SPRG, PERIR, PAEN, APOR
      IF(AOKM.LT.0.25) GO TO 888
C     COMPUTE TROPIC INTERVALS AND HEIGHTS
  666 AMKO = RMKO - 0.5*V
      CALL TWOPI(AMKO,1)
      BMKO = RMKO - 0.5*XW
      CALL RELATE(BMKO,90.0,BMKO)
      CALL TWOPI(BMKO,1)
      IF(AOKM.LT.2.0) GO TO 133
      TEA = RMKO*0.0174533
  133 CALL DEFECT( 6.2103006,AP ,AOKM,AMKO, 2.0, 1, 1,IWAVE,1.0)
      CALL DEFECT( 6.2103006,APP,AOKM,AMKO,-2.0, 1, 3,JWAVE,-1.0)
      CALL DEFECT( 6.2103006,BP ,AOKM,BMKO, 2.0, 1, 2,LWAVE,1.0)
      CALL DEFECT( 6.2103006,BPP,AOKM,BMKO,-2.0, 1, 4,KWAVE,-1.0)
      IF(IWAVE.EQ.1.AND.LWAVE.EQ.1) GO TO 777
      TCHHWI = GHWI - AP
      TCLLWI = GLWI - BP
      TCLHWI = GHWI - APP
      TCHLWI = GLWI - BPP
      IF(AMKO.GE.000.0.AND.AMKO.LE.090.0) GO TO 135
      IF(AMKO.GE.270.0.AND.AMKO.LE.360.0) GO TO 135
  134 TCHHWI = TCHHWI + 12.42
      GO TO 137
  135 GESS = GM2*0.0345 + 0.069*GONL
      ESS = ABS(GESS - TCHHWI)
      IF(ESS.LE.2.00) GO TO 136
      GO TO 134
  136 TCLHWI = TCLHWI + 12.42
  137 IF(BMKO.GE.000.0.AND.BMKO.LE.090.0) GO TO 140
      IF(BMKO.GE.270.0.AND.BMKO.LE.360.0) GO TO 140
  138 TCLLWI = TCLLWI + 12.42
      GO TO 142
  140 PROX = (GM2 + 180.0)/28.9841042 + 0.069*GONL
      ROX = ABS(PROX - TCLLWI)
      IF(ROX.LE.2.00) GO TO 141
      GO TO 138
  141 TCHLWI = TCHLWI + 12.42
  142 CONTINUE
      IF(JWAVE.EQ.1) TCLHWI = 0.0
      IF(KWAVE.EQ.1) TCHLWI = 0.0
C     DETERMINE TIDE SEQUENCE
      DMKO = RMKO
      CALL TWOPI(DMKO,1)
      IF( DMKO.EQ.000.0 ) IS = 1
      IF( DMKO.GT.000.0.AND.DMKO.LT.090.0 ) IS = 2
      IF( DMKO.EQ.090.0 ) IS = 3
      IF( DMKO.GT.090.0.AND.DMKO.LT.180.0 ) IS = 4
      IF( DMKO.EQ.180.0 ) IS = 5
      IF( DMKO.GT.180.0.AND.DMKO.LT.270.0 ) IS = 6
      IF( DMKO.EQ.270.0 ) IS = 7
      IF( DMKO.GT.270.0.AND.DMKO.LT.360.0 ) IS = 8
      CALL HITEF(AMKO,BMKO,AP,BP,APP,BPP,QA,QB,QC,QD,TE,TF,AOKM,IWAVE,JW
     1AVE,KWAVE,LWAVE)
      PRINT 39,AP,BP,APP,BPP,QA,QB,QC,QD,TE,TF
      QC = -1.0*QC
      QD = -1.0*QD
      HP1 = 0.193
      CALL PEFECT(HP1,TE,T3,CORR1,JWAVE)
      CALL PEFECT(HP1,TF,T4,CORR2,KWAVE)
      PRINT 40,CORR1,CORR2, T3,T4
      T1 = QA + 0.5*CORR1
      T2 = QC + 0.5*CORR2
      TCHHW = HM2*T1
      TCLLW = -1.0*HM2*T2
      TCLHW = TCHHW - HM2*T3
      TCHLW = TCLLW + HM2*T4
      IF(JWAVE.EQ.1) TCLHW = 0.0
      IF(KWAVE.EQ.1) TCHLW = 0.0
  176 PRINT 12, TCHHWI, TCLLWI, TCLHWI, TCHLWI
      PRINT 22, TCHHW, TCLLW, TCLHW, TCHLW
      PRINT 23
      GO TO (181,182,183,184,185,186,187,188),IS
  181 PRINT 14
      GO TO 189
  182 PRINT 15
      GO TO 189
  183 PRINT 16
      GO TO 189
  184 PRINT 17
      GO TO 189
  185 PRINT 18
      GO TO 189
  186 PRINT 19
      GO TO 189
  187 PRINT 20
      GO TO 189
  188 PRINT 21
  189 CONTINUE
      PRINT 43
      GO TO 888
C     COMPUTE MEAN HEIGHTS
      CALL TAB67(HO1,HK1,EKO,CFO,XX,EQUIB)
      AOKR = AOK/0.646
      IF(AOKR.GT.1.0) AOKR = 1.0/AOKR
      CAMKO = AMKO + XX
      CBMKO = BMKO + XX
      CALL TWOPI(CAMKO,1)
      CALL TWOPI(CBMKO,1)
      GOKM = AOKM*CFO
      IF(GOKM.GT.1.0) GOKM = 1.0/GOKM
      CALL HITEF(CAMKO,CBMKO,AP,BP,APP,BPP,HA,HB,HC,HD,TC,TD,GOKM,IWAVE,
     1JWAVE,KWAVE,LWAVE)
      PRINT 39,AP,BP,APP,BPP,HA,HB,HC,HD,TC,TD
      IF(JWAVE.EQ.1) TC = 0.0
      IF(KWAVE.EQ.1) TD = 0.0
      CALL PEFECT(HP1,TC,TP,COREC1,JWAVE)
      CALL PEFECT(HP1,TD,TQ,COREC2,KWAVE)
      PRINT 40,COREC1,COREC2,TP,TQ
      COKM = CFO*AOK
      IF(COKM.GT.1.0) COKM = 1.0/COKM
      CALL HITEF(CAMKO,CBMKO,AP,BP,APP,BPP,HA,HB,HC,HD,TA,TB,COKM,IWAVE,
     1JWAVE,KWAVE,LWAVE)
      PRINT 41, CFO,XX
      PRINT 39,AP,BP,APP,BPP,HA,HB,HC,HD,TA,TB
      HC = -1.0*HC
      HD = -1.0*HD
      IF(KWAVE.EQ.1) TB = 0.0
      IF(JWAVE.EQ.1) TA = 0.0
      CALL PEFECT(HP1,TA,TR,COREC1,JWAVE)
      CALL PEFECT(HP1,TB,TS,COREC2,KWAVE)
      PRINT 40, COREC1,COREC2,TR,TS
      COKM = AOK
      IF(COKM.GT.1.0) COKM = 1.0/COKM
      CALL HITEF( AMKO, BMKO,AP,BP,APP,BPP,HA,HB,HC,HD,TG,TF,COKM,IWAVE,
     1JWAVE,KWAVE,LWAVE)
      PRINT 39,AP,BP,APP,BPP,HA,HB,HC,HD,TG,TF
      IF(JWAVE.EQ.1) TG = 0.0
      IF(KWAVE.EQ.1) TF = 0.0
      CALL PEFECT(HP1,TG,TT,COREC1,JWAVE)
      CALL PEFECT(HP1,TF,TU,COREC2,KWAVE)
      PRINT 40, COREC1,COREC2,TT,TU
      TX = 0.646 - ABS(TT - TR)
      TY = 0.646 - ABS(TU - TS)
      DHQ = AOKK*TX
      DLQ = AOKK*TY
      XMHHW =  0.5*RANGE + XMTL + DHQ
      XMLLW = -0.5*RANGE + XMTL - DLQ
      XMLHW = XMHHW - 2.0*DHQ
      XMHLW = XMLLW + 2.0*DLQ
      PRINT 24,XMHHW, XMLLW, XMLHW, XMHLW
      PRINT 25, DHQ, DLQ
      GO TO 888
C       COMPUTE DIURNAL TIDE
  777 RECIP = HM2/AOKK
      AMKO = RMKO
      CALL RELATE(AMKO,090.0,BMKO)
      CALL EFFECT(DX ,RECIP,2.0,AMKO,-1.0,0.25,6.2103006,2,2.0,-1.0)
      CALL EFFECT(DXP,RECIP,2.0,BMKO, 1.0,0.25,6.2103006,2,2.0,-1.0)
      TCHHWI = (RGKO - DX)*0.0345
      TGKO = RGKO - DXP
      TCLLWI = TGKO*0.0345
      CALL RELATE(TCLLWI,12.42,TCLLWI)
      TD = DX*0.0174533
      TDD = DXP*0.0174533
      CALL TAB67(HO1,HK1,EKO,CFO,XX,EQUIB)
      COKM = RECIP*CFO
      CAMKO = AMKO + XX
      CBMKO = BMKO + XX
      PRINT 41, CFO,XX
      ALD = 2.0*(CAMKO + DX)
      CALL TWOPI(ALD,1)
      ALD = ALD*0.0174533
      ALDD = 2.0*(CBMKO + DXP)
      CALL TWOPI(ALDD,1)
      ALDD = ALDD*0.0174533
C     TABLE 17
      FACTH = COS(TD) + COKM*COS(ALD)
      FACTL =  1.0*COS(TDD) + COKM*COS(ALDD)
      XMLLW =-AOKK*FACTL
      BLD = 2.0*(AMKO + DX)
      XMHHW = AOKK*FACTH
      CALL TWOPI(BLD,1)
      BLD = BLD*0.0174533
      BLDD = 2.0*(BMKO + DXP)
      CALL TWOPI(BLDD,1)
      BLDD = BLDD*0.0174533
C         TABLE 16
      FACH = COS(TD) + RECIP*COS(BLD)
      FACL =  1.0*COS(TDD) + RECIP*COS(BLDD)
      TCHHW = AOKK*FACH
      TCLLW =-AOKK*FACL
      PRINT 40, FACTH,FACTL,FACH,FACL
      PRINT 26,TCHHWI, TCLLWI, TCHHW, TCLLW,XMHHW,XMLLW
  888 RETURN
      END
      SUBROUTINE PEFECT(HP1,TAA,T11T,CORR,MWAVE)
      IF(MWAVE.EQ.1) GO TO 9
      TLL = 2.0*HP1
      TVAL = 0.386*TAA
      AVL = TAA - TVAL
      VAL = ABS(AVL)
      IF(VAL.GE.TLL) GO TO 9
      XL = VAL/TLL
      IF(ABS(XL).GE.1.0) GO TO 9
      BERE = ACOS(XL)*57.29578
      T12 = 0.6366*(1.0 - XL**2)**0.5 - 0.0111*XL*BERE
      GO TO 11
   10 T12 = -0.0111*XL*BERE
      GO TO 11
    9 CONTINUE
      T12 = 0.0
      TLL = 0.0
   11 CORR =     T12
      T11T = TAA + CORR*TVAL
      RETURN
      END
      SUBROUTINE EFFECT(EVP,R,CONST,EPAR,SYGN,CRIT,ENTRE,IKE,FAC1,FAC2)
    5 FORMAT(//// '    ERROR IN COMPUTATION OF V PRIME   ',4F10.3)
    6 FORMAT(//// '    ERROR IN COMPUTATION OF W PRIME   ',4F10.3)
      EPAX = EPAR
      FACT = 1.0
      IF(EPAR.EQ.0.0.OR.EPAR.EQ.360.0) GO TO 100
      IF(R.EQ.0.0) GO TO 100
      GO TO (55,56),IKE
   55 IF(EPAR.EQ.180.0.AND.R.LE.CRIT) GO TO 100
      IF(EPAR.GT.000.0.AND.EPAR.LT.180.0) GO TO 65
      IF(EPAR.GT.180.0.AND.EPAR.LT.360.0) GO TO 66
      IF(EPAR.GE.180.0) GO TO 66
      GO TO 65
   56 IF(EPAR.EQ.180.0) GO TO 100
      IF(EPAR.EQ.090.0.AND.R.LE.CRIT) GO TO 100
      IF(EPAR.EQ.270.0.AND.R.LE.CRIT) GO TO 100
      IF(EPAR.EQ.090.0) GO TO 68
      IF(EPAR.EQ.270.0) GO TO 66
      IF(EPAR.GT.000.0.AND.EPAR.LT.090.0) GO TO 68
      IF(EPAR.GT.090.0.AND.EPAR.LT.180.0) GO TO 65
      IF(EPAR.GT.180.0.AND.EPAR.LT.270.0) GO TO 66
      IF(EPAR.GT.270.0.AND.EPAR.LT.360.0) GO TO 67
   65 ENTRI = 0.0
      XIT = ENTRE
      GO TO 69
   66 EPAR = 360.0 - EPAR
      FACT = -1.0
      GO TO 65
   67 EPAR = 360.0- EPAR
      FACT = 1.0
      GO TO 65
   68 FACT = -1.0
      GO TO 65
   69 CONTINUE
      EVP = ENTRI
      SOLN = R*CONST
   99 CONTINUE
      IF(EVP.LT.ENTRI.OR.EVP.GT.XIT) GO TO 111
      VP = EVP*14.4920521
      VV =  VP*0.0174533
      VX = FAC1*EPAX - FAC2*CONST*VP*FACT
      CALL TWOPI(VX,1)
      VX = VX*0.0174533
      IF(VX.LE.0.000087.OR.VX.GE.6.283101) GO TO 101
      IF(VX.GT.3.141507.AND.VX.LT.3.141681) GO TO 101
      VEX = SYGN*(SIN(VV)/SIN(VX))
      VT = ABS(VEX) - SOLN
      TEST1 = ABS(VT)
      TEST2 =.000180
      IF(EPAR.LE.030.0) TEST2 = .000500
      IF(EPAR.LT.015.0) TEST2 = .001000
      IF(TEST1.LT.TEST2 ) GO TO 101
   88 CONTINUE
      EVP = EVP +.00005
      GO TO 99
   89 EVP = EVP/CONST
      GO TO 101
  111 PRINT 5, EVP,R,EPAX,VEX
      STOP
  100 EVP = 0.0
  101 EVP = EVP*FACT*14.4920521
      EPAR = EPAX
      RETURN
      END
      SUBROUTINE TAB67(HO1,HK1,EKO,CFO,XX,EQUIB)
      PHI = EQUIB - 0.5*EKO
      CALL TWOPI(PHI,1)
      IF(PHI.GE.000.0.AND.PHI.LT.090.0) GO TO 200
      IF(PHI.GT.270.0.AND.PHI.LE.360.0) GO TO 200
      IF(PHI.EQ.090.0.OR.PHI.EQ.270.0) GO TO 202
      IPTRAN = 2
      GO TO 201
  200 IPTRAN = 1
  201 IF(IPTRAN.EQ.2) CALL RELATE(PHI,180.0,PHI)
  202 CONTINUE
      RPHI = PHI*0.0174533
      IF(HO1.LE.HK1) GO TO 203
      IF(HO1.GT.HK1) GO TO 204
  203 RX = HO1/HK1
      SYNE = 1.0
      GO TO 205
  204 RX = HK1/HO1
      SYNE = -1.0
  205 IF(PHI.EQ.000.0.OR.PHI.EQ.180.0) GO TO 208
      IF(PHI.EQ.360.0) GO TO 208
      CFO = (( 1.0 + RX**2 + 2.0*RX*COS(RPHI*2.0))**0.5)/(1.0 + RX)
      IF(PHI.EQ.090.0.OR.PHI.EQ.270.0) GO TO 209
      SX = ((1.0 - RX)*TAN(RPHI))/(1.0 + RX)
      XX = ATAN(SX)*SYNE*57.29578
      GO TO 206
  208 CFO = 1.0
      XX = 0.0
      GO TO 210
  209 XX = SYNE*90.0
  206 IF(PHI.GE.000.0.AND.PHI.LE.180.0) GO TO 207
      IF(PHI.GE.270.0.AND.PHI.LE.360.0) GO TO 207
      GO TO 210
  207 CONTINUE
  210 RETURN
      END
      SUBROUTINE  HITEF(AMKO,BMKO,AP,BP,APP,BPP,HA,HB,HC,HD,T1,T2,AOKM,I
     1WAVE,JWAVE,KWAVE,LWAVE)
      APX = AP*14.4920521
      BPX = BP*14.4920521
      APPX = APP*14.4920521
      BPPX = BPP*14.4920521
      PAP = (AMKO - APX)*0.0174533
      PBP = (BMKO - BPX)*0.0174533
      PAPP = (AMKO - APPX)*0.0174533
      PBPP = (BMKO - BPPX)*0.0174533
      APX = APX*0.0174533
      BPX = BPX*0.0174533
      APPX = APPX*0.0174533
      BPPX = BPPX*0.0174533
      IF(AMKO.GE.000.0.AND.AMKO.LE.090.0) GO TO 170
      IF(AMKO.GE.270.0.AND.AMKO.LE.360.0) GO TO 170
      IF(AMKO.GT.090.0.AND.AMKO.LT.270.0) GO TO 171
  170 HA = COS(2.0*APX) + AOKM*COS(PAP)
      IF(JWAVE.EQ.1) GO TO 50
      HB = COS(2.0*APPX) - AOKM*COS(PAPP)
      IF(AMKO.EQ.090.0.OR.AMKO.EQ.270.0) HB = HA
      T1 = HA - HB
      GO TO 172
  171 HA = COS(2.0*APX ) - AOKM*COS(PAP )
      IF(JWAVE.EQ.1) GO TO 50
      HB = COS(2.0*APPX)+ AOKM*COS(PAPP)
      IF(AMKO.EQ.090.0.OR.AMKO.EQ.270.0) HB = HA
      T1 = HA - HB
      GO TO 172
   50 HB = 0.0
      T1 = 0.0
  172 CONTINUE
      IF(BMKO.GE.000.0.AND.BMKO.LE.090.0) GO TO 174
      IF(BMKO.GE.270.0.AND.BMKO.LE.360.0) GO TO 174
      IF(BMKO.GT.090.0.AND.BMKO.LT.270.0) GO TO 175
  174 HC = -1.0*(COS(2.0*BPX ) + AOKM*COS(PBP ))
      IF(KWAVE.EQ.1) GO TO 60
      HD = -1.0*(COS(2.0*BPPX) - AOKM*COS(PBPP))
      IF(BMKO.EQ.090.0.OR.BMKO.EQ.270.0) HD = HC
      T2 = HD - HC
      GO TO 176
  175 HC = -1.0*(COS(2.0*BPX ) - AOKM*COS(PBP ))
      IF(KWAVE.EQ.1) GO TO 60
      HD = -1.0*(COS(2.0*BPPX) + AOKM*COS(PBPP))
      IF(BMKO.EQ.090.0.OR.BMKO.EQ.270.0) HD = HC
      T2 = HD - HC
      GO TO 176
   60 HD = 0.0
      T2 = 0.0
  176 CONTINUE
      IF(T1.EQ.-0.0) T1 = 0.0
      IF(T2.EQ.-0.0) T2 = 0.0
      RETURN
      END
      SUBROUTINE RELATE(AUG,TEST,RESULT)
      IF(AUG.GE.TEST ) GO TO 120
      RESULT = AUG + TEST
      GO TO 121
  120 RESULT = AUG - TEST
  121 RETURN
      END
      SUBROUTINE DEFECT(ENTRE,CP,R,EPAR,SYGN1,ISW,IPHASE,IWAVE,FATC)
    4 FORMAT(//// '  ERROR IN COMPUTATION OF EFFECTS  ', 4F10.5,I5)
    5 FORMAT(/   ' TIDE BECOMES DIURNAL BY MERGING OF THE TWO LOW WATERS
     1 OF EQUAL HEIGHT WITH THE LHW TO FORM SINGLE LOW WATER ' / ' P = '
     2, F7.3, ' R = ', F5.2, ' ARGUMENT = ', F7.3,I5)
    6 FORMAT(/   ' TIDE BECOMES DIURNAL BY MERGING OF THE TWO HIGH WATER
     1S OF EQUAL HEIGHT WITH THE HLW TO FORM SINGLE HIGH WATER ' / ' P =
     2 ', F7.3, ' R = ', F5.2, ' ARGUMENT = ', F7.3,I5)
    7 FORMAT(//  ' COMPOUND WAVE IS DIURNAL ',' P = ',F7.3,' R = ',F5.2,
     1' ARGUMENT = ',F7.3, I5 )
    8 FORMAT(//  ' WAVE  CHANGES FROM SEMI-DIURNAL TO DIURNAL ' )
      EPAX = EPAR
      IWAVE = 0
      FACT = 1.0*FATC
      IF(R.EQ.0.0) GO TO 87
      IF(R.LT.4.00) GO TO 77
      IF(EPAR.EQ.090.0.OR.EPAR.EQ.270.0) GO TO 85
   77 IF(EPAR.EQ.000.0.OR.EPAR.EQ.180.0) GO TO 87
      IF(EPAR.EQ.090.0) GO TO 65
      IF(EPAR.EQ.270.0) GO TO 66
      IF(EPAR.EQ.360.0) GO TO 87
      IF(EPAR.GT.000.0.AND.EPAR.LT.090.0) GO TO 65
      IF(EPAR.GT.090.0.AND.EPAR.LT.180.0) GO TO 67
      IF(EPAR.GT.180.0.AND.EPAR.LT.270.0) GO TO 68
      IF(EPAR.GT.270.0.AND.EPAR.LT.360.0) GO TO 66
   65 ENTRI = 0.0
      XIT = ENTRE
      GO TO 69
   66 EPAR = 360.0 - EPAR
      FACT = -1.0*FATC
      GO TO 65
   67 EPAR = 180.0 - EPAR
      FACT = -1.0*FATC
      GO TO 65
   68 EPAR = EPAR - 180.0
      FACT = 1.0*FATC
      GO TO 65
   69 CONTINUE
      IF(R.LE.2.0) GO TO 76
      GO TO (76,76,70,70),IPHASE
   70 IF(EPAX.EQ.090.0.OR.EPAX.EQ.270.0) GO TO 76
      PR = EPAR*0.0174533
      TN = TAN(PR)**(1./3.)
      AT = ATAN(TN)
      AR = AT*57.29578
      CP = AR*0.06900334
      QUOT = EPAR + AR
      TOUQ = QUOT*0.0174533
      IF(QUOT.EQ.000.0.OR.QUOT.EQ.180.0) GO TO 90
      GRX = (-2.0*SIN(2.0*AT))/SIN(TOUQ)
      RXX = ABS(GRX)*100.0 + 0.5
      RX = AINT(RXX)*0.01
      GO TO 91
   90 RX = 999.0
   91 IF(R.GE.RX) GO TO 80
   76 CONTINUE
      CP = ENTRI
      SOLN = R
      P = EPAX
   78 CONTINUE
      IF(CP.LT.ENTRI.OR.CP.GT.XIT) GO TO 112
      CV = CP*14.4920521
      CVV = 2.0*CV*FACT
      CALL TWOPI(CVV,1)
      CVV = CVV*0.0174533
      CX = P - CV*FACT
      CALL TWOPI(CX,1)
      GO TO (71,72),ISW
   71 IF(CX.EQ.000.0.OR.CX.EQ.180.0) GO TO 79
      IF(CX.EQ.360.0) GO TO 79
   99 CX = CX*0.0174533
      BETA = SIN(CX)
      GO TO 73
   72 IF(CX.EQ.090.0.OR.CX.EQ.270.0) GO TO 79
      IF(EPAX.EQ.090.0.OR.EPAX.EQ.270.0) GO TO 99
      CX = CX*0.0174533
      BETA = COS(CX)
   73 CEX = SYGN1*(SIN(CVV)/BETA)
      CT = ABS(CEX) - R
      TEST1 = ABS(CT)
      TEST2 = 0.0005
      IF(EPAX.LT.015.0) TEST2 = 0.0010
      IF(TEST1.LT.TEST2) GO TO 116
   79 CONTINUE
      CP = CP + 0.00005
      GO TO 78
   80 PRINT 8
      IF(EPAX.EQ.090.0.OR.EPAX.EQ.270.0) GO TO 81
      IF(EPAX.EQ.000.0.OR.EPAX.EQ.180.0) GO TO 82
      IF(EPAX.EQ.045.0.OR.EPAX.EQ.135.0) GO TO 83
      IF(EPAX.EQ.225.0.OR.EPAX.EQ.315.0) GO TO 83
      GO TO (83,83,82,81),IPHASE
   81 PRINT 6, EPAX, R, AR, IPHASE
      CP = 0.0
      GO TO 84
   82 PRINT 5, EPAX, R, AR, IPHASE
      CP = 0.0
      GO TO 84
   83 PRINT 7, EPAX, R, AR, IPHASE
   84 IWAVE = 1
      GO TO 117
   85 CP = 6.2103006
      GO TO 116
   87 CP = 0.0
      GO TO 117
  112 PRINT 4, CP , R , EPAX , CEX,IPHASE
      STOP
  116 IF(CP.EQ.0.0) GO TO 117
      IF(EPAR.EQ.090.0) FACT = ABS(FACT)
      CP = CP*FACT
  117 EPAR = EPAX
      RETURN
      END

      SUBROUTINE CARDS(W,D,AVE,NSPH,TFAC,IDENS,K,AZI)
C
      CHARACTER*10 IDENS(16)
      DIMENSION W(24),D(24),IW(25),ID(25)
      write(10,1)(IDENS(I),I=1,16)
      IAVE = (AVE+.0005) * 1000.
      IAZI = AZI
      NOS=0
      IF(K.EQ.4) NOS=1
      IF(K.EQ.3.AND.IAZI.EQ.0) NOS=1
      IF(K.EQ.3.AND.IAZI.NE.0) NOS=2
C     write(10,2)IAVE,NSPH,NOS,TFAC
      write(10,2)IAVE
C
      DO 5 I=1,24
      IW(I) = (W(I)+.05)*10.
      ID(I) = (D(I)+.0005)*1000.
    5 CONTINUE
      IW(25) = 0
      ID(25) = 0
      N=1
      write(10,6)N,ID(6),IW(6),ID(18),IW(18),ID(10),IW(10),ID(2),IW(2)
     *,ID(7),IW(7),ID(12),IW(12),ID(8),IW(8)
      N=2
      write(10,6)N,ID(25),IW(25),ID(19),IW(19),ID(25),IW(25),ID(23)
     *,IW(23),ID(20),IW(20),ID(25),IW(25),ID(11),IW(11)
      N=3
      write(10,6)N,ID(13),IW(13),ID(22),IW(22),ID(25),IW(25),ID(5),IW
     *(5),ID(1),IW(1),ID(25),IW(25),ID(25),IW(25)
      N=4
      write(10,6)N,ID(25),IW(25),ID(25),IW(25),ID(25),IW(25),ID(24)
     *,IW(24),ID(15),IW(15),ID(21),IW(21),ID(17),IW(17)
      N=5
      write(10,6)N,ID(16),IW(16),ID(14),IW(14),ID(25),IW(25),
     *ID(25),IW(25),ID(4),IW(4),ID(25),IW(25),ID(3),IW(3)
      N=6
      write(10,6)N,ID(9),IW(9),ID(25),IW(25)
C
      IEBB=IAZI+180
      if(iebb.gt.360)iebb=iebb-360
      IF(NOS.NE.1) write(10,7)IAZI,IEBB
      RETURN
    1 FORMAT(8A10)
    2 FORMAT(I6,I2,I2,2H 0,F10.8)
    6 FORMAT(7X,I1,7(I5,I4))
    7 FORMAT(2I5)
      END
      SUBROUTINE VELDIR1(CMAJ,CMIN,VEL,DIR,XMAJ)
      REAL*8 RADDEG
      DATA RADDEG /57.295779513802D0/
      IF(CMAJ.NE.0.0) GO TO 113
      IF(CMIN)110,111,112
  110 DIR =270.0
      GO TO 117
  111 DIR = 0.0
      GO TO 117
  112 DIR = 90.0
      GO TO 117
  113 DIR = ATAN(CMIN/CMAJ)
      DIR = RADDEG * DIR
      IF(CMAJ.GT.0.0) GO TO 114
      DIR = DIR + 180.0
      GO TO 117
  114 IF(CMIN.GT.0.0) GO TO 117
      DIR = DIR + 360.0
  117 VEL = SQRT(CMIN**2 + CMAJ**2)
      DIR = DIR + XMAJ
      IF(DIR.GE.360.0) DIR = DIR - 360.0
      RETURN
      END
      SUBROUTINE ELIPSE (D,SAVH,W,SAVK,AZI,STM,GONL)
C
C     DETERMINES AND PRINTS ELLIPSE PARAMETERS
      CHARACTER*10 NAMES(24),C(5,2)*6
!      DIMENSION D(25),SAVH(25),W(24),SAVK(25),C(5,2),SPCONS(24)
      DIMENSION D(25),SAVH(25),W(24),SAVK(25),SPCONS(24)
      DIMENSION A(24),AP(24),GSAVE(24),HSAVE(24),ISAVE(24)
      DATA(NAMES(I2),I2=1,24)/'J(1)      ','K(1)      ','K(2)      '
     1 ,'L(2)      ','M(1)      ','M(2)      ','M(4)      '
     2 ,'M(6)      ','M(8)      ','N(2)      ','2N(2)     '
     3 ,'O(1)      ','OO(1)     ','P(1)      ','Q(1)      '
     4 ,'2Q(1)     ','R(2)      ','S(2)      ','S(4)      '
     5 ,'S(6)      ','T(2)      ','LAMBDA(2) ','NU(2)     '
     7 ,'RHO(1)    '/
      DATA(C(I,1),I=1,5)/4H  CW,4H    ,4H    ,4H    ,4H    /
      DATA(C(I,2),I=1,5)/4H CCW,4H    ,4H    ,4H    ,4H    /
      DATA(SPCONS(I2),I2=1,24)/15.5854433,15.0410686,30.0821373,29.52847
     A89,14.4966939,28.9841042,57.9682084,86.9523127,115.9364169,28.4397
     B295,27.8953548,13.9430356,16.1391017,14.9589314,13.3986609,12.8542
     C862,30.0410667,30.0000000,60.0000000,90.0000000,29.9589333,29.4556
     D253,28.5125831,13.4715145/
       lout1=6
c      print*, '  ELLIPSE SUBROUTINE CALLED TO WRITE PARAMETERS'
      AZ = AZI*0.017453292
      I1=1
      DO 95444 I=1,24
      ISAVE(I)=0
      GSAVE(I)=0.0
      HSAVE(I)=0.0
95444 CONTINUE
      WRITE (lout1,40)
      WRITE(lout1,42)STM
      AZJ =AZI + 90.
      IF(AZJ.GE.360.)AZJ=AZJ-360.
      WRITE (lout1,43)AZI,AZJ
      WRITE (lout1,44)
      WRITE (lout1,45)
      DO 88 I=1,24
      ZETA1 = W(I)*0.017453292
      ZETA2 = SAVK(I)*0.017453292
      A1 = D(I)*COS(ZETA1)*COS(AZ)+SAVH(I)*COS(ZETA2)*SIN(AZ)
      B1 = D(I)*SIN(ZETA1)*COS(AZ)+SAVH(I)*SIN(ZETA2)*SIN(AZ)
      A2 = SAVH(I)*COS(ZETA2)*COS(AZ)-D(I)*COS(ZETA1)*SIN(AZ)
      B2 = SAVH(I)*SIN(ZETA2)*COS(AZ)-D(I)*SIN(ZETA1)*SIN(AZ)
      CC = A1*A1 + A2*A2
      DD = B1*B1 + B2*B2
      E  = A1*B1 + A2*B2
      F  = (CC-DD)/2.0
      G  = (CC+DD)/2.0
      AN = 0.0
      IF(F.NE.0.0) AN = ATAN2(E,F)
      IF(F.EQ.0.0.AND.E.GT.0.0 ) AN =  1.570796
      IF(F.EQ.0.0.AND.E.LT.0.0 ) AN = -1.570796
      IF (AN.LT.0.0)AN = 6.283185 + AN
      IF(AN.NE.0.0)GO TO 15
      H = F/COS(AN)
      GO TO 14
   15 H = E/SIN(AN)
   14 AN = AN/2.0
C
      W1 = SQRT(G+H)
      TEST=G-H
      IF(TEST.GT.0.0) GO TO 83
      W2=0.00001
      ISAVE(I1)=1
      GSAVE(I1)=G
      HSAVE(I1)=H
      GO TO 84
   83 W2 = SQRT(G-H)
   84 I1=I1+1
      WS1 = A2*COS(AN) + B2  *SIN(AN)
      WC1 = A1*COS(AN) + B1* SIN(AN)
      THETA1 = ATAN2(WS1,WC1)*57.2957812
      WS2 = B2*COS(AN) - A2 * SIN(AN)
      WC2 = B1*COS(AN) - A1  *SIN(AN)
      THETA2 = ATAN2(WS2,WC2)*57.2957812
      IF(THETA1.LT.0.0)THETA1=THETA1+360.0
      IF(THETA2.LT.0.0)THETA2=THETA2+360.0
      TEST1 = THETA1 - 89.0
      IF(TEST1.LT.0.0) TEST1 = TEST1 +360.0
      TEST1A= THETA1 - 91.0
      IF(TEST1A.LT.0.0)TEST1A= TEST1A+360.0
      IROT=2
      IF(TEST1.GE.THETA2.AND.TEST1A.LE.THETA2)IROT=1
      THETA1 = 90.0 - THETA1
      IF(THETA1.LT.0.0)THETA1 = THETA1 + 360.0
      THETA2 = 90.0 - THETA2
      IF(THETA2.LT.0.0)THETA2 = THETA2 + 360.0
      IF(IROT.EQ.2) THETA2=THETA2 + 180.
      IF(THETA2.GE.360.) THETA2 = THETA2 - 360.
      AN = AN*57.2957812
      HOUR = AN/SPCONS(I)
      IF(IROT.EQ.1)ANN = AN + 90.0
      IF(IROT.EQ.2) ANN = AN - 90.0
      IF(ANN.GT.360.) ANN = ANN - 360.0
      IF(ANN.LT.0.0)  ANN = ANN + 360.0
      HOUR1 = ANN/SPCONS(I)
      ECC = SQRT(W1**2-W2**2)/W1
      IF(IROT.EQ.1) ECC = - ECC
C
      WRITE (lout1,66)NAMES(I),THETA1,W1,AN,HOUR,THETA2,W2,ANN,HOUR1,
     1              (C(J,IROT),J=1,1),ECC,SAVH(I),SAVK(I),D(I),W(I)
      A(I) = W1
      AP(I) = AN
   88 CONTINUE
C
      WRITE (lout1,95)
      WRITE (lout1,96)
C     PRINT 98562
98562 FORMAT(1H1)
      DO 85 I2=1,24
      IF(ISAVE(I2).EQ.0) GO TO 85
      PRINT 98562
      PRINT 95652
95652 FORMAT(5X,'ATTENTION!!!,AXIS RESET')
      PRINT 95316, NAMES(I2),GSAVE(I2),HSAVE(I2)
95316 FORMAT(5X,'CONSTITUENT',2X,A10,2X,'W1 SET',5X,'G=',F12.9,5X,
     1'H=',F12.9)
      PRINT 98562
   85 CONTINUE
      RETURN
C
   40 FORMAT( ' Ellipse Parameters    (Right-Handed)')
   42 FORMAT( 80X,' Time Meridian of results =',F5.0/)
   43 FORMAT( 92X,F5.0,' Axis',4X,F5.0,' Axis')
   44 FORMAT(13X,'-----Major Ellipse Axis-----',5X,'-----Minor Ellipse A
     *xis-----',23X,'Kappa',9X,'Kappa')
   45 FORMAT(' Constituent  Dir      Ampl   Phase  Hour      Dir      Am
     *pl   Phase  Hour   Rot',4X,'Ecc',4X,'H(A)  Prime   H(A)  Prime'/)
   66 FORMAT(1X,A8 ,1X,2(3X,F5.1,2X,F7.3,3X,F5.1,1X,F5.2,2X),1A6,F8.2,2(
     *F8.3,F6.1))
   95 FORMAT( /20X,'NOTE --- Minor axis direction is considered as the m
     1ajor axis plus 90 degrees in all cases.  If'/29X,'rotation is cloc
     2kwise, phase will be major axis phase plus 90 degrees.  If counter
     3-'/29X,'clockwise, minus 90 degrees.')
   96 FORMAT(  20X,'NOTE --- Direction of major axis is usually flood di
     *rection but may be ebb sometimes.  If this'/29X,'occurs, to get th
     *e results for flood direction, add 180 degrees to all directions a
     *nd'/29X,'phases, and (0.5 X constituent period) to all hour values
     *.'/20X,'NOTE --- Ecc = Eccentricity of ellipse ( + if CCW rotation
     * )')
      END
       subroutine ntran(iunit,ifunc,ipar1,ipar2,ipar3,ipar4)
c
        dimension ipar2(700)
        ipar3=0
        go to (1,2,99,99,99,99,7,8,9,10,99,99,99,99,99,99,99,99,99,99,
     x  99,99,99),ifunc
c
        return
1       write(iunit,err=100)(ipar2(i),i=1,ipar1)
        return
c
c****       read function
c
2       read(iunit,end=101,err=100)(ipar2(i),i=1,ipar1)
        return
c
c       move block/record routine
c
7       isave=iabs(ipar1)
        if(ipar1.lt.0)then
        do 70 i=1,isave
        backspace iunit
70      continue
c
        elseif(ipar1.gt.0)then
        do 71 i=1,isave
        read(iunit,err=100,end=101)idum
71      continue
c
        endif
        return
c
c       routine to move files
c
8       isave=iabs(ipar1)
        if(ipar1.lt.0)then
        do 80 i=1,isave
        go to 100
80      continue
c
        elseif(ipar1.gt.0)then
        do 81 i=1,isave
82       read(iunit,end=81,err=100)idum
        go to 82
81      continue
c
        endif
        return
c
c       routine to write an eof on unit
c
9       endfile iunit
        return
c
c       routine to rewind the unit
c
10        rewind iunit
        return
c
c
c       catch all else
c
99      return
100      ipar3=-3
      print*,' !!!!!!!!!!! read err= ',ipar3
         return
101      ipar3=-2
         return
c
        end
C**********************************************
