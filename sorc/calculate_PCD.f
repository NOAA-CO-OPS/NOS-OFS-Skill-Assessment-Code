C lf95 
C This subroutine is used for principle component computation
C  CALL statement is as,
C        CALL prcmp(NMAX,U,V,pcd,RUV,RATIO)
C   NMAX is number of data
C   U is true eastern velocity component, input variable
C   V is true northern velocity component,input variable
C   PCD is principal current direction in degrees clockwise from north
C   RUV is Correlation coefficient of U and V 
C   RATIO is the ratio of minor axis variance to major axis variance
C  AUTHOR:  Aijun Zhang
C   Coast Survey Development Laboratory, NOS/NOAA
C  Language:  Fortran
C  
      PARAMETER(NT=100000)
      character*200 fname,FCTL,sigmaunits*20,sigmaname*40
      character*200 staid,fileshort,longnames,OCEAN_MODEL,ANAME
      character*200 filestation,outputf0,AQ,BUFFER 
      DIMENSION TIME(NT),VIN(NT),DIN(NT),SIYR(NT),SMON(NT)
     1          ,SDD(NT),SHH(NT),SMIN(NT),ifrst(NT),ilast(NT)
     2          ,XTMP(NT),YTMP(NT),U(NT),V(NT)
      FCTL='/home/net/azhang/SKILL_OFS/control_files/'
      FCTL=trim(FCTL)//'creofs_cu_station.ctl'
      OPEN(9,file=FCTL)
98    read(9,*,err=999,end=999)staid,fileshort,longnames
      read(9,*)alat,alon,dirflood,sdepth  !! read standard depth
!      write(*,*)'sdepth= ',sdepth

      OPEN(2,file=trim(fileshort)//'_nowcast.dat')
      I=1
5     READ(2,*,END=10)tday,SIYR(I),SMON(I),SDD(I),SHH(I),SMIN(I)
     1                   ,VIN(I),DIN(I)       
      I=I+1
      GOTO 5
10    CONTINUE
      close(2)
      NMAX=I-1
      DO I=1,NMAX
!          RDATA(I,1)=VIN(IBEGIN+I-1)
 !         RDATA(I,2)=DIN(IBEGIN+I-1)
          XTMP(I)=VIN(I)*sin(DIN(I)*3.14159265/180.)
          YTMP(I)=VIN(I)*cos(DIN(I)*3.14159265/180.)
      ENDDO 
C    calculate principle current direction
      CALL prcmp(NMAX,XTMP,YTMP,pcd,rxy,RATIO)
      print *,'PCD= ',trim(fileshort),PCD
      goto 98
999      STOP
      END
      SUBROUTINE prcmp(ndata,x,y,pcd,rxy,RATIO)
C This subroutine is used for principle component computation
      
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
        WRITE(*,*)'RATIO=',RATIO
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
C*****************************************************************************
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
