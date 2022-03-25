C     PROGRAM read_netcdf_fcst.f
C  Compile on opmarine
C   f90 read_netcdf_fcst.f -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include \
C    -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf -o read_netcdf_fcst.x
C  Compile on Linux 
C   lf95 read_netcdf_fcst.f -I/usr/Local/Linux/netcdf/include \
C    -L/usr/Local/Linux/netcdf/lib -lnetcdf -o ../binLinux/read_netcdf_fcst.x
C
C  Modified by Lianyuan Zheng on 03/10/2017

      PROGRAM READ_NETCDF_FCST
      include 'netcdf.inc'
      CHARACTER*50  SIGMANAME
      CHARACTER*200 FNAME,FCTL,STAID,LONGNAMES,CC,OCEAN_MODEL,ANAME
      CHARACTER*200 FILESHORT(1000)
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE
      REAL*8 JULIAN,YEARB,MONTHB,DAYB,HOURB
      INTEGER, ALLOCATABLE ::  ISTA(:)
      REAL*4, ALLOCATABLE ::  SDEPTH(:)
      REAL*4, ALLOCATABLE ::  TIME(:,:,:)
      REAL*4, ALLOCATABLE ::  ALON(:)
      REAL*4, ALLOCATABLE ::  ALAT(:)
      REAL*4, ALLOCATABLE ::  LON(:)
      REAL*4, ALLOCATABLE ::  LAT(:)
      REAL*4, ALLOCATABLE ::  DEPTH(:)
      REAL*4, ALLOCATABLE ::  SIGMA(:)
      REAL*4, ALLOCATABLE ::  SIGLAY(:,:)
      REAL*4, ALLOCATABLE ::  ZVAL(:,:,:)
      REAL*4, ALLOCATABLE ::  ZVAL0(:,:,:)
      REAL*4, ALLOCATABLE ::  ZSIGMA(:)
      REAL*4, ALLOCATABLE ::  UK(:)
      REAL*4, ALLOCATABLE ::  VK(:)

      REAL*4, ALLOCATABLE ::  TIMETMP(:)
      REAL*4, ALLOCATABLE ::  ZETA(:,:)
      REAL*4, ALLOCATABLE ::  UEL(:,:,:)
      REAL*4, ALLOCATABLE ::  VEL(:,:,:)
      REAL*4, ALLOCATABLE ::  TEMP(:,:,:)
      REAL*4, ALLOCATABLE ::  SALT(:,:,:)
      INTEGER BASE_DATE(4)
      INTEGER TZ_MODEL

      READ(5,*) IYRS, IMMS, IDDS, IHHS, IMNS
      READ(5,*) DELT0, NCYCLE, KINDAT, NFILES, NFDURATION
      DELT_MIN = DELT0
      DELT0 = DELT0/60.0  !! 11/30/2006 convert DELT0 from minutes into hours
      READ(5,'(A200)') FCTL
      READ(5,*) TZ_MODEL 
      READ(5,'(A200)') OCEAN_MODEL 

      OPEN(9,FILE = TRIM(FCTL))
      DO I = 1,10000
        READ(9,*, END = 3)
        READ(9,*, END = 3)
      ENDDO

3     NSTATION = I - 1
      IF(NSTATION .EQ. 0) THEN
        WRITE(*,*) 'Error in reading the sattion control file, STOP!'
        STOP
      END IF
      REWIND(9)

      ALLOCATE(ALON(NSTATION))
      ALLOCATE(ALAT(NSTATION))
      ALLOCATE(SDEPTH(NSTATION))
      ALLOCATE(ISTA(NSTATION))

      YEARB  = IYRS*1.0
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JDAY0  = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      DO I = 1, NSTATION
        READ(9,*,ERR = 9,END = 9) STAID, FILESHORT(I), LONGNAMES
        READ(9,*) ALAT(I), ALON(I), DIRFLOOD, SDEPTH(I), TDEPTH
      ENDDO
9     CLOSE(9)

      IUNIT = 50
      DO I = 1, NSTATION
        OPEN(IUNIT+I,FILE = TRIM(FILESHORT(I)))
      ENDDO  

      NCUT = INT(NFDURATION/DELT0+0.1) + 1
      NDAYS = INT(NFILES/NCYCLE)
      DO I = 1, NDAYS
        DO J = 1, NCYCLE
          READ(5,'(A200)') FNAME
          READ(5,*) IYR1, IMM1, IDD1, IHH1, IMIN1
          WRITE(*,*) ' FNAME=',TRIM(FNAME)
          CALL FIRST_READ_NETCDF(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,
     1      BASE_DATE,SIGMANAME,THETA_S,THETA_B,TCLINE,HC,OCEAN_MODEL)
          IF((I .EQ. 1) .AND. (J .EQ. 1)) THEN
            IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') THEN
              ALLOCATE(SIGLAY(NSTA,KB))
              STATUS = NF_OPEN(fname, NF_NOWRITE, NCID)
              STATUS = NF_INQ_VARID(NCID,'siglay',IDVAR)
              IF(STATUS .EQ. NF_NOERR) THEN
                STATUS = NF_GET_VAR_REAL(NCID,IDVAR,SIGLAY)
              ENDIF
              STATUS=NF_CLOSE(NCID)
            ENDIF	     	
	  ENDIF

          ALLOCATE(TIMETMP(NT))
          ALLOCATE(ZETA(NSTA,NT))
          ALLOCATE(UEL(NSTA,KB,NT))
          ALLOCATE(VEL(NSTA,KB,NT))
          ALLOCATE(TEMP(NSTA,KB,NT))
          ALLOCATE(SALT(NSTA,KB,NT))
          ALLOCATE(DEPTH(NSTA))
          ALLOCATE(SIGMA(KB))
          ALLOCATE(ZSIGMA(KB))
          ALLOCATE(UK(KB))
          ALLOCATE(VK(KB))
          ALLOCATE(LON(NSTA))
          ALLOCATE(LAT(NSTA))

          YEARB  = BASE_DATE(1)
          MONTHB = BASE_DATE(2)
          DAYB   = BASE_DATE(3)
          HOURB  = BASE_DATE(4)
          JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)

          YEARB  = IYR1
          MONTHB = IMM1
          DAYB   = IDD1
          HOURB  = IHH1 + IMIN1/60.0
          JDAY1  = JULIAN(YEARB,MONTHB,DAYB,HOURB)
          TSTART0 = JDAY1 - JBASE_DATE
          TEND    = TSTART0 + NFDURATION/24.0 !!! Forecast duration
          TSTART1 = JDAY1 - JDAY0 + 1

          IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
            IF(ALLOCATED(ZVAL)) DEALLOCATE(ZVAL)
            IF(ALLOCATED(ZVAL0)) DEALLOCATE(ZVAL0)
            ALLOCATE(ZVAL(NSTA,KB,NT))
            ALLOCATE(ZVAL0(NSTA,KB,NT))
            STATUS = NF_OPEN(FNAME, NF_NOWRITE, NCID)
            STATUS = NF_INQ_VARID(NCID,'ZVAL',IDVAR)
            IF(STATUS .EQ. NF_NOERR) THEN
              STATUS = NF_GET_VAR_REAL(NCID,IDVAR,ZVAL0)
            ENDIF
            STATUS = NF_CLOSE(NCID)

C  Reverse vertical coordinates (K=1 for bottom and K=KB for surface)
C  Change negative sign to positive sign
C  zeta=-9999.0 is dry cell for SELFE
            DO NS = 1, NSTA   
              DO N = 1, NT
                DO K = 1, KB
	          K1 = KB - K + 1
                  ZVAL(NS,K,N) = -ZVAL0(NS,K1,N)
	        ENDDO
	      ENDDO
	    ENDDO  
          ENDIF

   	  CALL READ_NETCDF(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,
     1      TIMETMP,ZETA,UEL,VEL,TEMP,SALT,DEPTH,SIGMA,
     2      LON,LAT,KINDAT,OCEAN_MODEL)

          DO NSS0 = 1, NSTATION
            DISTMIN = 99999999.0
            NODE_FIND = 0
            IF(ALON(NSS0) .GT. 180.0) ALON(NSS0) = ALON(NSS0)-360.0
            DO NSS = 1, NSTA
              IF(LON(NSS) .GT. 180.0) LON(NSS) = LON(NSS)-360.0	     
              IF(ABS(LAT(NSS)) .LE. 90.0) THEN
                CALL DIST(LAT(NSS),LON(NSS),ALAT(NSS0),ALON(NSS0),DIS)
                IF(DIS .LE. DISTMIN) THEN
   	 	  DISTMIN = DIS
   		  NODE_FIND = NSS
                ENDIF
              ENDIF
            ENDDO
            ISTA(NSS0) = NODE_FIND
          ENDDO  

          N = 0
          DO I1 = 1, NT
            IF((TIMETMP(I1) .GE. TSTART0-0.0001) .AND.
     1         (TIMETMP(I1) .LE. TEND+0.0001)) THEN
              N = N + 1
              JDAY = TIMETMP(I1) + JBASE_DATE
              CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
              IYEAR = INT(YEARB+0.001)
              ICM   = INT(MONTHB+0.001)
              ICD   = INT(DAYB+0.001)
              IHR   = INT(HOURB+0.001)
              IMN0  = INT((HOURB-IHR)*60+0.1)
	      IMN   = INT(IMN0/DELT_MIN+0.5)*DELT_MIN

              YEARB  = IYEAR
              MONTHB = ICM
              DAYB   = ICD
              HOURB  = IHR + IMN/60.0
              JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB)
              CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
              IYEAR = INT(YEARB+0.001)
              ICM   = INT(MONTHB+0.001)
              ICD   = INT(DAYB+0.001)
              IHR   = INT(HOURB+0.001)
              IMN0  = INT((HOURB-IHR)*60+0.1)
	      IMN   = INT(IMN0/DELT_MIN+0.5)*DELT_MIN
              DAYJ  = JDAY - JDAY0 + 1.0
	      
              DO NS = 1, NSTATION
                IF(TRIM(OCEAN_MODEL) .EQ. 'ADCIRC' )THEN
                  WRITE(IUNIT+NS,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,
     1              ZETA(ISTA(NS),I1)
                  CYCLE
                ENDIF

                ELE  = ZETA(ISTA(NS),I1)
                DTMP = DEPTH(ISTA(NS))
                IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') THEN
                  CALL SIGMA2Z_ROMS_FIX(SIGMA,DTMP,ELE,KB,ZSIGMA,
     1                 HC,THETA_S,THETA_B,TCLINE)
                ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'POM') THEN
                  DO K = 1, KB
                    IF(SIGMA(K) .GT. 0.0) SIGMA(K) = -SIGMA(K)
                  ENDDO
                  CALL SIGMA2Z_POM(SIGMA,DTMP,ELE,KB,ZSIGMA)
                ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') THEN
                  DO K = 1, KB
	            SIGMA(K) = SIGLAY(ISTA(NS),K)
                    IF(SIGMA(K) .GT. 0.0) SIGMA(K) = -SIGMA(K)
                  ENDDO
                  CALL SIGMA2Z_POM(SIGMA,DTMP,ELE,KB,ZSIGMA)  
                ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'SELFE' )THEN
                  DO K = 1, KB
	            ZSIGMA(K) = ZVAL(ISTA(NS),K,I1)
                  ENDDO
                ELSE
                  DO K = 1, KB
                    ZSIGMA(K) = SIGMA(K)
                  ENDDO
                ENDIF

                IF(KINDAT .EQ. 1) THEN
                  IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
                    DO KKK = 1, KB
                      UK(KKK) = ZSIGMA(KB-KKK+1)
                    ENDDO
                    DO KKK = 1, KB
                      ZSIGMA(KKK) = UK(KKK)
                    ENDDO
                    DO KKK = 1, KB
                      UK(KKK) = UEL(ISTA(NS),KB-KKK+1,I1)
                      VK(KKK) = VEL(ISTA(NS),KB-KKK+1,I1)
                    ENDDO
                  ELSE     
                    DO KKK = 1, KB
                     UK(KKK) = UEL(ISTA(NS),KKK,I1)
                     VK(KKK) = VEL(ISTA(NS),KKK,I1)
                    ENDDO
                  ENDIF

                  IF(SDEPTH(NS) .LT. ZSIGMA(1)) THEN
                     UTMP = UK(1)
                     VTMP = VK(1)
                  ELSEIF(SDEPTH(NS) .GT. ZSIGMA(KB)) THEN
                     UTMP = UK(KB)
                     VTMP = VK(KB)
                  ELSE
                    DO K9 = 1, KB-1
                      IF((SDEPTH(NS) .GE. ZSIGMA(K9)) .AND.
     1                   (SDEPTH(NS) .LT. ZSIGMA(K9+1))) THEN
                        X1 = ZSIGMA(K9)
                        X2 = ZSIGMA(K9+1)
                        Y1 = UK(K9)
                        Y2 = UK(K9+1)
                        CALL LINEAR(X1,Y1,X2,Y2,SDEPTH(NS),UTMP)
                        Y1 = VK(K9)
                        Y2 = VK(K9+1)
                        CALL LINEAR(X1,Y1,X2,Y2,SDEPTH(NS),VTMP)
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF

                  CALL VELDIR(VTMP,UTMP,AANGLE,AVEL)
                  IF(AVEL .LE. 0.0) AANGLE = 0.0
                  IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0

                  IF((I1 .EQ. 1) .AND. 
     1               (TIMETMP(1) .GT. TSTART0+0.0001) ) THEN
                    WRITE(IUNIT+NS,100)TSTART1,IYR1,IMM1,IDD1,IHH1,0,
     1                AVEL,AANGLE,UTMP,VTMP 
                  ENDIF
                  WRITE(IUNIT+NS,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,
     1                AVEL,AANGLE,UTMP,VTMP

                ELSE IF(KINDAT .eq. 2) THEN
                  IF((I1 .EQ. 1) .AND. 
     1               (TIMETMP(1) .GT. TSTART0+0.0001) ) THEN
                    WRITE(IUNIT+NS,100)TSTART1,IYR1,IMM1,IDD1,IHH1,0,
     1                ZETA(ISTA(NS),I1) 
                  ENDIF
                  WRITE(IUNIT+NS,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,
     1                 ZETA(ISTA(NS),I1) 

                ELSE IF(KINDAT .eq. 3) THEN
                  IF(zSIGMA(1) .GT. ZSIGMA(KB)) THEN
                    DO KKK = 1, KB
                      UK(KKK) = ZSIGMA(KB-KKK+1)
                    ENDDO
                    DO KKK = 1, KB
                      ZSIGMA(KKK) = UK(KKK)
                    ENDDO
                    DO KKK = 1, KB
                      UK(KKK) = TEMP(ISTA(NS),KB-KKK+1,I1)
                    ENDDO
                  ELSE     
                    DO KKK = 1, KB
                     UK(KKK) = TEMP(ISTA(NS),KKK,I1)
                    ENDDO
                  ENDIF

                  IF(SDEPTH(NS) .LT. ZSIGMA(1)) THEN
                    UTMP = UK(1)
                  ELSEIF(SDEPTH(NS) .GT. ZSIGMA(KB)) THEN
                    UTMP = UK(KB)
                  ELSE
                    DO K9 = 1, KB-1
                      IF((SDEPTH(NS) .GE. ZSIGMA(K9)) .AND.
     1                   (SDEPTH(NS) .LT. ZSIGMA(K9+1))) THEN
                        X1 = ZSIGMA(K9)
                        X2 = ZSIGMA(K9+1)
                        Y1 = UK(K9)
                        Y2 = UK(K9+1)
                        CALL LINEAR(X1,Y1,X2,Y2,SDEPTH(NS),UTMP)
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF

                  IF((I1 .EQ. 1) .AND. 
     1               (TIMETMP(1) .GT. TSTART0+0.0001) ) THEN
                    WRITE(IUNIT+NS,100)TSTART1,IYR1,IMM1,IDD1,IHH1,
     1                0,UTMP
                  ENDIF
                  WRITE(IUNIT+NS,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP

                ELSE IF(KINDAT .EQ. 4) THEN
                  IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
                    DO KKK = 1, KB
                      UK(KKK) = ZSIGMA(KB-KKK+1)
                    ENDDO
                    DO KKK = 1, KB
                      ZSIGMA(KKK) = UK(KKK)
                    ENDDO
                    DO KKK = 1, KB
                      UK(KKK) = SALT(ISTA(NS),KB-KKK+1,I1)
                    ENDDO
                  ELSE     
                    DO KKK = 1, KB
                      UK(KKK) = SALT(ISTA(NS),KKK,I1)
                    ENDDO
                  ENDIF

                  IF(SDEPTH(NS) .LT. ZSIGMA(1)) THEN
                    UTMP = UK(1)
                  ELSEIF(SDEPTH(NS) .GT. ZSIGMA(KB)) THEN
                    UTMP = UK(KB)
                  ELSE
                    DO K9 = 1 ,KB-1
                      IF((SDEPTH(NS) .GE. ZSIGMA(K9)) .AND.
     1                   (SDEPTH(NS) .LT. ZSIGMA(K9+1))) THEN
                        X1 = ZSIGMA(K9)
                        X2 = ZSIGMA(K9+1)
                        Y1 = UK(K9)
                        Y2 = UK(K9+1)
                        CALL LINEAR(X1,Y1,X2,Y2,SDEPTH(NS),UTMP)
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF

                  IF((I1 .EQ. 1) .AND. 
     1               (TIMETMP(1) .GT. TSTART0+0.0001) ) THEN
                    WRITE(IUNIT+NS,100)TSTART1,IYR1,IMM1,IDD1,IHH1,
     1                0,UTMP
                  ENDIF
                  WRITE(IUNIT+NS,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP
                ENDIF
              ENDDO
            ENDIF  
          ENDDO

          DEALLOCATE(TIMETMP)
          DEALLOCATE(ZETA)
          DEALLOCATE(UEL)
          DEALLOCATE(VEL)
          DEALLOCATE(TEMP)
          DEALLOCATE(SALT)
          DEALLOCATE(DEPTH)
          DEALLOCATE(SIGMA)
          DEALLOCATE(ZSIGMA)
          DEALLOCATE(UK)
          DEALLOCATE(VK)
          DEALLOCATE(LON)
          DEALLOCATE(LAT)
        ENDDO
      ENDDO
      DO I = 1, NSTATION
        CLOSE(IUNIT+I)
      ENDDO  
100   FORMAT(F10.5,I5,4I3,50F10.4)

      STOP
      END
