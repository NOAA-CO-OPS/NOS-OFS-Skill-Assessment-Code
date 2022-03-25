C     PROGRAM read_netcdf_nowcast.f
C f90 read_netcdf_now.f -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include \
C -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf -o read_netcdf_now.x
C     ifc read_netcdf.f  -I/disks/NASWORK/ngofs/oqctools/netcdflinux/include \
C       -L/disks/NASWORK/ngofs/oqctools/netcdflinux/lib -lnetcdf -o read_netcdf.x
C Compile:
C     lf95 read_netcdf_nowcast.f read_netcdf_subs.f utility.f -I/usr/local/include \
C       -L/usr/local/lib -lnetcdf -o read_netcdf_nowcast.x
C
C  Modified by Lianyuan Zheng on 03/10/2017

      PROGRAM READ_NETCDF_NOWCAST
      PARAMETER (NMAX=500000)
      include 'netcdf.inc'

      character*50  sigmaname
      character*100 staid(1000),fileshort(1000),longnames(1000)
      character*200 FNAME,FCTL,BUFFER,OCEAN_MODEL,ANAME
      real*4, allocatable ::  time(:)
      real*4, allocatable ::  lon(:)
      real*4, allocatable ::  lat(:)
      real*4, allocatable ::  depth(:)
      real*4, allocatable ::  sigma(:)
      real*4, allocatable ::  siglay(:,:)
      real*4, allocatable ::  zval(:,:,:)
      real*4, allocatable ::  zval0(:,:,:)
      real*4, allocatable ::  zsigma(:)
      real*4, allocatable ::  h(:,:)
      real*4, allocatable ::  u(:,:,:)
      real*4, allocatable ::  v(:,:,:)
      real*4, allocatable ::  t(:,:,:)
      real*4, allocatable ::  s(:,:,:)
      real*4, allocatable ::  lon_sta(:)
      real*4, allocatable ::  lat_sta(:)
      real*4, allocatable ::  dirflood_sta(:)
      real*4, allocatable ::  sdepth_sta(:)

      real*4, allocatable ::  timetmp(:)
      real*4, allocatable ::  lontmp(:)
      real*4, allocatable ::  lattmp(:)
      real*4, allocatable ::  depthtmp(:)
      real*4, allocatable ::  zeta(:,:)
      real*4, allocatable ::  uel(:,:,:)
      real*4, allocatable ::  vel(:,:,:)
      real*4, allocatable ::  temp(:,:,:)
      real*4, allocatable ::  salt(:,:,:)
      real*4, allocatable ::  UK(:)
      real*4, allocatable ::  VK(:)
      integer, allocatable :: INDEX_STA(:)
      real*4, allocatable ::  siglaytmp(:,:)
      real*8 JDAY, JDAY0, JDAY1, JBASE_DATE, JULIAN
      real*8 YEARB, MONTHB, DAYB, HOURB
      integer base_date(4)
      integer TZ_MODEL

      READ(5,*) IYRS, IMMS, IDDS, IHHS, IMNS
      READ(5,*) DELT, NCYCLE, KINDAT, NFILES
      READ(5,'(A200)') FCTL
      READ(5,*) TZ_MODEL 
      READ(5,'(A200)') OCEAN_MODEL 

      DELT_MIN = DELT
      DELT = DELT/60.0
      NCUT = INT(24/NCYCLE/DELT+0.1)

      YEARB  = IYRS*1.0
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JDAY0 = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      WRITE(*,*)'DELT=',DELT,'NCYCLE=',NCYCLE,'NMAX=',NMAX

      READ(5,'(a200)') FNAME
C      WRITE(*,*) ' FNAME = ', TRIM(FNAME)
      CALL FIRST_READ_NETCDF(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,
     1 BASE_DATE,SIGMANAME,THETA_S,THETA_B,TCLINE,HC,OCEAN_MODEL)

      I = 0
      OPEN(9,FILE = FCTL)
C10    READ(9,'(A200)',ERR=20,END=30) BUFFER
C      BUFFER = TRIM(ADJUSTL(BUFFER))
C      LEN1   = LEN_TRIM(BUFFER)
C      IF(LEN1 .LE. 0) GOTO 10
C      I = I + 1
C      GOTO 10
C20    CONTINUE 
C      WRITE(*,*) 'THERE IS ERROR IN READING SATTION CONTROL FILE'
C      STOP
C30    CONTINUE
C      NSTA_CTL = I/2   ! Each station hase 2 lines
      DO I = 1,10000
        READ(9,*, END = 3)
        READ(9,*, END = 3)
      ENDDO

3     NSTA_CTL = I - 1
      IF(NSTA_CTL .eq. 0) THEN
        WRITE(*,*) 'THERE IS ERROR IN READING SATTION CONTROL FILE'
        STOP
      END IF
      REWIND(9)

      I = 0
      ALLOCATE(LON_STA(NSTA_CTL) )
      ALLOCATE(LAT_STA(NSTA_CTL) )
      ALLOCATE(DIRFLOOD_STA(NSTA_CTL) )
      ALLOCATE(SDEPTH_STA(NSTA_CTL) )
      ALLOCATE(INDEX_STA(NSTA_CTL) )
C40    READ(9,'(A200)',ERR=42,END=46) BUFFER
C      BUFFER = TRIM(ADJUSTL(BUFFER))
C      LEN1   = LEN_TRIM(BUFFER)
C      IF(LEN1 .LE. 0) GOTO 40
C      I = I + 1
C      READ(BUFFER,*) STAID(I),FILESHORT(I),LONGNAMES(I)
C      READ(9,*) LAT_STA(I),LON_STA(I),DIRFLOOD_STA(I),SDEPTH_STA(I)
C      GOTO 40
C42    CONTINUE 
C      WRITE(*,*) 'THERE IS ERROR IN READING SATTION CONTROL FILE'
C      STOP
C46    CONTINUE
C      NSTA_CTL = I

      DO I = 1, NSTA_CTL
        READ(9,*,ERR=9,END=9) STAID(I),FILESHORT(I),LONGNAMES(I)
        READ(9,*) LAT_STA(I),LON_STA(I),DIRFLOOD_STA(I),SDEPTH_STA(I)
      ENDDO
9     CLOSE(9)

      ALLOCATE(TIME(NMAX))
      ALLOCATE(H(NSTA_CTL,NMAX))
      IF(KINDAT .EQ. 1) THEN
        ALLOCATE(U(NSTA_CTL,KB,NMAX))
        ALLOCATE(V(NSTA_CTL,KB,NMAX))
      ELSEIF(KINDAT .EQ. 3) THEN
        ALLOCATE(T(NSTA_CTL,KB,NMAX))
      ELSEIF(KINDAT .EQ. 4) THEN
        ALLOCATE(S(NSTA_CTL,KB,NMAX))
      ENDIF

      ALLOCATE(DEPTH(NSTA_CTL))
      ALLOCATE(SIGMA(KB))
      ALLOCATE(ZSIGMA(KB))
      ALLOCATE(UK(KB))
      ALLOCATE(VK(KB))
      ALLOCATE(LON(NSTA_CTL))
      ALLOCATE(LAT(NSTA_CTL))
      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
        ALLOCATE(ZVAL(NSTA_CTL,KB,NMAX))
      ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') THEN
        IF(ALLOCATED(SIGLAY)) DEALLOCATE(SIGLAY)
        ALLOCATE(SIGLAY(NSTA_CTL,KB))
      ENDIF

      REWIND(5)
      READ(5,*)
      READ(5,*)
      READ(5,*)
      READ(5,*)
      READ(5,*)
      N = 0
      TIME0 = 0.0
      DO I = 1, NFILES
        READ(5,'(A200)') FNAME
        WRITE(*,*) ' FNAME=',TRIM(FNAME)
        CALL FIRST_READ_NETCDF(NSTA1,NT1,KB1,NCHAR1,MESHDIM1,FNAME,
     1    BASE_DATE,SIGMANAME,THETA_S,THETA_B,TCLINE,HC,OCEAN_MODEL)

        ALLOCATE(TIMETMP(NT1))
        ALLOCATE(LONTMP(NSTA1))
        ALLOCATE(LATTMP(NSTA1))
        ALLOCATE(DEPTHTMP(NSTA1))
        ALLOCATE(ZETA(NSTA1,NT1))
        ALLOCATE(UEL(NSTA1,KB1,NT1))
        ALLOCATE(VEL(NSTA1,KB1,NT1))
        ALLOCATE(TEMP(NSTA1,KB1,NT1))
        ALLOCATE(SALT(NSTA1,KB1,NT1))
        ALLOCATE(SIGLAYTMP(NSTA1,KB1))
        YEARB  = BASE_DATE(1)
        MONTHB = BASE_DATE(2)
        DAYB   = BASE_DATE(3)
        HOURB  = BASE_DATE(4)
        JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)

        IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
          IF (ALLOCATED(ZVAL0)) DEALLOCATE(ZVAL0)
          ALLOCATE(ZVAL0(NSTA1,KB1,NT1))
          STATUS = NF_OPEN(TRIM(fname), NF_NOWRITE, NCID)
          STATUS = NF_INQ_VARID(NCID,'zval',IDVAR)
          IF(STATUS .EQ. NF_NOERR) THEN
            STATUS = NF_GET_VAR_REAL(NCID,IDVAR,zval0)
          ENDIF
          STATUS = NF_CLOSE(NCID)
        ENDIF

	CALL READ_NETCDF(NSTA1,NT1,KB1,NCHAR1,MESHDIM1,FNAME,
     1    TIMETMP,ZETA,UEL,VEL,TEMP,SALT,DEPTHTMP,SIGMA,
     2    LONTMP,LATTMP,KINDAT,OCEAN_MODEL)

        DO NSS = 1, NSTA_CTL
          DISTMIN = 99999999.0
          ALON = LON_STA(NSS)
          ALAT = LAT_STA(NSS)
          IF(ALON .GT. 180.0) ALON = ALON - 360.0
          DO ITMP = 1, NSTA1
            IF(LONTMP(ITMP) .GT. 180.0) LONTMP(ITMP)=LONTMP(ITMP)-360.0
            IF(ABS(LATTMP(ITMP)) .LE. 90.0) THEN
              CALL DIST(LATTMP(ITMP),LONTMP(ITMP),ALAT,ALON,DIS)    
              IF(DIS .LE. DISTMIN) THEN
                DISTMIN = DIS
                INDEX_STA(NSS) = ITMP
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        IF(I .EQ. 1) THEN
          DO I2 = 1, NSTA_CTL
            DEPTH(I2) = DEPTHTMP(INDEX_STA(I2))
            LON(I2) = LONTMP(INDEX_STA(I2))
            LAT(I2) = LATTMP(INDEX_STA(I2))
          ENDDO
          IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM' ) THEN
            IF(ALLOCATED(siglaytmp) ) DEALLOCATE(siglaytmp)
            ALLOCATE(siglaytmp(NSTA1,KB1))
            STATUS = NF_OPEN(TRIM(fname), NF_NOWRITE, NCID)
            STATUS = NF_INQ_VARID(NCID,'siglay',IDVAR)
            IF(STATUS .EQ. NF_NOERR) THEN
              STATUS = NF_GET_VAR_REAL(NCID,IDVAR,siglaytmp)
            ENDIF
            STATUS = NF_CLOSE(NCID)
            DO I2 = 1, NSTA_CTL
              DO I3 = 1,KB
                SIGLAY(I2,I3) = SIGLAYTMP(INDEX_STA(I2),I3)
              ENDDO
            ENDDO 
          ENDIF
	ENDIF   

        DO I1 = 1, NT1
          TIME9 = TIMETMP(I1) + JBASE_DATE - JDAY0
C  Zheng to include the first data in the first date
          IF((TIME9 .GE. TIME0 .AND. N .EQ. 0) .OR. 
     1       (TIME9 .GT. TIME0 .AND. N .GE. 1)) THEN
            IF(N .GT. NMAX) THEN
              WRITE(*,*)'Actual data points exceeds array upper bound'
              WRITE(*,*)'Need to reset NMAX OR cut off the data points'
              WRITE(*,*)'NMAX= ', NMAX, '  N= ', N  
              GOTO 97
            ENDIF 
            N = N + 1
            TIME(N) = TIME9
	    TIME0 = TIME9
            DO I2 = 1, NSTA_CTL
              H(I2,N) = ZETA(INDEX_STA(I2),I1)
              DO I3 = 1, KB1
                IF(KINDAT .EQ. 1) THEN
                  U(I2,I3,N) = UEL(INDEX_STA(I2),I3,I1)
                  V(I2,I3,N) = VEL(INDEX_STA(I2),I3,I1)
                ELSEIF (KINDAT .EQ. 3) THEN
                  T(I2,I3,N) = TEMP(INDEX_STA(I2),I3,I1)
                ELSEIF (KINDAT .EQ. 4) THEN
                  S(I2,I3,N) = SALT(INDEX_STA(I2),I3,I1)
                ENDIF
C  Reverse vertical coordinates (K=1 for bottom and K=KB for surface)
C  zeta=-9999.0 is dry cell for SELFE
                IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
  	          K1 = KB1 - I3 + 1
                  ZVAL(I2,I3,N) = -ZVAL0(INDEX_STA(I2),K1,I1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF(ALLOCATED(TIMETMP)) DEALLOCATE(TIMETMP)
        IF(ALLOCATED(LONTMP)) DEALLOCATE(LONTMP)
        IF(ALLOCATED(LATTMP)) DEALLOCATE(LATTMP)
        IF(ALLOCATED(DEPTHTMP)) DEALLOCATE(DEPTHTMP)
        IF(ALLOCATED(ZETA)) DEALLOCATE(ZETA)
        IF(ALLOCATED(UEL)) DEALLOCATE(UEL)
        IF(ALLOCATED(VEL)) DEALLOCATE(VEL)
        IF(ALLOCATED(TEMP)) DEALLOCATE(TEMP)
        IF(ALLOCATED(SALT)) DEALLOCATE(SALT)
        IF(ALLOCATED(SIGLAYTMP)) DEALLOCATE(SIGLAYTMP)
      ENDDO
97    WRITE(*,*) 'The total data points is ', N

      IF(ALLOCATED(TIMETMP)) DEALLOCATE(TIMETMP)
      IF(ALLOCATED(ZETA)) DEALLOCATE(ZETA)
      IF(ALLOCATED(UEL))  DEALLOCATE(UEL)
      IF(ALLOCATED(VEL))  DEALLOCATE(VEL)
      IF(ALLOCATED(TEMP))  DEALLOCATE(TEMP)
      IF(ALLOCATED(SALT))  DEALLOCATE(SALT)
      IF(ALLOCATED(ZVAL0))  DEALLOCATE(ZVAL0)
      IF(ALLOCATED(LONTMP)) DEALLOCATE(LONTMP)
      IF(ALLOCATED(LATTMP)) DEALLOCATE(LATTMP)
      IF(ALLOCATED(DEPTHTMP)) DEALLOCATE(DEPTHTMP)
      IF(ALLOCATED(SIGLAYTMP)) DEALLOCATE(SIGLAYTMP)

      NT = N
C  Switch time zone to GMT/UTC if it is not GMT 
      DO N = 1, NT
        TIME(N) = TIME(N) + float(TZ_MODEL)/15.0/24.0  ! change time zone from local into GMT
      ENDDO

      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
        DO N = 1, NT
          DO I2 = 1, NSTA1
            WRITE(44,'(2I5,F10.4,120F10.4)') N, I2, TIME(N),
     1      (zval(I2,K,N),K = 1, KB1)
          ENDDO
        ENDDO
        CLOSE(44)
      ENDIF

      DO NJ = 1, NSTA_CTL
        SDEPTH = SDEPTH_STA(NJ)    
        IF(DEPTH(NJ) .GT. 0.0) THEN
          IF(SDEPTH .GT. DEPTH(NJ)) SDEPTH = DEPTH(NJ)/2.0
        ENDIF	 

        OPEN(10,FILE = TRIM(FILESHORT(NJ)))
        DO N = 1, NT
          JDAY = TIME(N) + JDAY0
          CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
          IYEAR = INT(YEARB)
          ICM   = INT(MONTHB+0.001)
          ICD   = INT(DAYB+0.001)
          IHR   = INT(HOURB+0.001)
          IMN0  = INT((HOURB-IHR)*60+0.1)
	  IMN   = INT(IMN0/DELT_MIN+0.5)*DELT_MIN

          YEARB  = IYEAR
          MONTHB = ICM
          DAYB   = ICD
          HOURB  = IHR + IMN/60.0
          JDAY = JULIAN(YEARB,MONTHB,DAYB,HOURB)
          CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
          IYEAR = INT(YEARB+0.001)
          ICM   = INT(MONTHB+0.001)
          ICD   = INT(DAYB+0.001)
          IHR   = INT(HOURB+0.001)
          IMN0  = INT((HOURB-IHR)*60+0.1)
	  IMN   = INT(IMN0/DELT_MIN+0.5)*DELT_MIN
          DAYJ  = JDAY - JDAY0 + 1.0
          ELE = H(NJ,N)
          DTMP = DEPTH(NJ)

          IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') THEN
            CALL SIGMA2Z_ROMS_FIX(SIGMA,DTMP,ELE,KB,ZSIGMA,
     1           HC,THETA_S,THETA_B,TCLINE)

          ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'POM') THEN
            DO K = 1, KB
              IF(SIGMA(K) .GT. 0.0) SIGMA(K) = -SIGMA(K)
            ENDDO
            CALL SIGMA2Z_POM(SIGMA,DTMP,ELE,KB,ZSIGMA)  
          ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') THEN
            DO K = 1, KB
	      SIGMA(K) = SIGLAY(NJ,K)
              IF(SIGMA(K) .GT. 0.0) SIGMA(K) = -SIGMA(K)
            ENDDO
            CALL SIGMA2Z_POM(SIGMA,DTMP,ELE,KB,ZSIGMA)  
          ELSEIF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
            DO K = 1, KB
              ZSIGMA(K) = ZVAL(NJ,K,N)
            ENDDO
          ELSE
            DO K = 1, KB
              ZSIGMA(K) = SIGMA(K)
            ENDDO
          ENDIF

          IF(KINDAT .eq. 1) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK)=UK(KKK)
              ENDDO
              DO KKK = 1, KB
                UK(KKK) = U(NJ,KB-KKK+1,N)
                VK(KKK) = V(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK) = U(NJ,KKK,N)
                VK(KKK) = V(NJ,KKK,N)
              ENDDO
            ENDIF

            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
              VTMP = VK(1)
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP = UK(KB)
              VTMP = VK(KB)
            ELSE
              DO K9 = 1, KB-1
                IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1             (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                  X1 = ZSIGMA(K9)
                  X2 = ZSIGMA(K9+1)
                  Y1 = UK(K9)
                  Y2 = UK(K9+1)
                  CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
                  Y1 = VK(K9)
                  Y2 = VK(K9+1)
                  CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,VTMP)
                  EXIT
                ENDIF
              ENDDO
            END IF
            CALL VELDIR(VTMP,UTMP,AANGLE,AVEL)
            IF(AVEL .LE. 0.0) AANGLE = 0.0
            IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,
     1        AVEL,AANGLE,UTMP,VTMP

          ELSE IF(KINDAT .EQ. 2) THEN
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,H(NJ,N)
          ELSE IF(KINDAT .EQ. 3) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK) = UK(KKK)
              ENDDO
              DO KKK = 1, KB
                UK(KKK) = T(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK) = T(NJ,KKK,N)
              ENDDO
            ENDIF

            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP = UK(KB)
            ELSE
              DO K9 = 1, KB-1
                IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1             (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                  X1 = ZSIGMA(K9)
                  X2 = ZSIGMA(K9+1)
                  Y1 = UK(K9)
                  Y2=UK(K9+1)
                  CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
15          WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP

          ELSE IF(KINDAT .EQ. 4) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK) = UK(KKK)
              ENDDO
              DO KKK = 1, KB
                UK(KKK) = S(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK) = S(NJ,KKK,N)
              ENDDO
            ENDIF
            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP=UK(KB)
            ELSE
              DO K9 = 1, KB-1
                IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1             (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                  X1 = ZSIGMA(K9)
                  X2 = ZSIGMA(K9+1)
                  Y1 = UK(K9)
                  Y2 = UK(K9+1)
                  CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
                  EXIT
                ENDIF
              ENDDO
            END IF
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP
          ENDIF
        ENDDO
        CLOSE(10)
      ENDDO
      CLOSE(9)
      WRITE(*,*) 'Concatenation of nowcast files completed !!!'
100   FORMAT(F10.5,I5,4I3,50F10.4)

      STOP
      END

