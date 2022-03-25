C     PROGRAM read_netcdf.f
C f90 read_netcdf_modeltides.f -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include \
C -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf -o read_netcdf_modeltides.x
C lf95 read_netcdf_modeltides.f -I/disks/NASPUB/usr/Local/Linux/netcdf/include \
C    -L/disks/NASPUB/usr/Local/Linux/netcdf/lib -lnetcdf -o ../binLinux/read_netcdf_modeltides.x
C   call subroutines 
C      first_read_netcdf
C      read_netcdf
C      GREGORIAN
C      VELDIR
C      linear
C      spline
C  Modified by Lianyuan Zheng on 03/01/2017

      PROGRAM READ_NETCDF_MODELTIDES
      include 'netcdf.inc'
      CHARACTER*50  SIGMANAME
      CHARACTER*200 FNAME,FCTL,FILESHORT,LONGNAMES,STAID
      CHARACTER*200 OCEAN_MODEL,ANAME,BUFFER
      REAL*8 JDAY,JDAY0,JDAY1,JBASE_DATE,JULIAN
      REAL*8 YEARB,MONTHB,DAYB,HOURB
      REAL*4, ALLOCATABLE ::  TIME(:)
      REAL*4, ALLOCATABLE ::  LON(:)
      REAL*4, ALLOCATABLE ::  LAT(:)
      REAL*4, ALLOCATABLE ::  ANGLE(:)
      REAL*4, ALLOCATABLE ::  H(:,:)
      REAL*4, ALLOCATABLE ::  U(:,:,:)
      REAL*4, ALLOCATABLE ::  V(:,:,:)
      REAL*4, ALLOCATABLE ::  SPEED(:,:,:)
      REAL*4, ALLOCATABLE ::  DIR(:,:,:)
      REAL*4, ALLOCATABLE ::  T(:,:,:)
      REAL*4, ALLOCATABLE ::  S(:,:,:)
      REAL*4, ALLOCATABLE ::  DEPTH(:)
      REAL*4, ALLOCATABLE ::  SIGMA(:)
      REAL*4, ALLOCATABLE ::  ZSIGMA(:)
      REAL*4, ALLOCATABLE ::  UK(:)
      REAL*4, ALLOCATABLE ::  VK(:)
      REAL*4, ALLOCATABLE ::  SIGLAY(:,:)
      REAL*4, ALLOCATABLE ::  ZVAL(:,:,:)
      REAL*4, ALLOCATABLE ::  TMP3D(:,:,:)
      INTEGER BASE_DATE(4)

      CALL GETARG(1,BUFFER)
      READ(BUFFER,*) IYRS, IMMS, IDDS, IHHS, MNS
      CALL GETARG(2,FNAME)
      CALL GETARG(3,FCTL)
      CALL GETARG(4,BUFFER)
      READ(BUFFER,*) KINDAT
      CALL GETARG(5,OCEAN_MODEL)
      CALL GETARG(6,BUFFER)
      READ(BUFFER,*) DELT_MIN
C      IF(TRIM(OCEAN_MODEL) .EQ. 'ADCIRC') THEN
C        WRITE(*,*) 'ADCIRC is not included in this SA package yet!'
C        STOP
C      END IF
      
      YEARB  = IYRS
      MONTHB = 1.0
      DAYB   = 1.0
      HOURB  = 0.0
      JBASE_DATE = JULIAN(YEARB,MONTHB,DAYB,HOURB)

      YEARB  = IYRS
      MONTHB = IMMS
      DAYB   = IDDS
      HOURB  = IHHS + MNS/60.0
      JDAY0  = JULIAN(YEARB,MONTHB,DAYB,HOURB)

      WRITE(*,*) FNAME
      CALL FIRST_READ_NETCDF(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,
     1  BASE_DATE,SIGMANAME,THETA_S,THETA_B,TCLINE,HC,OCEAN_MODEL)
      WRITE(*,*) 'NSTA = ', NSTA,'NT = ', NT,'KB = ', KB

      ALLOCATE(TIME(NT))
      ALLOCATE(LON(NSTA))
      ALLOCATE(LAT(NSTA))
      ALLOCATE(DEPTH(NSTA))
      ALLOCATE(SIGMA(KB))
      ALLOCATE(H(NSTA,NT))
      ALLOCATE(U(NSTA,KB,NT))
      ALLOCATE(V(NSTA,KB,NT))
      ALLOCATE(T(NSTA,KB,NT))
      ALLOCATE(S(NSTA,KB,NT))

      ALLOCATE(SPEED(NSTA,KB,NT))
      ALLOCATE(DIR(NSTA,KB,NT))
      ALLOCATE(ZSIGMA(KB))
      ALLOCATE(UK(KB))
      ALLOCATE(VK(KB))

      YEARB  = DBLE(BASE_DATE(1)*1.0)
      MONTHB = DBLE(BASE_DATE(2)*1.0)
      DAYB   = DBLE(BASE_DATE(3)*1.0)
      HOURB  = DBLE(BASE_DATE(4)*1.0)
      JDAY1 = JULIAN(YEARB, MONTHB, DAYB, HOURB)
      IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM' ) THEN
        IF(ALLOCATED(SIGLAY)) DEALLOCATE(SIGLAY)
        ALLOCATE(SIGLAY(NSTA,KB))
        STATUS = NF_OPEN(TRIM(FNAME), NF_NOWRITE, NCID)
        STATUS = NF_INQ_VARID(NCID,'siglay',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,SIGLAY)
        ENDIF
	DO K = 1, KB
	  WRITE(55,'(I5,200F9.4)') K,(SIGLAY(N,K),N=1,20)
	ENDDO
	CLOSE(55)  
        STATUS = NF_CLOSE(NCID)
      ENDIF

      CALL READ_NETCDF(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,TIME,
     1  H,U,V,T,S,DEPTH,SIGMA,LON,LAT,KINDAT,OCEAN_MODEL)

      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE' ) THEN
        IF(ALLOCATED(ZVAL) ) DEALLOCATE(ZVAL)
        ALLOCATE(ZVAL(NSTA,KB,NT))
        IF(ALLOCATED(TMP3D) ) DEALLOCATE(TMP3D)
        ALLOCATE(TMP3D(NSTA,KB,NT))
        STATUS = NF_OPEN(TRIM(FNAME), NF_NOWRITE, NCID)
        STATUS = NF_INQ_VARID(NCID,'zval',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,TMP3D)
        ENDIF
        STATUS=NF_CLOSE(NCID)

C  Reverse vertical coordinates (K=1 for bottom and K=KB for surface)
C  zeta=-9999. is dry cell for SELFE
        DO NS = 1, NSTA   !!
          DO N = 1, NT
            DO K = 1, KB
	      K1 = KB - K + 1
              ZVAL(NS,K,N) = TMP3D(NS,K1,N)
            ENDDO
	  ENDDO
	ENDDO  
      ENDIF

      OPEN(9,FILE = TRIM(FCTL))
50    READ(9,*,ERR=999,END=999) STAID, FILESHORT, LONGNAMES
C  Read standard depth and total depth
      READ(9,*) ALAT, ALON, DIRFLOOD, SDEPTH, TDEPTH 

      DISTMIN = 99999999.0
      DO NSS = 1, NSTA
        IF(ALON .GT. 180.0) ALON = ALON - 360.0
        IF(LON(NSS) .GT. 180.0)LON(NSS) = LON(NSS) - 360.0
        YLAT = LAT(NSS)
        XLON = LON(NSS)
        IF((ABS(YLAT) .LT. 100.) .OR. (ABS(XLON) .LT. 360.) ) THEN
          CALL DIST(YLAT,XLON,ALAT,ALON,DIS)    
          IF(DIS .LE. DISTMIN) THEN
            DISTMIN = DIS
            NODE_FIND = NSS
          ENDIF
        ENDIF  
      ENDDO

      NJ = NODE_FIND
      WRITE(6,*)  TRIM(STAID),'  ', TRIM(FILESHORT),'  ',NJ,
     1   DEPTH(NJ), SDEPTH
      OPEN(10,FILE = TRIM(FILESHORT))
      DO N = 1, NT
        JDAY = TIME(N) + JDAY1
        IF(JDAY .GE. JDAY0) THEN
          CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)
          IYEAR  = INT(YEARB)
          ICM    = INT(MONTHB+0.001)
          ICD    = INT(DAYB+0.001)
          IHR    = INT(HOURB+0.001)
          IMN0   = INT((HOURB-IHR)*60+0.1)
	  IMN    = INT(IMN0/DELT_MIN+0.5)*DELT_MIN
          YEARB  = DBLE(IYEAR*1.0)
          MONTHB = DBLE(ICM*1.0)
          DAYB   = DBLE(ICD*1.0)
          HOURB  = DBLE(IHR*1.0 + IMN/60.0)
          JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB)
          CALL GREGORIAN(JDAY,YEARB,MONTHB,DAYB,HOURB)

          IYEAR  = INT(YEARB+0.001)
          ICM    = INT(MONTHB+0.001)
          ICD    = INT(DAYB+0.001)
          IHR    = INT(HOURB+0.001)
          IMN0   = INT((HOURB-IHR)*60+0.1)
	  IMN    = INT(IMN0/DELT_MIN+0.5)*DELT_MIN
          DAYJ   = JDAY - JBASE_DATE + 1.0

          IF(TRIM(OCEAN_MODEL) .EQ. 'ADCIRC' )THEN
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,H(NJ,N)
            CYCLE
          ENDIF          

          ELE  = H(NJ,N)
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
          END IF

          IF(KINDAT .EQ. 1) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK) = UK(KKK)
                UK(KKK) = U(NJ,KB-KKK+1,N)
                VK(KKK) = V(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK)=u(NJ,KKK,N)
                VK(KKK)=V(NJ,KKK,N)
              ENDDO
            ENDIF

            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
              VTMP = VK(1)
              GOTO 5
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP = UK(KB-1)
              VTMP = VK(KB-1)
              GOTO 5
            ENDIF

            DO K9 = 1, KB-1
              IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1           (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                GOTO 333
              ENDIF
            ENDDO

333         X1 = ZSIGMA(K9)
            X2 = ZSIGMA(K9+1)
            Y1 = UK(K9)
            Y2 = UK(K9+1)
            CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
            Y1 = VK(K9)
            Y2 = VK(K9+1)
            CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,VTMP)

5           CALL VELDIR (VTMP,UTMP,AANGLE,AVEL)
            IF(AVEL .LE. 0.0) AANGLE = 0.0
            IF(AANGLE .GT. 360.0) AANGLE = AANGLE - 360.0
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,
     1         AVEL,AANGLE,UTMP,VTMP

          ELSE IF(KINDAT .EQ. 2) THEN
            WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,h(NJ,N)

          ELSE IF(KINDAT .EQ. 3) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK) = UK(KKK)
                UK(KKK) = T(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK) = T(NJ,KKK,N)
              ENDDO
            ENDIF

            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
              GOTO 15
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP = UK(KB-1)
              GOTO 15
            ENDIF

            DO K9 = 1, KB-1
              IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1           (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                GOTO 334
              ENDIF
            ENDDO

334         X1 = ZSIGMA(K9)
            X2 = ZSIGMA(K9+1)
            Y1 = UK(K9)
            Y2 = UK(K9+1)
            CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
15          WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP

          ELSE IF(KINDAT .EQ. 4) THEN
            IF(ZSIGMA(1) .GT. ZSIGMA(KB)) THEN
              DO KKK = 1, KB
                UK(KKK) = ZSIGMA(KB-KKK+1)
              ENDDO
              DO KKK = 1, KB
                ZSIGMA(KKK) = UK(KKK)
                UK(KKK) = S(NJ,KB-KKK+1,N)
              ENDDO
            ELSE     
              DO KKK = 1, KB
                UK(KKK) = S(NJ,KKK,N)
              ENDDO
            ENDIF

            IF(SDEPTH .LT. ZSIGMA(1)) THEN
              UTMP = UK(1)
              GOTO 25
            ELSEIF(SDEPTH .GT. ZSIGMA(KB)) THEN
              UTMP = UK(KB-1)
              GOTO 25
            ENDIF

            DO K9 = 1, KB-1
              IF((SDEPTH .GE. ZSIGMA(K9)) .AND.
     1           (SDEPTH .LT. ZSIGMA(K9+1))) THEN
                GOTO 335
              ENDIF
            ENDDO

335         X1 = ZSIGMA(K9)
            X2 = ZSIGMA(K9+1)
            Y1 = UK(K9)
            Y2 = UK(K9+1)
            CALL LINEAR(X1,Y1,X2,Y2,SDEPTH,UTMP)
25          WRITE(10,100) DAYJ,IYEAR,ICM,ICD,IHR,IMN,UTMP
          ENDIF
        ENDIF
      ENDDO
      CLOSE(10)
      GOTO 50
999   CONTINUE
      CLOSE(9)
100   FORMAT(F10.5,I5,4I3,50F10.4)

      STOP
      END

