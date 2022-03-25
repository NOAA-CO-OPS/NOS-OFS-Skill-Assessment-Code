CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE first_read_netcdf(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,
     1  base_date,sigmaname,theta_s,theta_b,Tcline,hc,OCEAN_MODEL)

      include 'netcdf.inc'
      character*200 fname,buffer,OCEAN_MODEL,ANAME
      integer NCID,STATUS,NVARS,NGATTS,UNLIMDIMID
      integer IDVAR,COUNT(3)
      integer dimids(6),ndims,NT
      integer base_date(4)
      character*50 cc,sigmaunits,sigmaname
      LOGICAL ROMSMODEL

C  The following variables are used for standard ROMS sigma coordinate
      ROMSMODEL = .FALSE.
      THETA_S = 0.0
      THETA_B = 0.0
      TCLINE  = 0.0
      HC      = 0.0
      NSTA    = 0
      KB      = 0
      MESHDIM = 1
      NT      = 1
      NCHAR   = 1

      STATUS = NF_OPEN(FNAME, NF_NOWRITE, NCID)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_INQ(NCID, NDIMS, NVARS, NGATTS, UNLIMDIMID)
      ENDIF
      
C   NDIMS=number of dimension parameters
C   NVARS=number of total variables in the netcdf file
C   NGATTS= number of global attributes
C   UNLIMDIMID= dimension ID which is unlimited.   
      DO I = 1, NDIMS
        STATUS = NF_INQ_DIM(NCID, I, CC, ILATID)  !! extract dimension name
        STATUS = NF_INQ_DIMLEN(NCID, I, ILATID)
        IF(TRIM(CC) .EQ. 'station') THEN
          NSTA = ILATID
        ELSEIF((TRIM(CC) .EQ. 'sigma') 
     1    .OR. (TRIM(CC) .EQ. 's_rho')
     2    .OR. (TRIM(CC) .EQ. 'nvrt')
     3    .OR. (TRIM(CC) .EQ. 'siglay')) THEN
          KB = ILATID
        ELSEIF(TRIM(CC) .EQ. 'meshdim') THEN
          MESHDIM = ILATID
        ELSEIF((TRIM(CC) .EQ. 'time') .OR. 
     1         (TRIM(CC) .EQ. 'ocean_time')) THEN
          NT = ILATID
        ELSEIF((TRIM(CC) .EQ. 'charlength')
     1    .OR. (trim(cc) .eq. 'namelen')
     1    .OR. (TRIM(CC) .EQ. 'clen')) THEN
          NCHAR = ILATID
        ENDIF
      ENDDO

      ANAME='time'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') ANAME = 'ocean_time'
      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_ATT_TEXT(NCID,IDVAR,'units',BUFFER)
        IF(STATUS .EQ. NF_NOERR) THEN
	  CC   = TRIM(adjustL(BUFFER))
          LEN1 = LEN_TRIM(CC)
          LL   = INDEX(TRIM(CC),'since',BACK=.TRUE.)
          READ(CC(LL+6:LEN1),'(I4,3(1x,I2))') IYR, IMM, IDD, IHH
          IF((IYR .LT. 100) .AND. (IYR .GT. 50) ) THEN
            IYR = IYR + 1900
          ELSEIF((IYR .LT. 100) .AND. (IYR .LE. 50) ) THEN
            IYR = IYR + 2000
          ENDIF
          BASE_DATE(1) = IYR
          BASE_DATE(2) = IMM
          BASE_DATE(3) = IDD
          BASE_DATE(4) = IHH
        ENDIF
      ELSE 
        WRITE(*,*)'Reading base_date failed, and stop ...'
	STOP
      ENDIF

C      TIME_SCALE = 1.0
C      LEN1 = LEN_TRIM(CC)
C      LL = INDEX(TRIM(CC), 'days', BACK = .TRUE.)
C      IF(LL .GT. 0) TIME_SCALE = 1.0
C      LL = INDEX(TRIM(CC), 'hours')
C      IF(LL .GT. 0) TIME_SCALE = 1.0/24.0
C      LL=INDEX(TRIM(CC), 'seconds')
C      IF(LL .GT. 0) TIME_SCALE = 1.0/86400.0

      IF(trim(OCEAN_MODEL) .EQ. 'ADCIRC') GOTO 555 
      IF(TRIM(OCEAN_MODEL) .EQ. 'POM')   ANAME = 'sigma'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS')  ANAME = 's_rho'
      IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') ANAME = 'siglay'
      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') ANAME = 'zval'
      
      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_ATT_TEXT(NCID,IDVAR,'standard_name',sigmaname)
      ENDIF

      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') THEN
        STATUS = NF_INQ_VARID(NCID,'theta_s',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,theta_s)
        STATUS = NF_INQ_VARID(NCID,'theta_b',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,theta_b)
        STATUS = NF_INQ_VARID(NCID,'Tcline',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,Tcline)
        STATUS = NF_INQ_VARID(NCID,'hc',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,hc)
        ROMSMODEL = .TRUE.
      ENDIF 

      IF(TRIM(sigmaname(1:19)) == 'ocean_sz_coordinate') THEN
        STATUS = NF_INQ_VARID(NCID,'h_s',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,cw_1)
        STATUS = NF_INQ_VARID(NCID,'h_c',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,cw_2)
        STATUS = NF_INQ_VARID(NCID,'theta_b',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,h_cut)
        STATUS = NF_INQ_VARID(NCID,'theta_f',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,c_cut)
      ENDIF 

555   STATUS = NF_CLOSE(NCID)
      RETURN
      END       


C       This subroutine read cycle nowcast or forecast      
      SUBROUTINE read_netcdf(NSTA,NT,KB,NCHAR,MESHDIM,FNAME,time,zeta
     1 ,u,v,temp,salt,depth,sigma,lon,lat,KINDAT,OCEAN_MODEL)

      include 'netcdf.inc'
      character*200 fname,cc,sigmaname,OCEAN_MODEL,ANAME,buffer
      integer stationij(NSTA,meshdim)
      character stationnames(NCHAR,NSTA)
      real*4 time(NT),lon(NSTA),lat(NSTA),depth(NSTA),sigma(KB)
      real*4 zeta(NSTA,NT),u(NSTA,KB,NT),v(NSTA,KB,NT)
      real*4 temp(NSTA,KB,NT),salt(NSTA,KB,NT)
      real*4 angle(NSTA),tmp1(NT)
      real*4 xtmp(KB,NSTA,NT),tmp3d(NSTA,KB,NT)
      real*8 jday,jday0,jday1,jbase_date,JULIAN
      real*8 yearb,monthb,dayb,hourb

      integer NCID,STATUS,NVARS,NGATTS,UNLIMDIMID
      integer IDVAR,COUNT(3)
      integer dimids(6),ndims,NT
      integer base_date(4)
      LOGICAL ROMSMODEL

      ROMSMODEL = .FALSE.
      STATUS = NF_OPEN(fname, NF_NOWRITE, NCID)
      IF(STATUS .NE. NF_NOERR) THEN
        WRITE(6,*) 'Problem opening netcdf file'
        STOP
      ENDIF

      IF(TRIM(OCEAN_MODEL) .EQ. 'POM')   ANAME = 'stationnames'
      IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') ANAME = 'name_station'
      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') ANAME = 'name_station'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ADCIRC') ANAME = 'station_name'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') THEN
        ANAME = 'stationnames'
        ROMSMODEL = .TRUE.
      END IF

      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_VAR_TEXT(NCID,IDVAR,stationnames)
      ENDIF

      ANAME = 'time'
      IF(trim(OCEAN_MODEL) .EQ. 'ROMS') ANAME = 'ocean_time'
      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_ATT_TEXT(NCID,IDVAR,'units',BUFFER)
        IF(STATUS .EQ. NF_NOERR) THEN
	  CC = TRIM(adjustL(BUFFER))
        ENDIF
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,TIME)
      ENDIF

      TIME_SCALE = 1.0
      LEN1 = LEN_TRIM(CC)
      LL = INDEX(TRIM(CC),'days',BACK=.TRUE.)
      IF(LL .GT. 0) TIME_SCALE = 1.0
      LL = INDEX(TRIM(CC), 'hours')
      IF(LL .GT. 0) TIME_SCALE = 1.0/24.0
      LL = INDEX(TRIM(CC), 'minutes')
      IF(LL .GT. 0) TIME_SCALE = 1.0/1440.0
      LL=INDEX(TRIM(CC),'seconds')
      IF(LL .GT. 0) TIME_SCALE = 1.0/86400.0

      DO N = 1, NT
        TIME(N) = TIME(N) * TIME_SCALE
      ENDDO

      ANAME='zeta'
      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') ANAME = 'zeta_adj'
      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      STATUS = NF_GET_VAR_REAL(NCID,IDVAR,zeta)

      IF(TRIM(OCEAN_MODEL) .EQ. 'ADCIRC') THEN
        STATUS = NF_INQ_VARID(NCID,'x',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lon)
        END IF
        STATUS = NF_INQ_VARID(NCID,'y',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lat)
        END IF
        GOTO 560      
      ENDIF

      STATUS = NF_INQ_VARID(NCID,'lon',IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lon)
      ELSE
        STATUS = NF_INQ_VARID(NCID,'lon_rho',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lon)
      ENDIF

      STATUS = NF_INQ_VARID(NCID,'lat',IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lat)
      ELSE
        STATUS = NF_INQ_VARID(NCID,'lat_rho',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lat)
      ENDIF

      STATUS = NF_INQ_VARID(NCID,'depth',IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
         STATUS = NF_GET_VAR_REAL(NCID,IDVAR,depth)
      ELSE
         STATUS = NF_INQ_VARID(NCID,'h',IDVAR)
         STATUS = NF_GET_VAR_REAL(NCID,IDVAR,depth)
      ENDIF 

      STATUS = NF_INQ_VARID(NCID,'u',IDVAR)
C  Added by Zheng
      IF(STATUS .EQ. NF_NOERR) THEN
C  Added by Zheng
      IF((.NOT.ROMSMODEL) .OR. (TRIM(OCEAN_MODEL) .NE. 'ROMS')) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,tmp3d)
        IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
! For SELFE, reverse vertical coordinates because K=1 for bottom and K=KB for surface
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
  	        K1 = KB - K + 1
                u(NS,K,N) = tmp3d(NS,K1,N)
              ENDDO
            ENDDO
          ENDDO
	ELSE
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1,KB
                u(NS,K,N) = tmp3d(NS,K,N)
              ENDDO
            ENDDO
          ENDDO
	ENDIF          
      ELSE
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,xtmp)
        DO NS = 1, NSTA
          DO N = 1, NT
            DO K = 1, KB
              u(NS,K,N) = xtmp(K,NS,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C  Added by Zheng
      END IF
C  Added by Zheng

      STATUS = NF_INQ_VARID(NCID,'v',IDVAR)
C  Added by Zheng
      IF(STATUS .EQ. NF_NOERR) THEN
C  Added by Zheng
      IF((.NOT.ROMSMODEL) .OR. (trim(OCEAN_MODEL) .NE. 'ROMS')) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,tmp3d)
        IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
  	        K1 = KB - K + 1
                v(NS,K,N) = tmp3d(NS,K1,N)
              ENDDO
            ENDDO
          ENDDO
	ELSE
          DO NS = 1 ,NSTA
            DO N = 1, NT
              DO K = 1, KB
                v(NS,K,N) = tmp3d(NS,K,N)
              ENDDO
            ENDDO
          ENDDO
	ENDIF 
      ELSE
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,xtmp)
        DO NS = 1, NSTA
          DO N = 1, NT
            DO K = 1, KB
              v(NS,K,N) = xtmp(K,NS,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C  Added by Zheng
      END IF
C  Added by Zheng

      STATUS = NF_INQ_VARID(NCID,'temp',IDVAR)
C  Added by Zheng
      IF(STATUS .EQ. NF_NOERR) THEN
C  Added by Zheng
      IF((.NOT.ROMSMODEL) .OR. (trim(OCEAN_MODEL) .NE. 'ROMS'))THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,tmp3d)
        IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
  	        K1 = KB - K + 1
                temp(NS,K,N) = tmp3d(NS,K1,N)
              ENDDO
            ENDDO
          ENDDO
	ELSE
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
                temp(NS,K,N) = tmp3d(NS,K,N)
              ENDDO
            ENDDO
          ENDDO
	ENDIF 
      ELSE
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,xtmp)
        DO NS = 1, NSTA
          DO N = 1, NT
            DO K = 1, KB
              temp(NS,K,N) = xtmp(K,NS,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C  Added by Zheng
      END IF
C  Added by Zheng

      IF(TRIM(OCEAN_MODEL) .EQ. 'POM')   ANAME = 'salt'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS')  ANAME = 'salt'
      IF(TRIM(OCEAN_MODEL) .EQ. 'FVCOM') ANAME = 'salinity'
      IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') ANAME = 'salinity'
      STATUS = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
C  Added by Zheng
      IF(STATUS .EQ. NF_NOERR) THEN
C  Added by Zheng
      IF((.NOT.ROMSMODEL) .OR. (trim(OCEAN_MODEL) .NE. 'ROMS')) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,tmp3d)
        IF(TRIM(OCEAN_MODEL) .EQ. 'SELFE') THEN
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
                K1 = KB - K + 1
                salt(NS,K,N) = tmp3d(NS,K1,N)
              ENDDO
            ENDDO
          ENDDO
	ELSE
          DO NS = 1, NSTA
            DO N = 1, NT
              DO K = 1, KB
                salt(NS,K,N) = tmp3d(NS,K,N)
              ENDDO
            ENDDO
          ENDDO
	ENDIF
      ELSE
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,xtmp)
        DO NS = 1, NSTA
          DO N= 1, NT
            DO K= 1, KB
              salt(NS,K,N) = xtmp(K,NS,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C  Added by Zheng
      END IF
C  Added by Zheng

      IF(ROMSMODEL) THEN
        STATUS = NF_INQ_VARID(NCID,'angle',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,angle)
          DO NS = 1, NSTA
           DO N = 1, NT
            DO K = 1, KB
             Utmp = u(NS,K,N)*cos(angle(NS))-V(NS,K,N)*sin(angle(NS))
             Vtmp = u(NS,K,N)*sin(angle(NS))+V(NS,K,N)*cos(angle(NS))
             U(NS,K,N) = Utmp
             V(NS,K,N) = Vtmp
            ENDDO
           ENDDO
          ENDDO
        ENDIF
      END IF

! reverse vertical coordinates because K=1 for bottom and K=KB for surface
! zeta=-9999. is dry cell for SELFE
      IF(trim(OCEAN_MODEL) .EQ. 'SELFE') THEN
        DO NS = 1, NSTA
          DO N = 1, NT
            DO K = 1, KB
	      IF(ZETA(NS,N) .LE. -100.0) THEN
	        U(NS,K,N)=-9999.0 
	        V(NS,K,N)=-9999.0 
	      ENDIF   
            ENDDO
          ENDDO
        ENDDO

!! get ride of all -9999.0 for the deep layers 
        DO NS = 1, NSTA
         DO N = 1, NT
          DO K = 2, KB
	   if(u(NS,K,N) .LE. -99.9) u(NS,K,N) = u(NS,K-1,N)
	   if(v(NS,K,N) .LE. -99.9) v(NS,K,N) = v(NS,K-1,N)
	   if(temp(NS,K,N) .LE. -99.9) temp(NS,K,N) = temp(NS,K-1,N)
	   if(salt(NS,K,N) .LE. -99.9) salt(NS,K,N) = salt(NS,K-1,N)
	  ENDDO 
         ENDDO
        ENDDO
      ENDIF  

      STATUS = NF_INQ_VARID(NCID,'sigma',IDVAR)
      IF(STATUS .EQ. NF_NOERR) THEN
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,sigma)
      ELSE
        STATUS=NF_INQ_VARID(NCID,'s_rho',IDVAR)
        IF(STATUS .EQ. NF_NOERR) THEN
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,sigma)
        ENDIF
      ENDIF

560   STATUS=NF_CLOSE(NCID)
      RETURN
      END
