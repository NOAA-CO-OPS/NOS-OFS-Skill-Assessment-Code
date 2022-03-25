C   Program read_netcdf_fcst1.f
C  Compile on opmarine
C  f90 read_netcdf_fcst.f -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include \
C  -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf -o ../binUnix/read_netcdf_fcst.x
C  Compile on Linux 
C  lf95 read_netcdf_fcst1.f read_netcdf_subs.f -I/disks/NASPUB/usr/Local/Linux/netcdf/include \
C  -L/disks/NASPUB/usr/Local/Linux/netcdf/lib -lnetcdf  -o ../binLinux/read_netcdf_fcst1.x
C
C  Modified by Lianyuan Zheng on 03/10/2017

      PROGRAM READ_NETCDF_FCST1
      include 'netcdf.inc'
      CHARACTER*200 FNAME,FCTL,OCEAN_MODEL,ANAME,BUFFER,CC
      REAL*4, ALLOCATABLE ::  TIME(:,:,:)
      REAL*4, ALLOCATABLE ::  TIMETMP(:)
      REAL*8 JDAY,JDAY1,JBASE_DATE,JULIAN
      REAL*8 YEARB,MONTHB,DAYB,HOURB
      INTEGER BASE_DATE(4),RETVAL

      OPEN(5,FILE = 'filetmp.ctl')
      READ(5,*) DELT, NCYCLE, KINDAT, NFILES, NFDURATION
      READ(5,'(A200)') FCTL
      READ(5,'(A200)') OCEAN_MODEL 
      DELT = DELT/60.0
      NCUT = INT(NFDURATION/DELT+0.1)

      READ(5,'(A200)') FNAME
      READ(5,*) IYR1,IMM1,IDD1,IHH1,IMIN1

      RETVAL = NF_OPEN(FNAME, NF_NOWRITE, NCID)
      IF(RETVAL .EQ. NF_NOERR) THEN
        RETVAL = NF_INQ(NCID, NDIMS, NVARS, NGATTS, UNLIMDIMID)
      ENDIF

      DO I = 1, NDIMS
        RETVAL = NF_INQ_DIM(NCID, I, CC, ILATID)
        RETVAL = NF_INQ_DIMLEN(NCID, I, ILATID)
        IF((TRIM(CC) .EQ. 'time') .OR. 
     1     (TRIM(CC) .EQ. 'ocean_time')) THEN
          NT = ILATID
        ENDIF
      ENDDO
      ALLOCATE(TIMETMP(NT))

      ANAME = 'time'
      IF(TRIM(OCEAN_MODEL) .EQ. 'ROMS') ANAME = 'ocean_time'
      RETVAL = NF_INQ_VARID(NCID,TRIM(ANAME),IDVAR)
      IF(RETVAL .EQ. NF_NOERR) THEN
        RETVAL = NF_GET_VAR_REAL(NCID,IDVAR,TIMETMP)
        RETVAL = NF_GET_ATT_TEXT(NCID,IDVAR,'units',BUFFER)
        IF(RETVAL .EQ. NF_NOERR) THEN
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
          TIMETMP(N) = TIMETMP(N) * TIME_SCALE
        ENDDO
      ELSE 
        WRITE(*,*) 'Reading base_date failed, and stop ...'
	STOP
      ENDIF
      RETVAL = NF_CLOSE(NCID)

      YEARB  = DBLE(BASE_DATE(1)*1.0)
      MONTHB = DBLE(BASE_DATE(2)*1.0)
      DAYB   = DBLE(BASE_DATE(3)*1.0)
      HOURB  = DBLE(BASE_DATE(4)*1.0)
      JBASE_DATE =J ULIAN(YEARB,MONTHB,DAYB,HOURB)

      YEARB  = DBLE(IYR1*1.0)
      MONTHB = DBLE(IMM1*1.0)
      DAYB   = DBLE(IDD1*1.0)
      HOURB  = DBLE(IHH1*1.0 + IMIN1/60.0)
      JDAY   = JULIAN(YEARB,MONTHB,DAYB,HOURB)
      TSTART9 = JDAY-JBASE_DATE
      TEND = TSTART9 + NFDURATION/24.0 + 0.001
      TSTART9 = TSTART9-0.001

      N9 = 0
      DO IZ1 = 1, NT
        IF((TIMETMP(IZ1) .GE. TSTART9) .AND.
     1     (TIMETMP(IZ1) .LE. TEND)) THEN
          N9 = N9 + 1
        ENDIF  
      ENDDO
C      WRITE(*,*) ' N9 = ', N9, ' NT = ', NT, ' NCUT = ', NCUT

C  The following two lines have to be added to get correct N9, without the two lines 
C  N9 is not correct.  It is strange for this problem. 
      If (N9 .LT. NCUT) THEN
        WRITE(86,'(A1)') 'F'
      ELSE IF(N9 .GE. NCUT) THEN
        WRITE(86,'(A1)') 'T'
      ENDIF
      DEALLOCATE(TIMETMP)

      STOP
      END

