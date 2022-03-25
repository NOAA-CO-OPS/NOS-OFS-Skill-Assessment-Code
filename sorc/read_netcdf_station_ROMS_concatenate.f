C  lf95 read_netcdf_station_ROMS_concatenate.f utility.f -I/usr/local/include -L/usr/local/lib -lnetcdf -o read_netcdf_station_ROMS_concatenate.x
C  NOTES: vertical layer for ROMS is: 
C   sigma(k) < 0.0, k=1 is bottom, k=kb is for surface
C   call subroutines 
!      first_read_netcdf
!      read_netcdf
!      GREGORIAN
!      VELDIR
!      linear
!      spline
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      character*100 FNAME,FCTL,BUFFER*100,sigmaunits*20,sigmaname*20
      character*50 fileshort,longnames,staid,cc
      real*8 jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      real*4, allocatable ::  time(:)
      real*4, allocatable ::  lon(:)
      real*4, allocatable ::  lat(:)
      real*4, allocatable ::  angle(:)
      real*4, allocatable ::  h(:,:)
      real*4, allocatable ::  u(:,:,:)
      real*4, allocatable ::  v(:,:,:)
      real*4, allocatable ::  speed(:,:,:)
      real*4, allocatable ::  dir(:,:,:)
      real*4, allocatable ::  t(:,:,:)
      real*4, allocatable ::  s(:,:,:)
      real*4, allocatable ::  depth(:)
      real*4, allocatable ::  sigma(:)
      real*4, allocatable ::  zsigma(:)
      real*4, allocatable ::  UK(:)
      real*4, allocatable ::  VK(:)

      real*4, allocatable ::  timef(:)
      real*4, allocatable ::  hf(:,:)
      real*4, allocatable ::  uf(:,:,:)
      real*4, allocatable ::  vf(:,:,:)
      real*4, allocatable ::  speedf(:,:,:)
      real*4, allocatable ::  dirf(:,:,:)
      real*4, allocatable ::  tf(:,:,:)
      real*4, allocatable ::  sf(:,:,:)


      integer base_date(4)
!      CALL GETARG(1,BUFFER)
      NMAX=50000
      READ(5,*)IYRS,IMMS,IDDS,IHHS,MNS
      READ(5,'(a100)')FCTL
      READ(5,*)KINDAT,NFILES
      READ(5,'(a100)')FNAME
      print *,trim(FNAME)
      STATUS = NF_OPEN(trim(FNAME), NF_NOWRITE, NCID)
      STATUS = NF_INQ(NCID,NDIMS,NVARS,NGATTS,UNLIMDIMID)
      DO I=1,NDIMS
         STATUS = NF_INQ_DIM(NCID,i,cc,ILATID)  !! extract dimension name
         STATUS = NF_INQ_DIMLEN(NCID,i,ILATID)
         IF (trim(cc) .eq. 'station')then
            NSTA=ILATID
         ELSEIF (trim(cc) .eq. 's_rho')then
            KB=ILATID
         ELSEIF (trim(cc) .eq. 'ocean_time')then
            NT=ILATID
         ENDIF
      ENDDO
      allocate(time(NMAX))
      allocate(h(NSTA,NMAX))
      allocate(u(KB,NSTA,NMAX))
      allocate(v(KB,NSTA,NMAX))
      allocate(speed(KB,NSTA,NMAX))
      allocate(dir(KB,NSTA,NMAX))
      allocate(t(KB,NSTA,NMAX))
      allocate(s(KB,NSTA,NMAX))
      allocate(depth(NSTA))
      allocate(sigma(KB))
      allocate(zsigma(KB))
      allocate(UK(KB))
      allocate(VK(KB))
      allocate(lon(NSTA))
      allocate(lat(NSTA))
      allocate(angle(NSTA))
      STATUS = NF_INQ_VARID(NCID,'s_rho',IDVAR)
      STATUS = NF_GET_ATT_TEXT(NCID,IDVAR,'standard_name',sigmaname)
      IF(trim(sigmaname(1:18)) == 'ocean_s_coordinate')THEN
           STATUS = NF_INQ_VARID(NCID,'theta_s',IDVAR)
           STATUS = NF_GET_VAR_REAL(NCID,IDVAR,theta_s)
           STATUS = NF_INQ_VARID(NCID,'theta_b',IDVAR)
           STATUS = NF_GET_VAR_REAL(NCID,IDVAR,theta_b)
           STATUS = NF_INQ_VARID(NCID,'Tcline',IDVAR)
           STATUS = NF_GET_VAR_REAL(NCID,IDVAR,Tcline)
           STATUS = NF_INQ_VARID(NCID,'hc',IDVAR)
           STATUS = NF_GET_VAR_REAL(NCID,IDVAR,hc)
           STATUS = NF_INQ_VARID(NCID,'s_rho',IDVAR)
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,sigma)
      ENDIF 
      STATUS = NF_INQ_VARID(NCID,'lon_rho',IDVAR)
      STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lon)
      STATUS = NF_INQ_VARID(NCID,'lat_rho',IDVAR)
      STATUS = NF_GET_VAR_REAL(NCID,IDVAR,lat)
      STATUS = NF_INQ_VARID(NCID,'h',IDVAR)
      STATUS = NF_GET_VAR_REAL(NCID,IDVAR,depth)
      STATUS = NF_INQ_VARID(NCID,'angle',IDVAR)
      STATUS = NF_GET_VAR_REAL(NCID,IDVAR,angle)
      STATUS=NF_CLOSE(NCID)
      print *,'dimension NSTA, KB and TIME= ',NSTA,KB,NT
      yearb=IYRS
      monthb=1.
      dayb=1.
      hourb=0.
      jbase_date=JULIAN(yearb,monthb,dayb,hourb)
      yearb=IYRS
      monthb=IMMS
      dayb=IDDS
      hourb=IHHS+MNS/60.0
      jday0=JULIAN(yearb,monthb,dayb,hourb)
      REWIND(5)
      READ(5,*)
      READ(5,*)
      READ(5,*)
      NAJ=0
      time0=0.0
      DO IAJ=1,NFILES
        READ(5,'(a100)')FNAME
        write(*,*)trim(FNAME)
        STATUS = NF_OPEN(FNAME, NF_NOWRITE, NCID)
        STATUS = NF_INQ(NCID,NDIMS,NVARS,NGATTS,UNLIMDIMID)
        DO I=1,NDIMS
           STATUS = NF_INQ_DIM(NCID,i,cc,ILATID)  !! extract dimension name
           STATUS = NF_INQ_DIMLEN(NCID,i,ILATID)
           IF (trim(cc) .eq. 'station')then
              NSTA=ILATID
           ELSEIF (trim(cc) .eq. 's_rho')then
              KB=ILATID
           ELSEIF (trim(cc) .eq. 'ocean_time')then
              NT=ILATID
           ENDIF
        ENDDO
        print *,'dimension NSTA, KB and TIME= ',NSTA,KB,NT
        allocate(timef(NT))
        allocate(hf(NSTA,NT))
        allocate(uf(KB,NSTA,NT))
        allocate(vf(KB,NSTA,NT))
        allocate(tf(KB,NSTA,NT))
        allocate(sf(KB,NSTA,NT))
        STATUS = NF_INQ_VARID(NCID,'ocean_time',IDVAR)
        STATUS = NF_GET_ATT_TEXT(NCID,IDVAR,'units',cc)
        read(cc(15:18),*)base_date(1)
        read(cc(20:21),*)base_date(2)
        read(cc(23:24),*)base_date(3)
        read(cc(26:27),*)base_date(4)
        print *,'basedate= ',
     1       base_date(1),base_date(2),base_date(3),base_date(4)
        yearb=base_date(1)
        monthb=base_date(2)
        dayb=base_date(3)
        hourb=base_date(4)
        jday1=JULIAN(yearb,monthb,dayb,hourb)
        STATUS = NF_INQ_VARID(NCID,'ocean_time',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,timef)
        STATUS = NF_INQ_VARID(NCID,'zeta',IDVAR)
        STATUS = NF_GET_VAR_REAL(NCID,IDVAR,hf)
        IF (KINDAT .eq. 1)THEN
          STATUS = NF_INQ_VARID(NCID,'u',IDVAR)
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,uf)
          STATUS = NF_INQ_VARID(NCID,'v',IDVAR)
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,vf)
        ELSEIF (KINDAT .eq. 3)THEN
          STATUS = NF_INQ_VARID(NCID,'temp',IDVAR)
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,tf)
        ELSEIF (KINDAT .eq. 4)THEN
          STATUS = NF_INQ_VARID(NCID,'salt',IDVAR)
          STATUS = NF_GET_VAR_REAL(NCID,IDVAR,sf)
        ENDIF
        STATUS=NF_CLOSE(NCID)
        DO N=1,NT
             dayj=timef(N)/86400.0+jday1-jbase_date
          IF (dayj .GT. time0)then
             NAJ=NAJ+1
             time(NAJ)=dayj
             DO NS=1,NSTA
               h(NS,NAJ)=hf(NS,N)
               IF (KINDAT .eq. 1)THEN
                 DO K=1,KB
                   u(K,NS,NAJ)=uf(K,NS,N)
                   v(K,NS,NAJ)=vf(K,NS,N)
                 ENDDO
               ELSEIF (KINDAT .eq. 3)THEN
                 DO K=1,KB
                   t(K,NS,NAJ)=tf(K,NS,N)
                 ENDDO
               ELSEIF (KINDAT .eq. 4)THEN
                 DO K=1,KB
                   s(K,NS,NAJ)=sf(K,NS,N)
                 ENDDO
               ENDIF
             ENDDO
             time0=dayj
          ENDIF
        ENDDO
        deallocate(timef)
        deallocate(hf)
        deallocate(uf)
        deallocate(vf)
        deallocate(tf)
        deallocate(sf)
      ENDDO
      NT=NAJ
      DO NSS=1,NSTA
        write(*,666)NSS,lon(NSS),lat(NSS),depth(NSS)
      ENDDO
666   format(I6,3f12.5)      	
      OPEN(9,file=trim(FCTL))
50    read(9,*,err=999,end=999)staid,fileshort,longnames
      read(9,*)alat,alon,dirflood,sdepth,tdepth  !! read standard depth and total depth 
      distmin=99999999.
      DO NSS=1,NSTA
           ylat=lat(NSS)
           xlon=lon(NSS)
           if((abs(ylat) .lt. 100.) .or. (abs(xlon) .lt. 360.) )then
             call dist(ylat,xlon,alat,alon,dis)    
             if(dis.le.distmin) then
   		 distmin=dis
   		 node_find=NSS
             endif
           endif  
      ENDDO
      NJ= node_find
      write(6,*)trim(staid),'  ',trim(fileshort),'  ',NJ
     1 ,depth(NJ),sdepth,angle(NJ)
      open(10,file=trim(fileshort))
      DO N=1,NT
          jday=time(N)+jbase_date
          IF (jday .GE. jday0)then
           call GREGORIAN(jday,yearb,monthb,dayb,hourb)
           IYEAR=INT(yearb)
           ICM=int(monthb+0.001)
           ICD=INT(dayb+0.001)
           IHR=INT(hourb+0.001)
           IMN=INT((hourb-IHR)*60+0.1)
           dayj=jday-jbase_date+1.
           ele=h(NJ,N)
           dtmp=depth(NJ)
           IF(trim(sigmaname(1:18)) == 'ocean_s_coordinate')THEN  !! ROMS
              CALL sigma2Z_ROMS_FIX(sigma,dtmp,ele,KB,zsigma
     1              ,hc,theta_s,theta_b,Tcline)
           endif
           IF (KINDAT .eq. 1)THEN
              IF(zsigma(1) .GT. zsigma(KB) )THEN
                 DO KKK=1,KB
                    UK(KKK)=zsigma(KB-kkk+1)
                 ENDDO
                 DO KKK=1,KB
                   zsigma(kkk)=UK(kkk)
 !                  UK(KKK)=u(NJ,KB-KKK+1,N)
 !                  VK(KKK)=V(NJ,KB-KKK+1,N)
                   UK(KKK)=
     1 u(KB-KKK+1,NJ,N)*cos(angle(NJ))-V(KB-KKK+1,NJ,N)*sin(angle(NJ))
                   VK(KKK)=
     1 u(KB-KKK+1,NJ,N)*sin(angle(NJ))+V(KB-KKK+1,NJ,N)*cos(angle(NJ))
                 ENDDO
              ELSE     
                 DO KKK=1,KB
        !           UK(KKK)=u(NJ,KKK,N)
        !           VK(KKK)=V(NJ,KKK,N)
                   UK(KKK)=
     1 u(KKK,NJ,N)*cos(angle(NJ))-V(KKK,NJ,N)*sin(angle(NJ))
                   VK(KKK)=
     1 u(KKK,NJ,N)*sin(angle(NJ))+V(KKK,NJ,N)*cos(angle(NJ))
                 ENDDO
              ENDIF
              IF (sdepth .LT. zsigma(1) )THEN
                UTMP=UK(1)
                VTMP=VK(1)
                GOTO 5
              ELSEIF(sdepth .GT. zsigma(KB) )THEN
                UTMP=UK(KB)
                VTMP=VK(KB)
                GOTO 5
              ENDIF
!              IF (KB .LT. 5)THEN
                DO K9=1,KB-1
                  IF( (sdepth .GE. zsigma(K9)) .and.
     1                 (sdepth .LT. zsigma(K9+1)) )THEN
                    goto 333
                  ENDIF
                ENDDO
333             X1=zsigma(K9)
                X2=zsigma(K9+1)
                Y1=UK(K9)
                Y2=UK(K9+1)
                call linear(X1,Y1,X2,Y2,sdepth,UTMP)
                Y1=VK(K9)
                Y2=VK(K9+1)
                call linear(X1,Y1,X2,Y2,sdepth,VTMP)
!              ELSE IF (KB .GE. 5)THEN
!                CALL spline(KB,zsigma,UK,sdepth,UTMP)
!                CALL spline(KB,zsigma,VK,sdepth,VTMP)
!              ENDIF
5             CALL VELDIR (VTMP,UTMP,AANGLE,AVEL)
              IF (AVEL .LE. 0.0)   AANGLE = 0.
              IF (AANGLE.GT.360.)   AANGLE = AANGLE - 360.
              WRITE(10,100)dayj,IYEAR,ICM,ICD,IHR,IMN,
     1         AVEL,AANGLE,UTMP,VTMP
           ELSE IF(KINDAT .eq. 2)THEN
              WRITE(10,100)dayj,IYEAR,ICM,ICD,IHR,IMN,h(NJ,N)
           ELSE IF(KINDAT .eq. 3)THEN
              IF(zsigma(1) .GT. zsigma(KB) )THEN
                 DO KKK=1,KB
                    UK(KKK)=zsigma(KB-kkk+1)
                 ENDDO
                 DO KKK=1,KB
                   zsigma(kkk)=UK(KKK)
                   UK(KKK)=t(KB-KKK+1,NJ,N)
                 ENDDO
              ELSE     
                 DO KKK=1,KB
                   UK(KKK)=t(KKK,NJ,N)
                 ENDDO
              ENDIF
              IF (sdepth .LT. zsigma(1) )THEN
                UTMP=UK(1)
                GOTO 15
              ELSEIF(sdepth .GT. zsigma(KB) )THEN
                UTMP=UK(KB-1)
                GOTO 15
              ENDIF
!              IF (KB .LT. 5)THEN
                DO K9=1,KB-1
                  IF( (sdepth .GE. zsigma(K9)) .and.
     1                 (sdepth .LT. zsigma(K9+1)) )THEN
                    goto 334
                  ENDIF
                ENDDO
334             X1=zsigma(K9)
                X2=zsigma(K9+1)
                Y1=UK(K9)
                Y2=UK(K9+1)
                call linear(X1,Y1,X2,Y2,sdepth,UTMP)
!              ELSE IF (KB .GE. 5)THEN
!                CALL spline(KB,zsigma,UK,sdepth,UTMP)   !! commentted out 09/21/2006 since getting bad interpolated value from spline sometimes 
!              ENDIF
15          continue
!            if(dayj .GT. 212)
                 WRITE(10,100)dayj,IYEAR,ICM,ICD,IHR,IMN,UTMP
           ELSE IF(KINDAT .eq. 4)THEN
              IF(zsigma(1) .GT. zsigma(KB) )THEN
                 DO KKK=1,KB
                    UK(KKK)=zsigma(KB-kkk+1)
                 ENDDO
                 DO KKK=1,KB
                   zsigma(kkk)=UK(kkk)
                   UK(KKK)=s(KB-KKK+1,NJ,N)
                 ENDDO
              ELSE     
                 DO KKK=1,KB
                   UK(KKK)=s(KKK,NJ,N)
                 ENDDO
              ENDIF
              IF (sdepth .LT. zsigma(1) )THEN
                UTMP=UK(1)
                GOTO 25
              ELSEIF(sdepth .GT. zsigma(KB) )THEN
                UTMP=UK(KB-1)
                GOTO 25
              ENDIF
!              IF (KB .LT. 5)THEN
                DO K9=1,KB-1
                  IF( (sdepth .GE. zsigma(K9)) .and.
     1                 (sdepth .LT. zsigma(K9+1)) )THEN
                    goto 335
                  ENDIF
                ENDDO
335             X1=zsigma(K9)
                X2=zsigma(K9+1)
                Y1=UK(K9)
                Y2=UK(K9+1)
                call linear(X1,Y1,X2,Y2,sdepth,UTMP)
!              ELSE IF (KB .GE. 5)THEN
!                CALL spline(KB,zsigma,UK,sdepth,UTMP)
!              ENDIF
25          continue
!            if(dayj .GT. 212)
             WRITE(10,100)dayj,IYEAR,ICM,ICD,IHR,IMN,UTMP
           ENDIF
          endif
       ENDDO
       CLOSE(10)
       GOTO 50
999    continue
       CLOSE(9)
100    format(f10.5,I5,4I3,50f10.4)
      stop
      end
      
