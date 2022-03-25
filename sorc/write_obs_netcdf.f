C...  columncatnetcdf.f
C...   Reads in several files and adjusts their
C...   data onto a single time line 
C  f90 write_obs_netcdf.f Hydro_netcdfs_station.f  Hydro_netcdfs_grid.f  -o write_obs_netcdf.x \
C  -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf
CC   run:  write_obs_netcdf.x < columncat.input
C  lf95 write_obs_netcdf.f Hydro_netcdfs_station.f  Hydro_netcdfs_grid.f  -o write_obs_netcdf.x\
C  -I/disks/NASPUB/usr/Local/Linux/netcdf/include -L/disks/NASPUB/usr/Local/Linux/netcdf/lib -lnetcdf
      PARAMETER(NARRAY=120000,istation=99,nnvs=1)
      character*120 fileinput
      character*4 fileshort,suffix*20
      integer i,numtimes,ifile,numfiles,iyear
      real time(NARRAY),obs(NARRAY),jday(NARRAY),A(NARRAY,istation)
      real obsint(NARRAY),time_int(NARRAY)
      real*8 jday0,jday1,jbase_date,JULIAN,yr,month,day,hour
      integer nobs
      real lat,latm,lon,lonm
      character*1 lata,lona
      Character*80 scratchdir
      
C Netcdf arrays
c station variables:
c over dimensioning istation with nnvs=1 works,  might not if nnvs!=1
!      parameter (istation=99,nnvs=1)
      Character*80 netcdf_file_s
      integer stationij(istation)
      character stationnames(istation)*40
      real lons(istation),lats(istation),TOTDEPs(istation),zs(istation)
      real Us(istation,nnvs),Vs(istation,nnvs),Ws(istation,nnvs), 
     &  Ts(istation,nnvs),Ss(istation,nnvs)
      integer ibasedate(4)
      character globalstr(9)*40
	data globalstr/
     &  'Station_Observation','Surface','NWLON'
     & ,'NWLON_Observation_Data','MSL'
     & ,'cgi-bin/co-ops_qry_direct.cgi',
     &  'NOAA/CSDL/MMAP','getnwlon.sh','co-ops.nos.noaa.gov'/
           
      read(*,'(a80)') netcdf_file_s
      write(*,*) 'use netcdf filename=',netcdf_file_s
      write(*,*) ' ist_yr , ist_mon, ist_day, ist_hr?'
      read(*,*) ist_yr , ist_mon, ist_day, ist_hr
      ibasedate(1) = ist_yr
      ibasedate(2) = 1
      ibasedate(3) = 1
      ibasedate(4) = 0
       yr=ibasedate(1)
       month=ibasedate(2)
       day=ibasedate(3)
       hour=ibasedate(4)
       jbase_date=JULIAN(yr,month,day,hour)
       yr=ist_yr 
       month=ist_mon 
       day=ist_day 
       hour=ist_hr 
       jday0=JULIAN(yr,month,day,hour)
       dayfirst=jday0-jbase_date    !!! base date is Jan. 1:Julian day=1
      write(*,*)  ' columncat starttime ', 
     & ist_yr , ist_mon, ist_day, ist_hr, dayfirst

      write(*,*) ' lst_yr , lst_mon, lst_day, lst_hr?'
      read(*,*) lst_yr , lst_mon, lst_day, lst_hr
       yr=lst_yr 
       month=lst_mon 
       day=lst_day 
       hour=lst_hr 
       jday0=JULIAN(yr,month,day,hour)
       daylast=jday0-jbase_date+1.
       write(*,*) ' columncat endtime   ',
     & lst_yr , lst_mon, lst_day, lst_hr, daylast
      write(*,*) ' Input time interval (hr)?'
      read(*,*) dt
      write(*,*)  'dt = ',dt
      dtj = dt/24.
      numtimes = int((daylast-dayfirst)/dtj +1+0.1)
      write(*,*) ' numtimes =',numtimes
      do i = 1,numtimes
         time(i) = dayfirst + dtj*(i-1)
      enddo
      write(*,*)'time(1)=',time(1),'time(N)=',time(numtimes)           
       ifile=0
123    read(*,'(a120)',end=1010)fileinput
        write(*,*) 'fileinput=',fileinput
       call ncrght(fileinput,nct)
!        read(*,*)lats(ifile),lons(ifile),dirflood
!        write(*,*)  lats(ifile),'  ',lons(ifile)     
        open(10,file=fileinput(1:nct))
        do i=1,NARRAY
          read(10,*,end=1000)tday, year,rmonth,rday,rhour,rmin,obs(i)
          yr=year 
          month=rmonth 
          day=rday 
          hour=rhour 
          jday0=JULIAN(yr,month,day,hour)
          jday(i)=jday0+rmin/1440.-jbase_date
          if (jday(i) .gt. daylast )goto 1000
        enddo
 1000   nobs=i-1
        ifile=ifile+1
        do i=1,numtimes
           A(i,ifile)=-999.0
           do j=1,nobs
             if (abs(time(i)-jday(j)) .le. 0.0001)then
                A(i,ifile)=obs(j)
             endif
           enddo
        enddo
        goto     123        
!        write(*,*)'nobs=',nobs
!      enddo
 1010 write(*,*) 'last filenumber read attempt', ifile
      istations=ifile
!      do i=1,numtimes
!      time(i)=time(i)-1.0  !!! base date is Jan. 1:Julian day=1.0
!      enddo
       
c done reading all the files 
c initialize netcdf file
c      netcdf_file_s='nwlonstation.nc'
      do i=1,istations
       TOTDEPs(i) = 0.0
       zs(i) = 0.0
       Us(i,1) = 0.0
       Vs(i,1) = 0.0
       Ws(i,1) = 0.0
       Ts(i,1) = 0.0
       Ss(i,1) = 0.0
      enddo
      
      yday=time(1)
      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,1,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday,ibasedate,lons,lats,1.0,TOTDEPs,
     & 1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)
      
      write(*,*) ' start output',numtimes
      
C output the lines of data
      do i=1,numtimes
        yday=time(i)
!      write(22,77) IYR0,yday,(A(i,ifile),ifile=1,istations)
  77  format(i4, f12.6, 22f10.2)
  
      do ifile=1,istations
       zs(ifile) = A(i,ifile)
      enddo
      
       
!       write(*,*) 'call write_netcdf_station 2,  ',yday
      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,2,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday,ibasedate,lons,lats,1.0,TOTDEPs,
     & zs,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)

      enddo

      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,3,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday,ibasedate,lons,lats,1.0,TOTDEPs,
     & zs,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)

      stop      
      end
c----------------------------------------------------------------------
      SUBROUTINE CONJTC (IJD0,ICM,ICD,IYR) 
      DIMENSION   IDTBLE(12,2),ILTBLE(12,2) 
C 
      DATA ((IDTBLE(I,J),J=1,2),I=1,12)/1,31,32,59,60,90,91,120,121,151,
     1                                  152,181,182,212,213,243,244,273,
     2                                  274,304,305,334,335,365/
      DATA ((ILTBLE(I,J),J=1,2),I=1,12)/1,31,32,60,61,91,92,121,122,152,
     1                                  153,182,183,213,214,244,245,274,
     2                                  275,305,306,335,336,366/
      IJD=IJD0
      IDAYS=365 
      IF(MOD(IYR,4) .EQ. 0) IDAYS=366
      IF(IJD .GT. IDAYS) THEN
       IJD=IJD - IDAYS
       IYR=IYR+1
       END IF 
C 
C**** TEST FOR LEAP YEAR
C 
      ISW = 1 
      IF (MOD(IYR,4) .EQ. 0)   ISW = 2
C 
      L = IJD/30
      IF (L.EQ.0)   GO TO 4 
      GO TO (6,7)   ISW 
    6 IF (IJD.GT.IDTBLE(L,2))   GO TO 4 
      GO TO 8 
    7 IF (IJD.GT.ILTBLE(L,2))   GO TO 4 
    8 ICM = L 
      GO TO 5 
    4 ICM = L + 1 
    5 GO TO (9,10)   ISW
    9 ICD = IJD - IDTBLE(ICM,1) + 1 
      RETURN
   10 ICD = IJD - ILTBLE(ICM,1) + 1 
      RETURN
      END 
      SUBROUTINE SEC2DAY(TIME0,IYR,IJD,ICM,ICD,IHR,IMN,ISEC)
      TIME=TIME0
      ITDAY=INT(TIME/86400.)
      IDAYS=365 
      IF(MOD(IYR,4) .EQ. 0) IDAYS=366
      IF(ITDAY .GT. IDAYS) THEN
       ITDAY=ITDAY -IDAYS
       TIME=TIME-IDAYS*86400
       IYR=IYR+1
       END IF 
      IJD=ITDAY
      IHR=(TIME-IJD*86400)/3600
      IMN=(TIME-IJD*86400-IHR*3600)/60
      ISEC=(TIME-IJD*86400-IHR*3600-IMN*60)
      IF (ISEC .EQ. 60)THEN
         ISEC=0
         IMN=IMN+1
      ENDIF
      IF (IMN .EQ. 60)THEN
         IMN=0
         IHR=IHR+1
      ENDIF
      IF (IHR .EQ. 24)THEN
         IHR=0
         IJD=IJD+1
      ENDIF
      CALL CONJTC(IJD,ICM,ICD,IYR)
      RETURN
      END
      SUBROUTINE CONCTJ (IJD,IMON,IDAY,IYR) 
C 
C**** THIS SUBROUTINE CONVERTS CALENDER TO JULIAN DAY (IJD) 
C 
      DIMENSION   IDTBLE(12),ILTBLE(12) 
C 
      DATA (IDTBLE(I),I=1,12)/1,32,60,91,121,152,182,213,244, 
     1                        274,305,335/
      DATA (ILTBLE(I),I=1,12)/1,32,61,92,122,153,183,214,245, 
     1                        275,306,336/
C 
C**** TEST FOR LEAP YEAR
C 
      ISW = 1 
      IF (MOD(IYR,4).EQ.0)   ISW = 2
C 
      GO TO (9,10)   ISW
    9 IJD = IDTBLE(IMON) + IDAY - 1 
      RETURN
   10 IJD = ILTBLE(IMON) + IDAY - 1 
      RETURN
      END 
     
C**********************************************************************************************
      SUBROUTINE ncrght (line,nc)
c -
c - Returns the last non blank character position in a string
c -
      CHARACTER*(*) line
      ilim = LEN (line)
      DO 100 i = 1,ilim
        IF(line(i:i) .ne. ' ') THEN
          nc = i
        END IF
  100 CONTINUE
      RETURN
      END
      FUNCTION JULIAN(yr,month,day,hour)
      real*8 jday,yr,month,day,hour
      real*8 Y, m ,julian
c      integer yr

      Y = yr
      if (yr.lt.100. .and. yr.gt.50.) then
        Y = yr+1900.
      elseif (yr.lt.100. .and. yr.le.50.) then
        Y = yr+2000.
      endif

      if (month.le. 2.) then
        Y = Y-1
        m = month +12.
      else
        Y = Y
        m = month
      endif

      JULIAN = aint(365.25*Y) + aint(30.6001*(m+1))
     &        + day + hour/24. + 1720981.50

C  1980 epic....
c      JULIAN = aint(365.25*(Y-1980.)) + aint(30.6001*(m+1))
c     &        + day + hour/24.
c      write(*,'(''Julian ='',f12.4,4f8.2)') julian, yr,month,day,hour
c      END FUNCTION JULIAN
       return
       END
