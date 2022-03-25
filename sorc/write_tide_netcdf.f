C...  columncatnetcdf.f
C...   Reads in several files and adjusts their
C...   data onto a single time line 
C...   Then outputs a single file with 
C...   year yday h1 h2 h3 h4 h5 ....
C...   Reads input file with  
C...   year month day hour min data
C...     (unlike river version with includes seconds
C  f90 write_tide_netcdf.f  Hydro_netcdfs_station.f  Hydro_netcdfs_grid.f  -o write_tide_netcdf.x \
C  -I/disks/NASWORK/ngofs/oqctools/netcdfsgi/include -L/disks/NASWORK/ngofs/oqctools/netcdfsgi/lib -lnetcdf
CC   run:  write_tide_netcdf.x < columncat.input


      PARAMETER(NARRAY=200000)
      character*120 fileinput
      character*4 fileshort
      integer i,numtimes,ifile,numfiles,iyear
!      real*8 ist_yr, ist_mon, ist_day, ist_hr
!      real*8 lst_yr, lst_mon, lst_day, lst_hr
!      real*8 dt,dtj
!      real*8 jdayfirst,jdaylast
      real time(NARRAY),obs(NARRAY),jday(NARRAY),A(NARRAY,20)
      real obsint(NARRAY),time_int(NARRAY)
!      real*8 julian, dy ,year,day,month,hour,min,sec,yday   
!      dimension time_p(NARRAY),obs_p(NARRAY),day_p(NARRAY)
!     1 ,obsint_p(NARRAY)     
      integer nobs
      real lat,latm,lon,lonm
      character*1 lata,lona
      Character*80 scratchdir
      
C Netcdf arrays
c station variables:
c over dimensioning istation with nnvs=1 works,  might not if nnvs!=1
      parameter (istation=30,nnvs=1)
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
     & ,'NWLON_Tide_Prediction_Data','MSL'
     & ,'cgi-bin/co-ops_qry_direct.cgi',
     &  'NOAA/CSDL/MMAP','tide_prediction.sh','co-ops.nos.noaa.gov'/
           
      read(*,*) netcdf_file_s
      write(*,*) 'use netcdf filename=',netcdf_file_s
      read(*,*) scratchdir
      write(*,*) 'use scratchdir=',scratchdir
      nscratch=index(scratchdir,' ')-1
      write(*,*) 'nscratch=',nscratch,'>',scratchdir(1:nscratch),'<'
      
      write(*,*) ' ist_yr , ist_mon, ist_day, ist_hr?'
      read(*,*) ist_yr , ist_mon, ist_day, ist_hr
      IYR0=ist_yr
      DAYS=365.
      IF(MOD(IYR0,4) .EQ. 0)DAYS=366.
      CALL CONCTJ(IJD,ist_mon, ist_day,ist_yr)
      dayfirst=(ist_yr-IYR0)*DAYS+IJD+float(ist_hr)/24.
      write(*,*)  ' columncat starttime ', 
     & ist_yr , ist_mon, ist_day, ist_hr, dayfirst

      write(*,*) ' lst_yr , lst_mon, lst_day, lst_hr?'
      read(*,*) lst_yr , lst_mon, lst_day, lst_hr
      CALL CONCTJ(IJD,lst_mon, lst_day,lst_yr)
      daylast=(lst_yr-IYR0)*DAYS+IJD+float(lst_hr)/24.
      write(*,*) ' columncat endtime   ',
     & lst_yr , lst_mon, lst_day, lst_hr, daylast
      write(*,*) ' Output time interval (hr)?'
      read(*,*) dt
      write(*,*)  'dt = ',dt
      write(*,*) ' Input time interval (hr)?'
      read(*,*) dt0
      write(*,*)  'dt0 = ',dt0
      write(*,*) ' Interpolation Method 0:Spline; 1:SVD?'
      read(*,*) method
      write(*,*)  'method = ',method
      dtj = dt/24.
      numtimes = int((daylast-dayfirst)/dtj +1+0.1)
      write(*,*) ' numtimes =',numtimes
      do i = 1,numtimes
         time(i) = dayfirst + dtj*(i-1)
      enddo
      write(*,*)'time(1)=',time(1),'time(N)=',time(numtimes)           
c  Number of input files  netcdf parameter istation fixed
      numfiles=istation

      do ifile=1,numfiles
        write(*,*) ' Read next filename'
        read(*,*,err=1010,end=1010) stationij(ifile),fileshort,
     &   stationnames(ifile)
        write(*,*) '<>',stationij(ifile),'<>',fileshort,'<>',
     &   stationnames(ifile),'<>' 
      fileinput=scratchdir(1:nscratch)//'/'//fileshort//'_tide.6min'
        write(*,*) 'fileinput=',fileinput,'<>'
        read(*,*)lats(ifile),lons(ifile),dirflood
        write(*,*)  lats(ifile),'  ',lons(ifile)     
        open(10,file=fileinput)
!        read(10,*,end=1000)amsl,amllw
        amsl=0.
        amllw=0.0
        do i=1,NARRAY
          read(10,*,end=1000) year,rmonth,day,hour,rmin,obs0
          obs(i)=obs0-(amsl-amllw) 
          sec = 0.0d00
          IYEAR=YEAR
          imon=rmonth
          iday=day
          CALL CONCTJ(IJD,imon, iday,iyear)
          jday(i)=(iyear-IYR0)*DAYS+IJD+hour/24.+rmin/1440.
          if (jday(i) .gt. daylast )goto 1000
        enddo
 1000   nobs=i-1
        if (nobs.eq.0) then
c.   Empty file.  (Ice in Choptank) Fill with -999.0
           nobs = 2
           jday(1) = dayfirst
           jday(nobs) =daylast
           obs(1)=-999.0
           obs(nobs)=-999.0
        endif
      write(*,'('' obsdata(1)    ''f12.4,g12.4,'' m'')')  
     &    jday(1),obs (1)
      write(*,'('' obsdata('',i8,'') ''f12.4,g12.4,'' m'')')                    
     &    nobs,jday(nobs),obs (nobs)
      close(10)
!      call interp1(nobs,jday,obs,numtimes,time,obsint)
!      DO KI=1,nobs
!         write(*,*)'raw=',KI,jday(ki),obs(ki)
!      ENDDO
      write(*,*)'raw data number=',nobs,jday(1),jday(nobs)
      CALL gap_filling(dayfirst,daylast,dt,dt0,method,
     1  jday,obs,time_int,obsint,nobs,num_int)
      IF (num_int .ne. numtimes )then
        write(*,*)'number after interpolation is not 
     1  equal to numtimes',num_int,numtimes   
       stop
!      else
!        do KI=1,numtimes
!         obsint(KI)=obsint_p(KI)
!        ENDDO
      endif  
      do i=1,numtimes
      A(i,ifile)=obsint(i)
      enddo

C loop on files
      enddo
 1010 write(*,*) 'last filenumber read attempt', ifile
      istations=ifile-1
 
c done reading all the files 
c initialize netcdf file
c      netcdf_file_s='nwlonstation.nc'
      ibasedate(1) = ist_yr
      ibasedate(2) = 1
      ibasedate(3) = 1
      ibasedate(4) = 0
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
      write(*,*) 'columncatnetcdf.x ',netcdf_file_s
      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,1,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday,ibasedate,lons,lats,1.0,TOTDEPs,
     & 1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)
      
      write(*,*) ' start output',numtimes, numfiles
      
C output the lines of data
      do i=1,numtimes
        yday=time(i)
      write(22,77) IYR0,yday,(A(i,ifile),ifile=1,istations)
  77  format(i4, f12.6, 22f10.2)
  
      do ifile=1,istations
       zs(ifile) = A(i,ifile)
      enddo
      
       yday1 = yday-1.
       write(*,*) 'call write_netcdf_station 2,  ',yday1
      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,2,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday1,ibasedate,lons,lats,1.0,TOTDEPs,
     & zs,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)

      enddo

      call write_netcdf_Hydro_station(netcdf_file_s,ncidst,3,
     & globalstr,istations,stationnames,stationij,1,nnvs,
     & yday1,ibasedate,lons,lats,1.0,TOTDEPs,
     & zs,-1.,-1.,-1.,-1.,-1.,-1.,-1.)
c     & zs,Us,Vs,Ws,Ts,Ss,wx,wy)

      stop      
      end
      
      
      
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
      END FUNCTION JULIAN
      
c      SUBROUTINE JULIANSUB(jday,yr,month,day,hour)
c      
c      jday = JULIAN(yr,month,day,hour)
c      END
            
      SUBROUTINE GREGORIAN(jday,yr,month,day,hour)
      real*8 jday,yr,month,day,hour,dayoweek,week
       
      a = aint(jday+.5)
      b = a+1537
      c = aint((b-122.1)/365.25)
      d = aint(365.25*c)
      e = aint((b-d)/30.6001)
      day = b-d-aint(30.6001*e)+mod(jday+0.5d00,1.0d00)
      month = e-1-12*aint(e/14)
      yr = c-4715-aint((7+month)/10)
      hour = mod(jday+0.5d00,1.0d00)*24.
      day = aint(day)
      
!      dayoweek = mod(aint(jday+0.5),7)
!      week = aint((jday-2444244.5)/7)
!  2000 leap year crashed.  
      if(month.eq.2. .and. day .eq. 31. ) then
       yr=yr-1.
       day = 29.
      endif
      
      END

      SUBROUTINE interp1(nk,x,y,ni,xi,yi)
c.....  x,y are filled data vectors  
c.....   xi  is vector of new xi,yi
c.....   yi will be filled with interpolated values at xi inside x,y
c.....  x should be monotonic
c.....   It won't bomb otherwise, but don't trust the method!
      real *8 x(*),y(*),xi(*),yi(*)
      integer nk, ni, kf, kl, i
      
c      write(*,*)' interp1 x(1),y(1),x(nk),y(nk)',x(1),y(1),nk,x(nk),y(nk)
      kf = 1
      kl = 2
      do i = 1,ni
      if (xi(i).lt.x(1)) then
c.....  first interpolates are before the data starts
           yi(i) = y(kf)
      else
        do while(xi(i).gt.x(kl) .and. kl.lt.nk)
c.....  Reposition so that x(kf)<= xi(i) < x(kl)
c.....   When xi(i)>x(nk) it runs up till  kl =nk
	  kf = kf+1
	  kl = kl+1
	enddo
	if (xi(i).gt.x(kl)) then
c.....  last interpolates are after the data ends
	    yi(i) = y(kl)
	else
c.....  linear interpolate for the range  x(kf)<= xi(i) < x(kl)  
            yi(i) = y(kf) + (xi(i)-x(kf))*(y(kl)-y(kf))/(x(kl)-x(kf))
	endif 
      endif
      enddo
c      write(*,*)' interp1 xi(1),yi(1),xi(ni),yi(ni)',xi(1),yi(1),ni,xi(ni),yi(ni)
      return
      end
CCC***************************************************************      
      subroutine gap_filling(day_begin,day_end,delta,delt0,method,
     1  time,wl,time_new,wl_new,num,m_new)
CCCC Program gap_filling.f  is used to fill data gaps and form a continuous time series
CCCC if gap <2 hours, using linear interpolated values to fill the data gap, 
CCCC if 2< gap <=6 hours, using Spline or SVD interpolated values to fill the data gap, 
CCCC if gap >6 hours, using -999.0 to fill the data gap
CCCC method=0, using Cubic Spline interpolation method
CCCC method=1, using Singular Value Decomposition method
CC        delta= time interval of gap-filt output time series
CC        delt0= time interval of original input time series in hours 
C   compile    F90 gap_filling.f -o gap_filling.x
C   run        gap_filling.x < gap_filling.ctl 

      parameter(ndata=20,ma=5)
      dimension time(num),wl(num),wl_new(90000),time_new(90000)
     1  ,XTMP(200),YTMP(200),X(ndata),Y(ndata),a(ma),afunc(ma)
      hour_begin=day_begin*24.
      hour_end=day_end*24.0
      write(*,*)'running gap_filling...'
      IF (NUM .LE. 1)STOP
      diff=0.
      DO I=1,num
        time(i)=time(i)*24.0
!        if (i .ge. 2)diff=time(i)-time(i-1)
!        write(*,*) i,time(i),diff
      ENDDO
      write(*,*)'begin and end time=',hour_begin,hour_end 
      write(*,*)'raw data number=',num,time(1),time(num)
cccc calculate maximum gap        
      time0=hour_begin
      Y0=wl(1)
      gp_max=0.0
      gp=1.5*delta     !!! gap > 0.15*60 = 9 minutes
      M_new=0
       IF (time(1) .gt. hour_begin)then
         print *,time(1),hour_begin
         K=0
30       T1=hour_begin+K*delta
         IF (T1 .LT. time(1) )then
          m_new=M_new+1
          time_new(M_new)=T1
          wl_new(m_new)=-999.0
          K=K+1
	  goto 30
         ELSE
          TIME0=TIME(1)
          Y0=wl(1)
          goto 40
         ENDIF     
      ENDIF
40    CONTINUE
      do i=1,NUM
         if (time(i) .gt.hour_end)goto 255 
         gap=time(i)-time0
         if( (I .LE. 10) .or. (I .GE. NUM-10) ) then
             if (gap .le. gp)then
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             else if((gap.ge.gp) .and. (gap.le. 3.0))then
               WRITE(*,*)time0,time(i),gap
               TIMEb=TIME0
               TIMEE=time(i)
               NN=NINT(gap/delta)-1
               YB=Y0
               YE=wl(i)
               write(*,*)'call Linear interpolation..'
               CALL linear(TIMEb,yB,TIMEE,yE,delta,XTMP,YTMP,NN)
               DO KK=1,NN
                m_new=M_new+1
                time_new(M_new)=XTMP(kk)
                wl_new(m_new)=YTMP(KK)
               ENDDO
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             ELSE IF (gap .gt. 3.0)then
               TIMEb=TIME0
               TIMEE=time(i)
               NN=NINT(gap/delta)-1
               write(*,*)TIMEb,TIMEE,NN
               DO KK=1,NN
                 m_new=M_new+1
                 time_new(M_new)=TIMEB+kk*delta
                 wl_new(m_new)=-999.0
               ENDDO
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             ENDIF
         ENDIF         
         if( (I .gt. 10) .and. (I .lt. NUM-10) )then
             if (gap .le. gp)then
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             else if((gap.gt.gp) .and. (gap.le. 2.0))then
               WRITE(*,*)'gap= ',time0,time(i),gap
               TIMEb=TIME0
               TIMEE=time(i)
               NN=NINT(gap/delta)-1
               YB=Y0
               YE=wl(i)
               write(*,*)'call Linear interpolation..'
               CALL linear(TIMEb,yB,TIMEE,yE,delta,XTMP,YTMP,NN)
               DO KK=1,NN
                m_new=M_new+1
                time_new(M_new)=XTMP(kk)
                wl_new(m_new)=YTMP(KK)
               ENDDO
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             else if( (gap.gt.2.0) .and. (gap.le. 6.0))then
               WRITE(*,*)'gap= ',time0,time(i),gap
               NFIT=ndata
               IF (delt0 .ge. 0.5)NFIT=6
               DO J=1,NFIT
                 X(J)=time(I-NFIT/2+j)
                 Y(J)=wl(I-NFIT/2+j)
               ENDDO
               TIMEb=TIME0
               TIMEE=time(i)
               NN=NINT(gap/delta)-1
               write(*,*)TIMEb,TIMEE,TIMEE-TIMEB,NN
               IF ((x(NFIT)-x(1)) .GT. 15*delt0)THEN
                 DO KK=1,NN
                   XTMP(kk)=TIMEb+KK*delta
                   YTMP(KK)=-999.0
                 ENDDO
                 goto 444 
               ENDIF 
               IF (method .eq. 0)then      
                write(6,*)'call Spline, method=',method    
                CALL spline(NFIT,x,y,TIMEB,TIMEE,delta,XTMP,YTMP,NN)
               ELSE  
                write(6,*)'call SVD, method=',method    
                CALL SVD(NFIT,ma,x,y,TIMEB,TIMEE,delta,XTMP,YTMP,NN)
               ENDIF
444            CONTINUE
               DO KK=1,NN
                m_new=M_new+1
                time_new(M_new)=XTMP(kk)
                wl_new(m_new)=YTMP(KK)
               ENDDO
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             ELSE IF (gap .gt. 6.0)then
               TIMEb=TIME0
               TIMEE=time(i)
               NN=NINT(gap/delta)-1
               write(*,*)TIMEb,TIMEE,NN
               DO KK=1,NN
                 m_new=M_new+1
                 time_new(M_new)=TIMEB+kk*delta
                 wl_new(m_new)=-999.0
               ENDDO
               m_new=M_new+1
               time_new(M_new)=time(i)
               wl_new(m_new)=wl(i)
             ENDIF 
         ENDIF
         TIME0=TIME(i)
         Y0=wl(i)
      ENDDO
255   continue
      write(6,*)'time_new(M_new)=',M_new,time_new(M_new)
      IF (time_new(M_new) .LT. hour_end)THEN
         gap=hour_end-time_new(M_new)
         TIMEB=time_new(M_new)
         TIMEE=hour_end
         NN=NINT(gap/delta)
         write(*,*)TIMEb,TIMEE,gap,NN
         if( gap.le. 2.0)then
             NFIT=ndata
             IF (delt0 .ge. 0.5)NFIT=6
             DO J=1,NFIT
                X(J)=time_new(M_new-NFIT+j)
                Y(J)=wl_new(M_new-NFIT+j)
             ENDDO
!             IF ((x(NFIT)-x(1)) .GT. 15*delt0)THEN
!               write(6,*)'X(N)-X(1) > ',15*delt0
!               DO KK=1,NN
!                 XTMP(kk)=TIMEB+KK*delta
!                 YTMP(KK)=-999.0
!               ENDDO
!               goto 555 
!             ENDIF 
             IF (method .eq. 0)then      
               CALL spline(NFIT,x,y,TIMEB,TIMEE,delta,XTMP,YTMP,NN)
             ELSE  
               write(6,*)'call SVD, method=',method    
               CALL SVD(NFIT,ma,x,y,TIMEB,TIMEE,delta,XTMP,YTMP,NN)
             ENDIF
555          CONTINUE
             DO KK=1,NN
               m_new=M_new+1
               time_new(M_new)=XTMP(kk)
               wl_new(m_new)=YTMP(KK)
             ENDDO
          else
             DO KK=1,NN
               m_new=M_new+1
               time_new(M_new)=TIMEB+KK*delta
               wl_new(m_new)=-999.0
             ENDDO
          ENDIF
      ENDIF
!      DO K=1,M_NEW
!        write(12,'(2F12.4)')time_new(K)/24.,wl_new(K)
!      ENDDO        
  999 return
      end

       subroutine linear(TIMEB,YB,TIMEE,YE,delta,X_NEW,Y_NEW,NNEW)
        DIMENSION X_NEW(NNEW),Y_NEW(NNEW)
        slope=(YE-YB)/(TIMEE-TIMEB)
        DO J=1,NNEW
          XA=TIMEB+J*delta
          YA=YB+slope*(XA-TIMEB)
          X_NEW(J)=XA
          Y_NEW(J)=YA
        ENDDO
        RETURN
        END
       SUBROUTINE spline(num,x,y,day_begin,day_end,delta,
     1   X_NEW,Y_NEW,NNEW)
	DIMENSION X(NUM),Y(NUM),Y2(NUM),U(NUM)
        DIMENSION X_NEW(NNEW),Y_NEW(NNEW)
        YP1=(Y(2)-Y(1))/(X(2)-X(1) )
        YPN=(Y(NUM)-Y(NUM-1))/(X(NUM)-X(NUM-1) ) 
!!!IBC=1 FOR NATURAL CUBIC-SPLINE BOUNDARY CONDITION , WHICH HAS ZERO SECOND DERIVATIVE ON ONE OR BOTH END POINTS
!  IBC=2 SET EITHER OF Y2(1) AND Y2(N) TO VALUES CALCULATED FROM EQUATION (3.3.5) SO AS TO MAKE THE FIRST DERIVATIVE OF THE INTERPOLATING
                !  FUNCTION HAVE A SPECIFIED VALUE ON EITHER OR BOTH END POINTS.
      IBC=1
!	IF (YP1 .GT. 1.0E30) THEN
	IF (IBC .EQ. 1)THEN
	  Y2(1)=0.0
	  U(1)=0.0
	 ELSE
	  Y2(1)=-0.5
	  U(1)=(3.0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
	ENDIF
	DO I=2,NUM-1
	 SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
	 P=SIG*Y2(I-1)+2.0
	 Y2(I)=(SIG-1.0)/P
	 U(I)=(6.0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1   /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
	END DO
!	IF (YPN .GT. 1.0E30)THEN
	IF (IBC .EQ. 1) THEN
	  QN=0.0
	  UN=0.0
	 ELSE
	  QN=0.5
	  UN=(3.0/(X(NUM)-X(NUM-1)))*(YPN-(Y(NUM)-Y(NUM-1))
     1       /(X(NUM)-X(NUM-1)))
	END IF
	Y2(NUM)=(UN-QN*U(NUM-1))/(QN*Y2(NUM-1)+1.0)
	DO K=NUM-1,1,-1
	  Y2(K)=Y2(K)*Y2(K+1)+U(K)
	END DO
        DO J=1,NNEW
          XA=day_begin+J*delta
      	  print *, 'interpolating with Spline at time ',xa
	  KLO=1
 	  KHI=NUM
1	  IF (KHI-KLO .GT. 1)THEN
	     K=(KHI+KLO)/2
	     IF (X(K) .GT. XA)THEN
	       KHI=K
             ELSE
               KLO=K
	     END IF
             GOTO 1
	  END IF
  	  H=X(KHI)-X(KLO)
	  IF (ABS(H) .LT. 1.0E-5) THEN
            PRINT *, 'Bad XA input in SPLINE'
            STOP
          END IF
	  A=(X(KHI)-XA)/H
	  B=(XA-X(KLO))/H
	  YA=A*Y(KLO)+B*Y(KHI)+((A**3-A)*Y2(KLO)+(B**3-B)*Y2(KHI))*
     1     (H**2)/6.0 
          X_new(J)=XA
          Y_New(J)=YA
        ENDDO 
!         DO I=1,NN
!           write(*,*)X_new(I)/1440.,Y_NEw(I)
!	 ENDDO
!        write(12,121) ttt,ityr(nt),itmm(nt),itdd(nt),ithh(nt),
!     1           itmn(nt),wltot,tide(nt),subwl
!        write(13,121) ttt,ityr(nt),itmm(nt),itdd(nt),ithh(nt),
!     1           itmn(nt),wltot_1,tide(nt),subwl
!        write(120,*)xa,wltot
  121   format(f9.5,i5,4i3,3f9.4)
99      CONTINUE
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccccc
c              ndata = number of data points
c         x(1:ndata) = data times
c         y(1:ndata) = data values
c       sig(1:ndata) = std. dev. of data
c                  a = coefficients of fitting function
c                      y = a(i)*afunc(i)
c                 ma = number of coefficients
c       u(1:mp,1:np) = work array
c       v(1:np,1:np) = work array
c            w(1:np) = work array
c                 mp = array dimension. ge ndata
c                 np = array dimension. ge ma
c              chisq = chi squared value

c     subroutine funcs(x,afuncs,ma,nmax))
c     SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
c     SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
c     SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
c     SUBROUTINE svdvar(v,ma,np,w,cvm,ncvm)
c     FUNCTION pythag(a,b)
      subroutine SVD(ndata,ma,x,y,TIMEB,TIMEE,delta,XTMP,YTMP,NNEW)
!      parameter (ndata=20,ma=5,mp=ndata,np=ma)
!      parameter (ma=5)
      real x(ndata),y(ndata),sig(ndata)
      real a(ma),u(ndata,ma),v(ma,ma),w(ma)
      dimension afunc(ma),XTMP(NNEW),YTMP(NNEW)
       mp=ndata
       np=ma 
        do n=1,ndata
          x(n)=x(n)-TIMEB
          sig(n)=1.0
        enddo
c     call svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
        call svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq)
!        write(*,*) 'the coefficients'
      do n=1,ma
        write(6,*)n,a(n)
      enddo
c        filling the gap from timeb to timee using SVD predictions 
        do n=1,NNEW
          XA=N*delta
          call funcs(xa,afunc,ma)
          sum=0.
          do 15 j=1,ma
            sum=sum+a(j)*afunc(j)
15        continue
          XTMP(n)=XA+TIMEB
          YTMP(n)=sum
!         write(*,*)x(n)+timeb,sum
        enddo

        RETURN    
        end
c------------------------------------------------------------------------
      subroutine funcs(x,afunc,ma)
      dimension afunc(ma)
      do m=1,ma
        afunc(m)=x**(m-1)
      enddo
      return
      end
c------------------------------------------------------------------------
c     SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
      SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq)
      INTEGER ma,mp,ndata,np,NMAX,MMAX
      REAL chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np),x(ndata),
     *y(ndata),TOL
c     EXTERNAL funcs
      PARAMETER (TOL=1.e-5)
!      PARAMETER (NMAX=1000,MMAX=50,TOL=1.e-5)
CU    USES svbksb,svdcmp
      INTEGER i,j
      REAL sum,thresh,tmp,wmax,afunc(ma),b(ndata)
!      REAL sum,thresh,tmp,wmax,afunc(MMAX),b(NMAX)
      do 12 i=1,ndata
        call funcs(x(i),afunc,ma)
        tmp=1./sig(i)
        do 11 j=1,ma
          u(i,j)=afunc(j)*tmp
11      continue
        b(i)=y(i)*tmp
12    continue
      call svdcmp(u,ndata,ma,mp,np,w,v)
      wmax=0.
      do 13 j=1,ma
        if(w(j).gt.wmax)wmax=w(j)
13    continue
      thresh=TOL*wmax
      do 14 j=1,ma
        if(w(j).lt.thresh)w(j)=0.
14    continue
      call svbksb(u,w,v,ndata,ma,mp,np,b,a)
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),afunc,ma) ! +++
        sum=0.
        do 15 j=1,ma
          sum=sum+a(j)*afunc(j)
15      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
16    continue
      return
      END
c----------------------------------------------------------------------------
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL a(mp,np),v(np,np),w(np)
!      PARAMETER (NMAX=mp)
CU    USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(mp),pythag
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
c-------------------------------------------------------------------------
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
CCC
CC   Solves A*X=B for a vector X, where A is specified by the arrays u,w,v as returned by
CC   svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices.  mp 
CC   and np are the physical dimensions of a. b(1:m) is the input right-hand side. x(1:n) is 
CC   the output solution vector.  No input quantities are destroyed, so the routine may be called
CC   sequentially wiwth different b's.

      INTEGER m,mp,n,np,NMAX
      REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
!      PARAMETER (NMAX=MP)
      INTEGER i,j,jj
      REAL s,tmp(N)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
c----------------------------------------------------------------------
      SUBROUTINE svdvar(v,ma,np,w,cvm,ncvm)
      INTEGER ma,ncvm,np,MMAX
      REAL cvm(ncvm,ncvm),v(np,np),w(np)
!      PARAMETER (MMAX=ma)
!      PARAMETER (MMAX=20)
      INTEGER i,j,k
      REAL sum,wti(MA)
      do 11 i=1,ma
        wti(i)=0.
        if(w(i).ne.0.) wti(i)=1./(w(i)*w(i))
11    continue
      do 14 i=1,ma
        do 13 j=1,i
          sum=0.
          do 12 k=1,ma
            sum=sum+v(i,k)*v(j,k)*wti(k)
12        continue
          cvm(i,j)=sum
          cvm(j,i)=sum
13      continue
14    continue
      return
      END
c----------------------------------------------------------------------
      FUNCTION pythag(a,b)
      REAL a,b,pythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END
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
     
