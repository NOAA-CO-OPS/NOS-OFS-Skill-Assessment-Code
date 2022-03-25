ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Hydro_netcdfs_station.f  Writes standardized NetCDF files for Hydro-models.
c    Copyright (C) 2003,  Thomas F. Gross 
c
c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program; if not, write to the Free Software
c    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c  
c  writes the station oriented netcdf file with reals
c   
c      call write_netcdf_Hydro_station(netcdf_file,ncid,imode,
c     & globalstr,istation,stationnames,stationij,meshdim,l,
c     & time,ibasedate,lon,lat,sigma,depth,zeta,u,v,w,temp,salt)
c   
c
c Inputs:
c netcdf_file    char*80  filename for the netcdf output
c ncid                    netcdf id; generated on initialization 
c imode  1 for initialization, 2 for writing, 3 for closing file
c globalstr Global Attributes.  Set in a data statement like
c        data globalstr/
c     &  'grid_type','z_type','model'
c     & ,'title','comment','source',
c     &  'institution','history','references'/
c 
c istation  number of output stations
c stationnames  char stationnames(istation)*20   Ascii station labels
c stationij(istation,meshdim)   indices of main mesh of the stations
c meshdim   dimension of the main mesh  2 for u(i,j),  1 for fem u(inode)
c                       possibly 3 for three surrounding nodes of fem
c l     dimension of vertical outputs. may be =1
c time    time in days
c ibasedate(4)   iyear, imonth, iday, ihour of base date (time = 0)
c lon(istation) ,lat(istation)    longitude, latitude of stations
c sigma(l)  sigma values for vertical outputs
c zeta(istation)  sea surface displacement
c u(istation,l),v(istation,l),w(istation,l) Velocities
c temp(istation,l),salt(istation,l)         Temperature, Salinity 
c wx(istation), wy(istation)   wind velocities (blowing toward)
c   
C optional variable writing:
c                  Upon initialization set a writing variable =1
C                  If the variable is negative then
c                  no variable is created or written later written.
c             Only options are:  zeta,u,v,w,temp,salt,wx,wy
C   Example first call
c      call write_netcdf_Hydro_station(netcdf_file,ncid,1,
c     & globalstr,istation,stationnames,stationij,meshdim,l,
c     & time,ibasedate,lon,lat,sigma,depth,1.,1.,1.,-1.,-1.,-1.)
c  This will only create and write variables zeta,u,v   
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Hydro_netcdfs_station_options.f 
C  july 23,2002    include globalstr variable
c  july 24,2002    optional variables.  
c                  Upon initialization set a writing variable =1
C                  If the variable is negative then
c                  no variable is created or written later written.
c  july 25, 2002   Add wx,wy for wind components
c  Nov 11, 2002    Changes to bring into COARDS compliance
c                     time:units, sigma:positive, 
c                     global:conventions        :missing_value
c                     Allow Spaces in Strings
c  dec 27,2002  addition of nc_sync to save files even when model crashes
c
c  nov 26, 2003  addition of standard_name for COARDS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_netcdf_Hydro_station(netcdf_file,ncid,imode,
     & globalstr,istation,stationnames,stationij,meshdim,l,
     & time,ibasedate,lon,lat,sigma,depth,
     & zeta,u,v,w,temp,salt,wx,wy)

       implicit none

       include 'netcdf.inc'
c input variables
      integer ncid,imode,istation,l,meshdim
      integer ibasedate(4),stationij(istation,meshdim)
      real time,lon(istation),lat(istation),sigma(l),depth(istation),
     & zeta(istation),u(istation,l),v(istation,l),w(istation,l),
     & temp(istation,l),salt(istation,l),wx(istation),wy(istation)
      character stationnames(istation)*40
      character globalstr(9)*40
      Character netcdf_file*80

C Netcdf internal variables
      integer num, itime
      integer iret,intval(4),CORNER(4),COUNT(4)
      integer station_dim, sigma_dim, time_dim, mesh_dim
      integer  time_id, stat_id, lon_id, lat_id, sigma_id
      integer  depth_id, zeta_id,    sname_id, char_dim
      integer  u_id, v_id, w_id, salt_id, temp_id,wx_id,wy_id
      logical lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      
c date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
      
      save station_dim, sigma_dim, time_dim
      save  time_id, stat_id, lon_id, lat_id, sigma_id
      save  depth_id, zeta_id, sname_id, wx_id,wy_id
      save  u_id, v_id, w_id, salt_id, temp_id
      save lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      
      
      if (imode.eq. 1) then
      write(*,*) 'write_netcdf_Hydro_station imode =',imode
      
C Set optional variable flags
      lz = .TRUE.
      lu = .TRUE.
      lv = .TRUE.
      lw = .TRUE.
      lsalt = .TRUE.
      ltemp = .TRUE.
      lwx = .TRUE.
      lwy = .TRUE.
      if (zeta(1).le.0.0) lz=.FALSE.
      if (u(1,1).le.0.0) lu=.FALSE.
      if (v(1,1).le.0.0) lv=.FALSE.
      if (w(1,1).le.0.0) lw=.FALSE.
      if (salt(1,1).le.0.0) lsalt=.FALSE.
      if (temp(1,1).le.0.0) ltemp=.FALSE.
      if (wx(1).le.0.0) lwx=.FALSE.
      if (wy(1).le.0.0) lwy=.FALSE.
            
C Initialize
c         open file
      iret = nf_create(netcdf_file, NF_CLOBBER, ncid)
      call check_err(iret)
c define dimensions
      iret = nf_def_dim(ncid, 'station', istation, station_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'sigma', l, sigma_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'meshdim', meshdim, mesh_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'charlength', 40, char_dim)
      call check_err(iret)

c	 define variables
C time
      iret = nf_def_var(ncid, 'time', NF_REAL, 1, time_dim, time_id)
      iret = nf_put_att_text(ncid, time_id, 'long_name', 4, 'Time')
      call check_err(iret)
      iret = nf_put_att_text(ncid, time_id, 'units', 4, 'days')
      call check_err(iret)
      iret = nf_put_att_int(ncid, time_id, 'base_date', NF_INT, 4, 
     & ibasedate)
      write(date_string,70) ibasedate(1),ibasedate(2),ibasedate(3)
     & ,ibasedate(4)
 70   format('days since ',I4,'-',I2.2,'-',i2.2,' ',i2,':00:00'
     &  ,' 00:00')
      iret = nf_put_att_text(ncid, time_id, 'units'
     &  ,len_trim(date_string), date_string)
      call check_err(iret)
       iret = nf_put_att_text(ncid, time_id, 'standard_name', 4, 
     & 'time')

C stationij
      intval(1) = station_dim
      intval(2) = mesh_dim
      iret = nf_def_var(ncid, 'stationij', NF_INT, 2,intval, stat_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, stat_id, 'long_name',9,'StationIJ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, stat_id, 'units', 7, 'indices')
      call check_err(iret)

C stationnames
      intval(2) = station_dim
      intval(1) = char_dim
      iret = nf_def_var(ncid, 'stationnames',NF_CHAR,2,intval,sname_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, sname_id, 'long_name'
     &  ,13,'Station Names')
      call check_err(iret)
      iret = nf_put_att_text(ncid, sname_id, 'units', 7, 'char*40')
      call check_err(iret)
C
C lon
      intval(1) = station_dim
      iret = nf_def_var(ncid, 'lon', NF_REAL, 1,intval, lon_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9, 'Longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lon_id, 'standard_name', 9, 
     & 'longitude')
C lat
      iret = nf_def_var(ncid, 'lat', NF_REAL, 1,intval, lat_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8, 'Latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lat_id, 'standard_name', 8, 
     & 'latitude')
C depth
      iret = nf_def_var(ncid, 'depth', NF_REAL, 1,intval, depth_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid,depth_id,'long_name',11,'Bathymetry')
      call check_err(iret)
      iret = nf_put_att_text(ncid, depth_id,'units',6,'meters')
      iret = nf_put_att_text(ncid, depth_id,'positive',4,'down')
      call check_err(iret)
       iret = nf_put_att_text(ncid, depth_id, 'standard_name', 5, 
     & 'depth')
C sigma
      intval(1) = sigma_dim
      iret = nf_def_var(ncid, 'sigma', NF_REAL, 1,intval, sigma_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, sigma_id, 'long_name', 44, 
     & 'Sigma Stretched Vertical Coordinate at Nodes')
      call check_err(iret)
      iret = nf_put_att_text(ncid, sigma_id,'units',11,'sigma_level')
      iret = nf_put_att_text(ncid, sigma_id,'positive',4,'down')
      call check_err(iret)
       iret = nf_put_att_text(ncid, sigma_id, 'standard_name', 22, 
     & 'ocean_sigma_coordinate')
       iret = nf_put_att_text(ncid, sigma_id, 'formula_terms', 35, 
     & 'sigma: sigma eta: zeta depth: depth')
C zeta
      intval(2) = time_dim
      intval(1) = station_dim
      if (lz) then
      iret = nf_def_var(ncid, 'zeta', NF_REAL, 2,intval, zeta_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'long_name', 23, 
     & 'Water Surface Elevation')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 6, 'meters')
      call check_err(iret)
      iret = nf_put_att_real(ncid, zeta_id, 
     & 'missing_value', NF_REAL, 1,-999.0)
      iret = nf_put_att_text(ncid, zeta_id,'positive',2,'up')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 21, 
     & 'sea_surface_elevation')
      endif
c wind velocity East
      if (lwx) then
      iret = nf_def_var(ncid, 'air_u', NF_REAL, 2,intval, wx_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, wx_id, 'long_name', 21, 
     & 'Eastward Air Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, wx_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, wx_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, wx_id, 'standard_name', 13, 
     & 'eastward_wind')
      endif
c wind velocity North
      if (lwy) then
      iret = nf_def_var(ncid, 'air_v', NF_REAL, 2,intval, wy_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, wy_id, 'long_name', 22, 
     & 'Northward Air Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, wy_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, wy_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, wy_id, 'standard_name', 14, 
     & 'northward_wind')
      endif
C u
      intval(3) = time_dim
      intval(2) = sigma_dim
      intval(1) = station_dim
      if (lu) then
      iret = nf_def_var(ncid, 'u', NF_REAL, 3,intval, u_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'long_name', 23, 
     &'Eastward Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, u_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'standard_name', 27, 
     & 'eastward_sea_water_velocity')
      endif
C v
      if (lv) then
      iret = nf_def_var(ncid, 'v', NF_REAL, 3,intval, v_id)
      iret = nf_put_att_text(ncid, v_id, 'long_name', 24, 
     &'Northward Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, v_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'standard_name', 28, 
     & 'northward_sea_water_velocity')
      endif
C w
      if (lw) then
      iret = nf_def_var(ncid, 'w', NF_REAL, 3,intval, w_id)
      iret = nf_put_att_text(ncid, w_id, 'long_name', 23, 
     & 'Vertical Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, w_id, 'units', 3, 'm/s')
      call check_err(iret)
      iret = nf_put_att_real(ncid, w_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, w_id, 'standard_name', 25, 
     & 'upward_sea_water_velocity')
      endif
C temp
      if (ltemp) then
      iret = nf_def_var(ncid, 'temp', NF_REAL, 3,intval, temp_id)
      iret = nf_put_att_text(ncid, temp_id, 'long_name', 11, 
     & 'Temperature')
      call check_err(iret)
      iret = nf_put_att_text(ncid, temp_id, 'units', 7, 'Celsius')
      call check_err(iret)
      iret = nf_put_att_real(ncid, temp_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, temp_id, 'standard_name', 21, 
     & 'sea_water_temperature')
      endif
C salt
      if (lsalt) then
      iret = nf_def_var(ncid, 'salt', NF_REAL, 3,intval, salt_id)
      iret = nf_put_att_text(ncid, salt_id, 'long_name', 8, 
     & 'Salinity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, salt_id, 'units', 3, 'ppt')
      call check_err(iret)
      iret = nf_put_att_real(ncid, salt_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      call check_err(iret)
      iret = nf_put_att_text(ncid, salt_id, 'standard_name', 18, 
     & 'sea_water_salinity')
      endif
      write(*,*) 'end variables '
      
C Global Attributes
      iret = nf_put_att_text(ncid, NF_GLOBAL,'file_type'
     & ,7,'Station')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'Conventions'
     & ,6,'COARDS')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'grid_type'
     & ,len_trim(globalstr(1)),globalstr(1))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'z_type'
     & ,len_trim(globalstr(2)),globalstr(2))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'model'
     & ,len_trim(globalstr(3)),globalstr(3))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'title'
     & ,len_trim(globalstr(4)),globalstr(4))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'comment'
     & ,len_trim(globalstr(5)),globalstr(5))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'source'
     & ,len_trim(globalstr(6)),globalstr(6))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'institution'
     & ,len_trim(globalstr(7)),globalstr(7))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'history'
     & ,len_trim(globalstr(8)),globalstr(8))
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'references'
     & ,len_trim(globalstr(9)),globalstr(9))
      call check_err(iret)

      call date_and_time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
      call check_err(iret)


      write(*,*) 'end globals '
            

c         END DEFINATIONS 
      iret = nf_enddef(ncid)
      call check_err(iret)

c	 write lon,lat,mask,depth,sigma
      CORNER(1) = 1
      COUNT(1)=istation
      iret=nf_put_vara_real(ncid,lon_id,CORNER,COUNT,lon)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,lat_id,CORNER,COUNT,lat)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,depth_id,CORNER,COUNT,depth)
      call check_err(iret)
      COUNT(1)=l
      iret=nf_put_vara_real(ncid,sigma_id,CORNER,COUNT,sigma)
      call check_err(iret)
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=istation
      COUNT(2)=meshdim
      iret=nf_put_vara_int(ncid,stat_id,CORNER,COUNT,stationij)
      call check_err(iret)
      write(*,*) 'stationnames:',stationnames,':stationnames'
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=40
      COUNT(2)=istation 
      iret=nf_put_vara_text(ncid,sname_id,CORNER,COUNT,stationnames)
      call check_err(iret)



      write(*,*) 'Initialize station nc done',COUNT(1)


C time index counter
      itime=0
C other indices
      
      elseif (imode.eq.2) then
c         write time,u,v,....
c Inquire of this file what the last itime written was
C  should be equal to the dimension of time
      iret= nf_inq_dimlen(ncid,time_dim,itime)      
      itime=itime+1
!      write(*,*) 'write_netcdf_Hydro_station imode =',imode,'  itime='
!     &  ,itime,'  ',netcdf_file
c scalars
      CORNER(1) = itime
      iret=nf_put_var1_real(ncid,time_id,CORNER,time)
      call check_err(iret)
c single values at istation locations
      CORNER(1) = 1
      CORNER(2) = itime
      COUNT(1)=istation
      COUNT(2)=1
            if (lz) then
       iret=nf_put_vara_real(ncid,zeta_id,CORNER,COUNT,zeta)
       call check_err(iret)
      endif
      if (lwx) then
       iret=nf_put_vara_real(ncid,wx_id,CORNER,COUNT,wx)
       call check_err(iret)
      endif
      if (lwy) then
       iret=nf_put_vara_real(ncid,wy_id,CORNER,COUNT,wy)
       call check_err(iret)
      endif
c profiles of length l at istation locations
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = itime
      COUNT(1)=istation
      COUNT(2)=l
      COUNT(3)=1
      if (lu) then
       iret=nf_put_vara_real(ncid,u_id,CORNER,COUNT,u)
       call check_err(iret)
      endif
      if (lv) then
       iret=nf_put_vara_real(ncid,v_id,CORNER,COUNT,v)
       call check_err(iret)
      endif
      if (lw) then
       iret=nf_put_vara_real(ncid,w_id,CORNER,COUNT,w)
       call check_err(iret)
      endif
      if (ltemp) then
       iret=nf_put_vara_real(ncid,temp_id,CORNER,COUNT,temp)
       call check_err(iret)
      endif
      if (lsalt) then
       iret=nf_put_vara_real(ncid,salt_id,CORNER,COUNT,salt)
       call check_err(iret)
      endif

c  dec 27,2002  addition of sync to save files even when model crashes
       iret=nf_sync(ncid)
       call check_err(iret)

      else
c         close file
      iret = nf_close(ncid)
      call check_err(iret)
      write(*,*) 'imode =',imode,' close ',netcdf_file 
      endif
      
      return
      end

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c      subroutine check_err(iret)
c      integer iret
c      include 'netcdf.inc'
c      if (iret .ne. NF_NOERR) then
c      print *, nf_strerror(iret)
c      stop
c      endif
c      end

