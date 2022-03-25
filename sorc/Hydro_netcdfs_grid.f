ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Hydro_netcdfs_grid.f  Writes standardized NetCDF files for Hydro-models.
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
c  Writes the grid oriented netcdf file with reals
c   
c      subroutine write_netcdf_Hydro(netcdf_file,ncid,imode,
c     &  globalstr,m,n,l,
c     & time,ibasedate,lon,lat,mask,sigma,depth
c    & ,zeta,u,v,w,temp,salt,winde,windn)
c
c And version which scales to integers using prescribed ranges:
c      subroutine write_netcdf_Hydro_scale(netcdf_file,ncid,imode,
c     & globalstr,m,n,l,
c     & time,ibasedate,lon,lat,mask,sigma,depth,
c     & zeta,u,v,w,temp,salt,wx,wy)
c
c
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
c m     dimension of the X coordinate (Longitude)
c n     dimension of the Y coordinate (Latitude)
c l     dimension of vertical outputs. may be =1
c time    time in days
c ibasedate(4)   iyear, imonth, iday, ihour of base date (time = 0)
c lon(m,n) ,lat(m,n)    longitude, latitude of stations
c sigma(l)  sigma values for vertical outputs 0:-1 0 surface, -1 seabed
c zeta(m,n)  sea surface displacement
c u(m,n,l),v(m,n,l),w(m,n,l) Velocities
c temp(m,n,l),salt(m,n,l)         Temperature, Salinity 
c we(m,n),wn(m,n) Wind Velocity Vectors (wind toward)
c   
C optional variable writing:
c                  Upon initialization set a writing variable =1
C                  If the variable is negative then
c                  no variable is created or written later written.
c             Only options are:  zeta,u,v,w,temp,salt,we,wn
C   Example first call
c      subroutine write_netcdf_Hydro(netcdf_file,ncid,imode,
c     &  globalstr,m,n,l,
c     & time,ibasedate,lon,lat,mask,sigma,depth,1.,1.,1.,0.,0.,0.,0.,0.)
c  This will only create and write variables zeta,u,v   
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Hydro_netcdfs_grid_options.f 
c  july 24,2002    Changes to old version to use the options and
c                  global variables.  Also the following comments.
c  Nov 11, 2002    Changes to bring into COARDS compliance
c                     time:units, sigma:positive=down, 
c                     global:conventions        :missing_value
c                     Allow Spaces in Strings
c
c  dec 27,2002  addition of nc_sync to save files even when model crashes
c  may 20,2003  little error in the wx,wy dimensions fixed by moving 
c               2D fields definations together
c  nov 26, 2003  addition of standard_name for COARDS
c
c feb 6, 2004  Rename dimensions of grid from lon,lat to nx,ny
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Real values in 3D 
      subroutine write_netcdf_Hydro(netcdf_file,ncid,imode,
     &  globalstr,m,n,l,
     & time,ibasedate,lon,lat,mask,sigma,depth,
     & zeta,u,v,w,temp,salt,wx,wy)

      implicit none

      include 'netcdf.inc'

c input variables
      real time,lon(m,n),lat(m,n),mask(m,n),sigma(l),depth(m,n),
     & zeta(m,n),u(m,n,l),v(m,n,l),w(m,n,l),
     & temp(m,n,l),salt(m,n,l),wx(m,n),wy(m,n)
      integer imode,ibasedate(4)
      character globalstr(9)*40
      character netcdf_file*80

C Netcdf internal variables
      integer m,n,l,num,itime
      
      integer iret,ncid,intval(4),CORNER(4),COUNT(4)
      integer nx_dim, ny_dim, sigma_dim, time_dim
      integer  time_id, lon_id, lat_id, mask_id, sigma_id
      integer  depth_id, zeta_id
      integer  u_id, v_id, w_id, salt_id, temp_id
      integer  wx_id, wy_id
      logical lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      
c date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
      
      save nx_dim, ny_dim, sigma_dim, time_dim
      save  time_id, lon_id, lat_id, mask_id, sigma_id
      save  depth_id, zeta_id
      save  u_id, v_id, w_id, salt_id, temp_id
      save  wx_id, wy_id
      save lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      
      
      if (imode.eq. 1) then
      write(*,*) 'write_netcdf_Hydro imode =',imode
C Set optional variable flags
      lz = .TRUE.
      lu = .TRUE.
      lv = .TRUE.
      lw = .TRUE.
      lsalt = .TRUE.
      ltemp = .TRUE.
      lwx = .TRUE.
      lwy = .TRUE.
      if (zeta(1,1).le.0.0) lz=.FALSE.
      if (u(1,1,1).le.0.0) lu=.FALSE.
      if (v(1,1,1).le.0.0) lv=.FALSE.
      if (w(1,1,1).le.0.0) lw=.FALSE.
      if (salt(1,1,1).le.0.0) lsalt=.FALSE.
      if (temp(1,1,1).le.0.0) ltemp=.FALSE.
      if (wx(1,1).le.0.0) lwx=.FALSE.
      if (wy(1,1).le.0.0) lwy=.FALSE.

      write(*,*) 'next =',imode

C Initialize

c         open file
      iret = nf_create(netcdf_file, NF_CLOBBER, ncid)
      call check_err(iret)
c define dimensions
      iret = nf_def_dim(ncid, 'nx', m, nx_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'ny', n, ny_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'sigma', l, sigma_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

c         define variables
C time
      iret = nf_def_var(ncid, 'time', NF_REAL, 1, time_dim, time_id)
      iret = nf_put_att_text(ncid, time_id, 'long_name', 4, 'Time')
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

C lon
      intval(2) = ny_dim
      intval(1) = nx_dim
      iret = nf_def_var(ncid, 'lon', NF_REAL, 2,intval, lon_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9, 'Longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lon_id, 'standard_name', 9, 
     & 'longitude')
C lat
      iret = nf_def_var(ncid, 'lat', NF_REAL, 2,intval, lat_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8, 'Latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lat_id, 'standard_name', 8, 
     & 'latitude')
C mask
      iret = nf_def_var(ncid, 'mask', NF_REAL, 2,intval, mask_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_id, 'long_name', 9,'Land Mask')
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_id,'units',14,'nondimensional')
      call check_err(iret)
C depth
      iret = nf_def_var(ncid, 'depth', NF_REAL, 2,intval, depth_id)
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
      call check_err(iret)
      iret = nf_put_att_text(ncid, sigma_id,'positive',4,'down')
      call check_err(iret)
       iret = nf_put_att_text(ncid, sigma_id, 'standard_name', 22, 
     & 'ocean_sigma_coordinate')
       iret = nf_put_att_text(ncid, sigma_id, 'formula_terms', 35, 
     & 'sigma: sigma eta: zeta depth: depth')
CC zeta
      intval(3) = time_dim
      intval(2) = ny_dim
      intval(1) = nx_dim
      if (lz) then
      iret = nf_def_var(ncid, 'zeta', NF_REAL, 3,intval, zeta_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'long_name', 23, 
     & 'Water Surface Elevation')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 6, 'meters')
      call check_err(iret)
      iret = nf_put_att_real(ncid, zeta_id, 
     & 'missing_value', NF_REAL, 1,-99999.0)
      iret = nf_put_att_text(ncid, zeta_id,'positive',2,'up')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 21, 
     & 'sea_surface_elevation')
      endif
C wind velocity East
      if (lwx) then
      iret = nf_def_var(ncid, 'air_u', NF_REAL, 3,intval, wx_id)
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
C wind velocity North
      if (lwy) then
      iret = nf_def_var(ncid, 'air_v', NF_REAL, 3,intval, wy_id)
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
      intval(4) = time_dim
      intval(3) = sigma_dim
      intval(2) = ny_dim
      intval(1) = nx_dim
      if (lu) then
      iret = nf_def_var(ncid, 'u', NF_REAL, 4,intval, u_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'long_name',
     & 23, 'Eastward Water Velocity')
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
      iret = nf_def_var(ncid, 'v', NF_REAL, 4,intval, v_id)
      iret = nf_put_att_text(ncid, v_id, 'long_name', 
     & 24, 'Northward Water Velocity')
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
      iret = nf_def_var(ncid, 'w', NF_REAL, 4,intval, w_id)
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
      iret = nf_def_var(ncid, 'temp', NF_REAL, 4,intval, temp_id)
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
      iret = nf_def_var(ncid, 'salt', NF_REAL, 4,intval, salt_id)
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
      
C Global Attributes
C Global Attributes
      iret = nf_put_att_text(ncid, NF_GLOBAL,'file_type'
     & ,9,'Full_Grid')
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

c         write lon,lat,mask,depth,sigma
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=m
      COUNT(2)=n
      iret=nf_put_vara_real(ncid,lon_id,CORNER,COUNT,lon)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,lat_id,CORNER,COUNT,lat)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,mask_id,CORNER,COUNT,mask)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,depth_id,CORNER,COUNT,depth)
      call check_err(iret)
      COUNT(1)=l
      iret=nf_put_vara_real(ncid,sigma_id,CORNER,COUNT,sigma)
      call check_err(iret)


      write(*,*) 'Initialize done',COUNT(1),COUNT(2)


C time index counter
      itime=0
C other indices
      
      elseif (imode.eq.2) then
c         write time,u,v,....
c Inquire of this file what the last itime written was
C  should be equal to the dimension of time
      iret= nf_inq_dimlen(ncid,time_dim,itime)      
      itime=itime+1
      write(*,*) 'write_netcdf_Hydro imode =',imode,
     & '  itime=',itime,'  ',netcdf_file

c scalars
      CORNER(1) = itime
      iret=nf_put_var1_real(ncid,time_id,CORNER,time)
      call check_err(iret)
c 2-D fields
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = itime
      COUNT(1)=m
      COUNT(2)=n
      COUNT(3)=1
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
c 3-D fields
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = 1
      CORNER(4) = itime
      COUNT(1)=m
      COUNT(2)=n
      COUNT(3)=l
      COUNT(4)=1
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

c  dec 27,2002  addition to keep intact files even when model crashes
       iret=nf_sync(ncid)
       call check_err(iret)
       


      else
c         close file
      iret = nf_close(ncid)
      call check_err(iret)
      write(*,*) 'write_netcdf_Hydro imode =',
     & imode,' close ',netcdf_file 
      endif
      
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3D files scaled to short integers 
      subroutine write_netcdf_Hydro_scale(netcdf_file,ncid,imode,
     & globalstr,m,n,l,
     & time,ibasedate,lon,lat,mask,sigma,depth,
     & zeta,u,v,w,temp,salt,wx,wy)

      implicit none

      include 'netcdf.inc'
c input variables
      integer imode,ibasedate(4)
      real time,lon(m,n),lat(m,n),mask(m,n),sigma(l),depth(m,n),
     & zeta(m,n),u(m,n,l),v(m,n,l),w(m,n,l),
     & temp(m,n,l),salt(m,n,l),wx(m,n),wy(m,n)
     
      character globalstr(9)*40
      character netcdf_file*80

C Netcdf internal variables
      integer m,n,l,num,itime
      integer iret,ncid,intval(4),CORNER(4),COUNT(4)
      integer nx_dim, ny_dim, sigma_dim, time_dim
      integer  time_id, lon_id, lat_id, mask_id, sigma_id
      integer  depth_id, zeta_id
      integer  u_id, v_id, w_id, salt_id, temp_id
      integer  wx_id, wy_id
      logical lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      real zeta_offset, zeta_scale, u_offset, u_scale
      real v_offset, v_scale, w_offset, w_scale, temp_offset
      real temp_scale, salt_offset, salt_scale, wx_offset
      real wx_scale, wy_offset, wy_scale
      
c date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)

      save nx_dim, ny_dim, sigma_dim, time_dim
      save  time_id, lon_id, lat_id, mask_id, sigma_id
      save  depth_id, zeta_id
      save  u_id, v_id, w_id, salt_id, temp_id
      save  wx_id, wy_id
      save lz,lu,lv,lw,lsalt,ltemp,lwx,lwy
      
c      write(*,*) 'scale options start imode =',imode
c  Total range of integers will be -scale/2:+scale/2 +offset  >> -32767:+32767 
c    integer = 32767* (real-offset)/(scale/2)
      zeta_offset = 0.0
      zeta_scale = 3.0
      u_offset = 0.0
      u_scale = 6.0
      v_offset = 0.0
      v_scale = 6.0
      w_offset = 0.0
      w_scale =  1.0
      temp_offset = 20.0
      temp_scale = 50.
      salt_offset = 20.
      salt_scale = 50.
      wx_offset = 0.0
      wx_scale = 60.0
      wy_offset = 0.0
      wy_scale = 60.0
                  
      
      if (imode.eq. 1) then
      write(*,*) 'write_netcdf_Hydro scale imode =',imode
C Set optional variable flags
      lz = .TRUE.
      lu = .TRUE.
      lv = .TRUE.
      lw = .TRUE.
      lsalt = .TRUE.
      ltemp = .TRUE.
      lwx = .TRUE.
      lwy = .TRUE.
      if (zeta(1,1).le.0.0) lz=.FALSE.
      if (u(1,1,1).le.0.0) lu=.FALSE.
      if (v(1,1,1).le.0.0) lv=.FALSE.
      if (w(1,1,1).le.0.0) lw=.FALSE.
      if (salt(1,1,1).le.0.0) lsalt=.FALSE.
      if (temp(1,1,1).le.0.0) ltemp=.FALSE.
      if (wx(1,1).le.0.0) lwx=.FALSE.
      if (wy(1,1).le.0.0) lwy=.FALSE.
      write(*,*) 'after Trues imode =',imode

C Initialize
c         open file
      iret = nf_create(netcdf_file, NF_CLOBBER, ncid)
      call check_err(iret)
c define dimensions
      iret = nf_def_dim(ncid, 'lon', m, nx_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'lat', n, ny_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'sigma', l, sigma_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

c         define variables
C time
      iret = nf_def_var(ncid, 'time', NF_REAL, 1, time_dim, time_id)
      iret = nf_put_att_text(ncid, time_id, 'long_name', 4, 'Time')
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

C lon
      intval(2) = ny_dim
      intval(1) = nx_dim
      iret = nf_def_var(ncid, 'lon', NF_REAL, 2,intval, lon_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9, 'Longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lon_id, 'standard_name', 9, 
     & 'longitude')
C lat
      iret = nf_def_var(ncid, 'lat', NF_REAL, 2,intval, lat_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8, 'Latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
       iret = nf_put_att_text(ncid, lat_id, 'standard_name', 8, 
     & 'latitude')
C mask
      iret = nf_def_var(ncid, 'mask', NF_INT, 2,intval, mask_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_id, 'long_name', 9,'Land Mask')
      call check_err(iret)
      iret = nf_put_att_text(ncid, mask_id,'units',14,'nondimensional')
      call check_err(iret)
C depth
      iret = nf_def_var(ncid, 'depth', NF_REAL, 2,intval, depth_id)
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
      call check_err(iret)
      iret = nf_put_att_text(ncid, sigma_id,'positive',4,'down')
      call check_err(iret)
       iret = nf_put_att_text(ncid, sigma_id, 'standard_name', 22, 
     & 'ocean_sigma_coordinate')
       iret = nf_put_att_text(ncid, sigma_id, 'formula_terms', 35, 
     & 'sigma: sigma eta: zeta depth: depth')
C zeta
      intval(3) = time_dim
      intval(2) = ny_dim
      intval(1) = nx_dim
      if (lz) then
      iret = nf_def_var(ncid, 'zeta', NF_SHORT, 3,intval, zeta_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'long_name', 23, 
     & 'Water Surface Elevation')
      call check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 6, 'meters')
      iret = nf_put_att_text(ncid, zeta_id,'positive',2,'up')
      call check_err(iret)
      call scale_attributes(ncid,zeta_id,zeta_scale,zeta_offset)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 21, 
     & 'sea_surface_elevation')
      endif
C wind velocity East
      if (lwx) then
      iret = nf_def_var(ncid, 'air_u', NF_SHORT, 3,intval, wx_id)
      iret = nf_put_att_text(ncid, wx_id, 'long_name', 21, 
     & 'Eastward Air Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, wx_id, 'units', 3, 'm/s')
      call check_err(iret)
      call scale_attributes(ncid,wx_id,wx_scale,wx_offset)
      iret = nf_put_att_text(ncid, wx_id, 'standard_name', 13, 
     & 'eastward_wind')
      endif
C wind velocity North
      if (lwy) then
      iret = nf_def_var(ncid, 'air_v', NF_SHORT, 3,intval, wy_id)
      iret = nf_put_att_text(ncid, wy_id, 'long_name', 22, 
     & 'Northward Air Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, wy_id, 'units', 3, 'm/s')
      call check_err(iret)
      call scale_attributes(ncid,wy_id,wy_scale,wy_offset)
      iret = nf_put_att_text(ncid, wy_id, 'standard_name', 14, 
     & 'northward_wind')
      endif
C u
      intval(4) = time_dim
      intval(3) = sigma_dim
      intval(2) = ny_dim
      intval(1) = nx_dim
      if (lu) then
      iret = nf_def_var(ncid, 'u', NF_SHORT, 4,intval, u_id)
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'long_name', 
     & 23, 'Eastward Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_id, 'units', 3, 'm/s')
      call check_err(iret)
c offset by 0.0 m/s and scale = 6.0m/s   = -3.0 : 3.0 m/s   
      call scale_attributes(ncid,u_id,u_scale,u_offset)
      iret = nf_put_att_text(ncid, u_id, 'standard_name', 27, 
     & 'eastward_sea_water_velocity')
      endif

C v
      if (lv) then
      iret = nf_def_var(ncid, 'v', NF_SHORT, 4,intval, v_id)
      iret = nf_put_att_text(ncid, v_id, 'long_name', 
     & 24, 'Northward Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'units', 3, 'm/s')
      call check_err(iret)
      call scale_attributes(ncid,v_id,v_scale,v_offset)
      iret = nf_put_att_text(ncid, v_id, 'standard_name', 28, 
     & 'northward_sea_water_velocity')
      endif
C w
      if (lw) then
      iret = nf_def_var(ncid, 'w', NF_SHORT, 4,intval, w_id)
      iret = nf_put_att_text(ncid, w_id, 'long_name', 23, 
     & 'Vertical Water Velocity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, w_id, 'units', 3, 'm/s')
      call check_err(iret)
      call scale_attributes(ncid,w_id,w_scale,w_offset)
      iret = nf_put_att_text(ncid, w_id, 'standard_name', 25, 
     & 'upward_sea_water_velocity')
      endif
C temp
      if (ltemp) then
      iret = nf_def_var(ncid, 'temp', NF_SHORT, 4,intval, temp_id)
      iret = nf_put_att_text(ncid, temp_id, 'long_name', 11, 
     & 'Temperature')
      call check_err(iret)
      iret = nf_put_att_text(ncid, temp_id, 'units', 7, 'Celsius')
      call check_err(iret)
      call scale_attributes(ncid,temp_id,temp_scale,temp_offset)
      iret = nf_put_att_text(ncid, temp_id, 'standard_name', 21, 
     & 'sea_water_temperature')
      endif
C salt
      if (lsalt) then
      iret = nf_def_var(ncid, 'salt', NF_SHORT, 4,intval, salt_id)
      iret = nf_put_att_text(ncid, salt_id, 'long_name', 8, 
     & 'Salinity')
      call check_err(iret)
      iret = nf_put_att_text(ncid, salt_id, 'units', 3, 'ppt')
      call check_err(iret)
      call scale_attributes(ncid,salt_id,salt_scale,salt_offset)
      iret = nf_put_att_text(ncid, salt_id, 'standard_name', 21, 
     & 'sea_water_salinity')
      endif

C Global Attributes
      iret = nf_put_att_text(ncid, NF_GLOBAL,'file_type'
     & ,9,'Full_Grid')
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

c         write lon,lat,mask,depth,sigma
      CORNER(1) = 1
      CORNER(2) = 1
      COUNT(1)=m
      COUNT(2)=n
      iret=nf_put_vara_real(ncid,lon_id,CORNER,COUNT,lon)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,lat_id,CORNER,COUNT,lat)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,mask_id,CORNER,COUNT,mask)
      call check_err(iret)
      iret=nf_put_vara_real(ncid,depth_id,CORNER,COUNT,depth)
      call check_err(iret)
      COUNT(1)=l
      iret=nf_put_vara_real(ncid,sigma_id,CORNER,COUNT,sigma)
      call check_err(iret)


      write(*,*) 'Initialize done',COUNT(1),COUNT(2)


C time index counter
      itime=0
C other indices
      
      elseif (imode.eq.2) then
c         write time,u,v,....
c Inquire of this file what the last itime written was
C  should be equal to the dimension of time
      iret= nf_inq_dimlen(ncid,time_dim,itime)      
      itime=itime+1
      write(*,*) 'write_netcdf_Hydro_scale imode =',imode,
     & '  itime=',itime,'  ',netcdf_file

c scalars
      CORNER(1) = itime
      iret=nf_put_var1_real(ncid,time_id,CORNER,time)
      call check_err(iret)
c 2-D fields
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = itime
      COUNT(1)=m
      COUNT(2)=n
      COUNT(3)=1
      if (lz) then
      call scale_writer_2d(ncid,zeta_id,CORNER,COUNT,
     &  zeta,m,n,zeta_scale,zeta_offset)
      endif
      if (lwx) then
      call scale_writer_2d(ncid,wx_id,CORNER,COUNT,
     &  wx,m,n,wx_scale,wx_offset)
      endif
      if (lwy) then
      call scale_writer_2d(ncid,wy_id,CORNER,COUNT,
     &  wy,m,n,wy_scale,wy_offset)
      endif

c 3-D fields
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = 1
      CORNER(4) = itime
      COUNT(1)=m
      COUNT(2)=n
      COUNT(3)=l
      COUNT(4)=1
c offset by 0.0 m/s and scale = 6.0m/s   = -3.0 : 3.0 m/s   
c      call scale_attributes(ncid,u_id,6.0,0.0)
c      iret=nf_put_vara_real(ncid,u_id,CORNER,COUNT,u)
c      call check_err(iret)
      if (lu) then
      call scale_writer_3d(ncid,u_id,CORNER,COUNT,
     &  u,m,n,l, u_scale,u_offset)
      endif
      if (lv) then
      call scale_writer_3d(ncid,v_id,CORNER,COUNT,
     &  v,m,n,l, v_scale,v_offset)
      endif
      if (lw) then
      call scale_writer_3d(ncid,w_id,CORNER,COUNT,
     &  w,m,n,l, w_scale,w_offset)
      endif
      if (ltemp) then
      call scale_writer_3d(ncid,temp_id,CORNER,COUNT,
     &  temp,m,n,l, temp_scale,temp_offset)
      endif
      if (lsalt) then
      call scale_writer_3d(ncid,salt_id,CORNER,COUNT,
     &  salt,m,n,l, salt_scale,salt_offset)
      endif


c  dec 27,2002  addition of nc_sync to save files even when model crashes
       iret=nf_sync(ncid)
       call check_err(iret)

      else
c         close file
      iret = nf_close(ncid)
      call check_err(iret)
      write(*,*) 'write_netcdf_Hydro_scale imode =',
     & imode,' close ',netcdf_file 
      endif
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c scaling attribute subroutines
c  add scaling attributes to variable attributes
      subroutine scale_attributes(ncid,var_id,scalefactor,offset)
      include 'netcdf.inc'
      integer var_id
      real offset,scalefactor
      range = 32767.-(-32767.)
      
      iret = nf_put_att_real(ncid, var_id, 'scale_factor', NF_REAL, 1, 
     & scalefactor/range)
      iret = nf_put_att_real(ncid, var_id, 'add_offset', NF_REAL, 1, 
     & offset)      

c.    missing value = minimum of range
c.    integer = -32767 converts to this value   
      iret = nf_put_att_real(ncid, var_id, 
     &  'missing_value', NF_REAL,1,(-scalefactor/2+offset))
      
      call check_err(iret)

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c write a 3D floating point as an integer using the scaling attributes
      subroutine scale_writer_3d(ncid,var_id,CORNER,COUNT,
     &  VAR,m,n,l,scalefactor,offset)
      real VAR(m,n,l)
      real offset,scalefactor,minvalue
      integer IVAR(m,n,l),CORNER(4),COUNT(4),var_id

      range = 32767.-(-32767.)
      minvalue = -scalefactor/2+offset      
      do i=1,m
        do j=1,n
          do k=1,l
           if (VAR(i,j,k).LE. minvalue) then
            IVAR(i,j,k) = -32767
           else
            IVAR(i,j,k) = NINT((VAR(i,j,k)-offset)/(scalefactor/range))
           endif
          enddo
        enddo
      enddo

      iret=nf_put_vara_int(ncid,var_id,CORNER,COUNT,IVAR)
      call check_err(iret)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c write a 2D floating point as an integer using the scaling attributes
      subroutine scale_writer_2d(ncid,var_id,CORNER,COUNT,
     &  VAR,m,n,scalefactor,offset)
      real VAR(m,n)
      real offset,scalefactor,minvalue
      integer IVAR(m,n),CORNER(4),COUNT(4),var_id

      range = 32767.-(-32767.)
      minvalue = -scalefactor/2+offset      
      
      do i=1,m
        do j=1,n
           if (VAR(i,j).LE.minvalue) then
            IVAR(i,j) = -32767
           else
            IVAR(i,j) = NINT((VAR(i,j)-offset)/(scalefactor/range))
           endif
        enddo
      enddo

      iret=nf_put_vara_int(ncid,var_id,CORNER,COUNT,IVAR)
      call check_err(iret)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     

