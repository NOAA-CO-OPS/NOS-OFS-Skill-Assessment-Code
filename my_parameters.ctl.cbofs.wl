################################################################################
#   BEGINDATE:   the begining date: "yyyy mm dd hh mn"
#   ENDDATE:     the end date:    "yyyy mm dd hh mn"
#   OS           =1, run in Unix, fortran codes are compiled using F90.
#                =0, run in Linux, Frotran codes are compiled using LF95. 
#   NTYPE        =0 for non-tidal region; =1 for tidal regions.
#   FACTOR       used for event selection of non-tidal region.
#                =<0.0 hupper=-factor,hlower=factor
#                >=0.0 hupper=mean +factor * SD ,hlower=mean -factor * SD.
#   DBASE        Four options: "NWLON", "PORTS", "GLAKES", and "USGS"
#  KINDAT      =1 for vector data (current speed and direction);
#              =2 for scalar data (water level)
#              =3 for temperature
#              =4 for salinity
# TZ_MODEL:    transfer time meridian of model outputs from local to GMT/UTC if necessary.
#              positive for west longitude, and negative for east longitude.
#              =0 if model output is in GMT/UTCl 
#              =75 if model output is in standard eastern time of US
#   NCYCLE_T:    number of tide simulation cycles per day, =0 read from a single file
#   NCYCLE_H:    number of hindcast cycles per day, =0 read from a single file
#   NCYCLE_N:    number of nowcast cycles per day, >=1
#   NCYCLE_F:    number of forecast cycles per day, >=1
#   NFDURATION   Forecast duration hours for skill assessment (in hours)
#   DELT:       desired time interval (in minutes) of observation, tide prediction and model outputs.
#   DELT_O:     actual time interval (in minutes) of observation, tide prediction
#   DELT_M:     actual time interval (in minutes) of model outputs.
#   CUTOFF:     CUTOFF period (in hours) for Fourier filtering, =0 no filtering
#   IGAPFILL:   control switch of gap filling with interpolation 
#               =0, filling with missing value -999.0;
#               =1, filling with interpolation value
#   METHOD:     index of interpolation method
#               0: cubic spline  1:Singular Value Decomposition(SVD); 
#   CRITERIA1:  (in hours)means using linear and cubic spline interpolation
#                when gap is less than criteria1    
#   CRITERIA2:  (in hours) means using cubic spline or SVD interpolation method
#                when criteria1 < gap < criteria2.
#               fill gaps using missing value -999.0 while gap > criteria2
#   IS:         control model run scenorios. =0,off; =1, on
#              IS(1): Tidal simulation only    
#              IS(2): model hindcast
#              IS(3): semi-operational nowcast
#              IS(4): semi-operational forecast
#              IS(5): forecast method comparison:persistence forecast
#              IS(6): forecast method comparison:tidal prediction
#  ISG:          Storm surge output switch. =0 no output; =1 output storm surg or non-tidal results
#  IPRT:       print switch. =0, no screen output; =1 screen output 
#  DELHR:      maximum allowed time difference between high and low (in hours),
#              if small than delhr, eliminate both high and low.
#  DELAMP:     maximum allowed amplitude difference between high and low (in meters),
#              if smaller, eliminate both
#  DELPCT:     maximum allowed fraction of amplitude difference between high and low
#  IOPTA:      option for selecting amplitude criterion
#              =1, delamp=delamp
#              =2, delamp=delpct*(hmax-hmin)
#              =3,delamp=delpct*(average hmax)
#  X1          accepted error criteria for water level (0.15 m),
#              current (0.26 m/s), temperature (7.5 c),salinity (3.5 psu)        
#              if X1 < 0, then X1 is defined by tide range in skill.f as 
#              X1 = ABS(TIDAL_RANGE*X1)/100.0
#              IF(X1 .LT. 0.15) X1 = 0.15, e.g. x1=-10, tide_range=3m, then x1=0.3m
#  X2          accepted error criteria for time (in hours)
#  X11         accepted error criteria for phase (in degrees)
# NCON        = number of constituents to be analyzed by H.A., maximum=37
# ofs         OFS name to be assessed (cbofs,dbofs,ngofs,etc.)
# variable    parameter to be assessed (e.g. wl,cu,temp, salt)
# FINCLUDED0     T: forecast station output include f000
#                F: forecast station output does not include f000 (POM_GLOFS or ADCIRC)
#****************************************************************
#           Specify project path names
#****************************************************************  
#   HOME1=/usr/global/people/zhang/SKILL08
#   HOME1=/home/net/azhang/SKILL_OFS
#****************************************************************
#          Specify required file names
#****************************************************************
ofs=cbofs
variable=wl

variable=`echo $variable | tr '[A-Z]' '[a-z]'  `
#echo $variable

#  station control file name

 #    STATIONDATA=$HOME1/control_files/${ofs}_temp_station.ctl
    STATIONDATA=$HOME1/control_files/${ofs}_${variable}_station.ctl

# Tide Simulations
   MODELTIDES=$HOME1/archive/ocean_sta_Nowd_ver18_tide_test_0.nc
#   MODELTIDES=$HOME1/archive/ocean_sta_ver18_tide_EC2001_adjusted_0.nc
# or NCYCLE_T >0
  ARCHIVE_DIR_T=../../TAMPA_BAY/SKILL_TB/archive
  ARCHIVE_DIR_T=$HOME1/archive
   NAME_TIDECAST="ocean_sta_constdensity"

# Hindcast Simulation
   HINDCAST=$HOME1/archive/ocean_sta_Nowd_ver18_hindcast_ADCIRC9_0.nc
# or NCYCLE_H >0
   ARCHIVE_DIR_H=../archive/TBOFS/ROMS/hindcast_EC2001_NoADD_FSOBC
#   ARCHIVE_DIR_H=$HOME1/archive
   NAME_HINDCAST="ocean_sta_Nowd_ver18_hindcast_EC2001_"

# Nowcast and Forecast Runs

ARCHIVE_DIR=/opt/ofs-bdp/${ofs}/netcdf
#ARCHIVE_DIR=/opt/archive/historical/${ofs}/netcdf
#ARCHIVE_DIR=/opt/archive/dev/${ofs}/netcdf
#ARCHIVE_DIR=/opt/prod/historical/glofs/${ofs}/netcdf
#ARCHIVE_DIR=/opt/ofsdev/users/azhang/SKILL_OFS/archive/cbofs/netcdf

NAME_NOWCAST="%Y%m/nos.${ofs}.stations.nowcast."%Y%m%d".t"%H"z.nc"
NAME_FORECAST="%Y%m/nos.${ofs}.stations.forecast."%Y%m%d".t"%H"z.nc"
#For POM-version GLOFS
#NAME_NOWCAST="%Y%m/glofs.${ofs}.stations.nowcast."%Y%m%d".t"%H"z.nc"
#NAME_FORECAST="%Y%m/glofs.${ofs}.stations.forecast."%Y%m%d".t"%H"z.nc"

if [ $variable = "CU" -o $variable  = "cu" ]
then
  KINDAT=1
  X1=0.26
  X2=0.50
elif [ $variable = "WL" -o $variable  = "wl" ]
then
  KINDAT=2
#  X1=-10 # 10% of tide range
  X1=0.15
  X2=0.50
elif [ $variable = "TEMP" -o $variable  = "temp" ]
then
  KINDAT=3
  X1=3.0
  X2=0.50
elif [ $variable = "SALT" -o $variable  = "salt" ]
then
  KINDAT=4
  X1=3.5
  X2=0.50
fi

#****************************************************************
#          Specify required parameters
#****************************************************************

   BEGINDATE="2021 12 01 00 00" 
     ENDDATE="2021 12 31 00 00"

   OS=0
   OCEAN_MODEL=ROMS
   NTYPE=1
   FACTOR=2.0
   DBASE=PORTS
   TZ_MODEL=0
   NCYCLE_T=1
   NCYCLE_H=1
   NCYCLE_N=4
   NCYCLE_F=4
   NFDURATION=48
   DELT=6
   DELT_O=6
   DELT_M=6
   DELT_T=6
   CUTOFF=0.0
   IGAPFILL=1
   METHOD=1
   CRITERIA1=2
   CRITERIA2=5
   IS[1]=0
   IS[2]=0
   IS[3]=1
   IS[4]=1
   IS[5]=0
   IS[6]=0
   ISG=0
   IPRT=0 
   DELHR=2.0
   DELAMP=0.2
   DELPCT=0.03
   IOPTA=3
   X11=22.5 
   NCON=23
  FINCLUDE0=F   
#****************************************************************





