
# This script is used to compile all of the fortran and C programs
# While OS=0 for Linux, and OS=1 for Unix
#  Run as: COMPILE.sh 0
#  Aijun Zhang, CSDL
#  Jan. 15, 2005
BIN=../bin
if test ! -r $BIN
then
   mkdir -p $BIN
fi   
FC=ifort 
CC=gcc
NETCDF_ROOT=/usr/local

#NETCDFLIB="-I${NETCDF_ROOT}/include -L${NETCDF_ROOT}/lib -lnetcdff -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1 -lgcc"
NETCDFLIB="-I${NETCDF_ROOT}/include -L${NETCDF_ROOT}/lib -lnetcdff"
$FC -extend-source READUSGS.f -o $BIN/READUSGS.x
$FC -extend-source READCDEC.f -o $BIN/READCDEC.x
$FC -extend-source -O gaps0.f utility.f -o $BIN/gaps0.x
$FC -extend-source -O fillnan.f utility.f -o $BIN/fillnan.x
$FC -extend-source -O harm15d.f svd.f equal_interval.f utility.f -o $BIN/harm15d.x
$FC -extend-source -O harm29d.f svd.f equal_interval.f utility.f -o $BIN/harm29d.x
$FC -extend-source -O lsqha.f svd.f equal_interval.f utility.f -o $BIN/lsqha.x

$FC -extend-source -O harm15d.notimelimits.f svd.f equal_interval.f utility.f -o $BIN/harm15d.notimelimits.x
$FC -extend-source -O harm29d.notimelimits.f svd.f equal_interval.f utility.f -o $BIN/harm29d.notimelimits.x
$FC -extend-source -O lsqha.notimelimits.f svd.f equal_interval.f utility.f -o $BIN/lsqha.notimelimits.x

$FC -extend-source -O persistence.f utility.f  -o $BIN/persistence.x 
$FC -extend-source -O persistence1.f utility.f  -o $BIN/persistence1.x 
$FC -extend-source -O pred.f utility.f -o $BIN/pred.x 
$FC -extend-source -O read_netcdf_fcst1.f read_netcdf_subs.f utility.f  -o $BIN/read_netcdf_fcst1.x  $NETCDFLIB 
$FC -extend-source -O read_netcdf_fcst.f  read_netcdf_subs.f utility.f -o $BIN/read_netcdf_fcst.x $NETCDFLIB
$FC -extend-source -O read_netcdf_modeltides.f read_netcdf_subs.f utility.f   -o $BIN/read_netcdf_modeltides.x  $NETCDFLIB 
$FC -extend-source -O read_netcdf_nowcast.f read_netcdf_subs.f utility.f  -o $BIN/read_netcdf_nowcast.x $NETCDFLIB
$FC -extend-source -O reformat_ha.f  -o $BIN/reformat_ha.x 
$FC -extend-source -O refwl.f  utility.f -o $BIN/refwl.x  
$FC -extend-source -O reformat_PORTS.f utility.f   -o $BIN/reformat_PORTS.x  
$FC -extend-source -O reformat_USGS.f  utility.f  -o $BIN/reformat_USGS.x  
$FC -extend-source -O skills.f  foufil.f extremes.f svd.f equal_interval1.f slack.f table_tide.f table_nontide.f utility.f -o $BIN/skills.x  
$FC -extend-source -O table_Harmonic_C.f  -o $BIN/table_Harmonic_C.x  
$FC -extend-source -O write_obs_netcdf.f  Hydro_netcdfs_station.f  Hydro_netcdfs_grid.f -o $BIN/write_obs_netcdf.x  $NETCDFLIB  
$FC -extend-source -O filter.f foufil.f utility.f  -o $BIN/filter.x
$CC datemath.c -o $BIN/datemath -lm
$CC dateformat.c -o $BIN/dateformat
chmod 755 $BIN/*
rm -f *.o
