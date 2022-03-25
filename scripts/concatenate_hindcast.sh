#!/bin/sh
#Name:                concatenate_nowcast.sh
#purpose:             concatenate all netCDF files of the nowcast cycles between BEGINDATE and ENDDATE
# Author:             Aijun Zhang
#   Coast Survey Development LaboratorySDL, NOS of NOAA
# Date:               11/20/2004  
#Language:            Korn Shell Script
#input parameters:    BEGINDATE,ENDDATE,ARCHIVE_DIR,NAME_NOWCAST,NCYCLE_N,DELT,KINDAT,STATIONDATA
# Programs Called:
#          Name               Location                            Description
#  read_netcdf_now.x          $BIN               FORTRAN program to read nowcasts form a netCDF file
#*********************************************************************************************** 

  DIR=$1
  NAME=$2
  CYCLES=$3
  BEGINDATE1=$BEGINDATE
  ENDDATE1=`$BIN/datemath $ENDDATE + 0 0 0 1 0`
  rm -f filename.ctl

  FILE=$DIR'/*'$NAME'*'
  NFILES=`ls -al $FILE | wc -l`
  echo $BEGINDATE > filename.ctl
  echo $DELT_M $CYCLES $KINDAT $NFILES >> filename.ctl
  echo $STATIONDATA >> filename.ctl
  echo $TZ_MODEL >> filename.ctl
  echo $OCEAN_MODEL >> filename.ctl
  ls -al $FILE | awk '{ print $9 }' >> filename.ctl
  $BIN/read_netcdf_nowcast.x < filename.ctl
  echo ' hindcast data were concatenated'
  exit

