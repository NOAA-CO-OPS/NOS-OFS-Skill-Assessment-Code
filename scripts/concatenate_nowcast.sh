#!/bin/sh
# Name:                concatenate_nowcast.sh
# Purpose:             concatenate all netCDF files of the nowcast cycles 
#                        between BEGINDATE and ENDDATE
# Author:              Aijun Zhang
#                      Coast Survey Development LaboratorySDL, NOS of NOAA
# Date:                11/20/2004  
# Language:            Korn Shell Script
# Input parameters:    BEGINDATE,ENDDATE,ARCHIVE_DIR,NAME_NOWCAST,NCYCLE_N,
#                      DELT,KINDAT,STATIONDATA
# Programs Called:
#       Name             Location                      Description
# read_netcdf_now.x       $BIN         FORTRAN program to read nowcasts from a netCDF file
#********************************************************************************************

  DIR=$1
  NAME=$2
  CYCLES=$3
  BEGINDATE1=$BEGINDATE
  ENDDATE1=`$BIN/datemath $ENDDATE + 0 0 0 1 0`
  filename0='blank'

  rm -f filename.ctl
  while [ `$BIN/dateformat $BEGINDATE1 "%Y%m%d%H"` -le `$BIN/dateformat $ENDDATE1 "%Y%m%d%H"` ]
  do
  #   total length of filename cannot exceed 200 characters, otherwise, there will be an error
      filename=$DIR/`$BIN/dateformat $BEGINDATE1 $NAME`
    if [ -s $filename  -a  $filename != $filename0  ]; then
       echo $filename >> filename.ctl
       filename0=$filename
    fi
      BEGINDATE1=`$BIN/datemath $BEGINDATE1 + 0 0 0 0 10`
  done

  wc -l filename.ctl > junk
  read N nn < junk
  echo $BEGINDATE > file.ctl
  echo $DELT_M $CYCLES $KINDAT $N >> file.ctl
  echo $STATIONDATA >> file.ctl
  echo $TZ_MODEL >> file.ctl
  echo $OCEAN_MODEL >> file.ctl
  cat file.ctl filename.ctl > tmp1
  cp tmp1 now_filename.ctl
  rm -f tmp1
  $BIN/read_netcdf_nowcast.x < now_filename.ctl
  exit

