#!/bin/sh
# Name:                concatenate_forecast.sh
# Purpose:             concatenate all netCDF files of the forecast cycles 
#                        between BEGINDATE and ENDDATE
# Author:              Aijun Zhang
# Date:                11/20/2004  
# Language:            Korn Shell Script
# Input parameters:    BEGINDATE,ENDDATE,ARCHIVE_DIR,SUBNAME_FORECAST,NCYCLE_F,
#                      DELT,KINDAT,STATIONDATA
# Programs Called:
#        Name            Location                     Description
# read_netcdf_fcst1.x     $BIN     FORTRAN program to read netCDF file and check whether  
#                                    it contains 24-h forecasts.
# read_netcdf_fcst.x      $BIN     FORTRAN program to read 24-h forecasts from a netCDF file
#********************************************************************************************* 

  ENDDATE1=`$BIN/datemath $ENDDATE + 0 0 0 1 0`
  BEGINDATE1=`$BIN/dateformat $BEGINDATE "%Y %m %d %H %M"`
  filename0='blank' 

  if [ -s filename.ctl ]; then
    rm -f filename.ctl
  fi

  while [ `$BIN/dateformat $BEGINDATE1 "%Y%m%d%H"` -le `$BIN/dateformat $ENDDATE1 "%Y%m%d%H"` ]
  do 
     BEGINDATE2=$BEGINDATE1
     cycle=0
     count=0

     while (( count < 144))    # incremental 10 minutes each loop
     do
       filename1=$ARCHIVE_DIR/`$BIN/dateformat $BEGINDATE2 $NAME_FORECAST`
       if [ -s $filename1  -a  $filename1 != $filename0 ]; then
         filename[cycle]=$filename1
         time0[cycle]=$BEGINDATE2
         echo $DELT_M $NCYCLE_F $KINDAT 1 $NFDURATION> filetmp.ctl
         echo $STATIONDATA >> filetmp.ctl
         echo $OCEAN_MODEL >> filetmp.ctl
         echo ${filename[cycle]} >> filetmp.ctl
         echo `$BIN/dateformat $BEGINDATE2 "%Y %m %d %H %M" ` >> filetmp.ctl

         $BIN/read_netcdf_fcst1.x < filetmp.ctl
         read dummy < fort.86
         if [ $dummy = 'F' ]; then 
           echo $filename1
           echo ' The file does not contain correct data.'
           break
         fi
         (( cycle = cycle + 1 ))
         filename0=$filename1
       fi
       (( count = count + 1 ))
       BEGINDATE2=`$BIN/datemath $BEGINDATE2 + 0 0 0 0 10`
     done

     if [ $cycle -eq $NCYCLE_F ]; then
        cycle=0
        while (( cycle < $NCYCLE_F))
        do
           echo ${filename[cycle]} >> filename.ctl
           echo `$BIN/dateformat ${time0[cycle]} "%Y %m %d %H %M" ` >> filename.ctl
          (( cycle = cycle + 1 ))
        done
     fi 
     BEGINDATE1=`$BIN/datemath $BEGINDATE1 + 0 0 1 0 0`
  done

  wc -l filename.ctl > junk
  read N nn < junk
  (( N = N / 2 ))
  echo $BEGINDATE > filetmp.ctl
  echo $DELT_M $NCYCLE_F $KINDAT $N $NFDURATION >> filetmp.ctl
  echo $STATIONDATA >> filetmp.ctl
  echo $TZ_MODEL >> filetmp.ctl
  echo $OCEAN_MODEL >> filetmp.ctl
  cat filetmp.ctl filename.ctl > tmp1
  cp tmp1 fore_filename.ctl
  rm -f tmp1
  $BIN/read_netcdf_fcst.x < fore_filename.ctl
  exit

