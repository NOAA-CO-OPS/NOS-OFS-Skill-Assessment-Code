#!/bin/sh
#Name:                harmonic_analysis.sh 
#purpose:             perform harmonic analysis for the model simulated tidal time series.
# Author:             Aijun Zhang
# Date:               11/20/2004  
#Language:            Korn Shell Script
#input parameters:    BEGINDATE,ENDDATE,ARCHIVE_DIR,SUBNAME_NOWCAST,NCYCLE_N,DELT,KINDAT,STATIONDATA
# Programs Called:
#          Name               Location                            Description
#         lsqha.x               $BIN                FORTRAN program for least squares harmonic analysis
#         harm29d.x             $BIN               FORTRAN program for Fourier harmonic analysis for 29 days time series
#     table_Harmonic_C.x        $BIN               FORTRAN program to create constituents comparison tables 
#                                                  between the observed and modeled values  
#set -x
cd $WRK_DIR
rm -f cons.out
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel 
do
    read LAT LONGITUDE XMAJOR SDEPTH
    FILEIN=$stationname'_modeltides.dat'
    echo $FILEIN
    if [ -s $FILEIN ]
    then
      $BIN/gaps0.x $FILEIN  $KINDAT $DELT_M
      read max_time < tmp.dat
      echo $max_time
      if [ $max_time -ge 40 ]
      then
          $BIN/lsqha.x $KINDAT  $NCON $DELT_M $LONGITUDE $FILEIN "$BEGINDATE" "$ENDDATE" 
      elif [ $max_time -lt 40 -a $max_time -ge 29 ]
      then
        $BIN/harm29d.x $KINDAT  $NCON $DELT_M $LONGITUDE $FILEIN "$BEGINDATE" "$ENDDATE" 
      elif [ $max_time -lt 29 -a $max_time -ge 15 ]
      then
        $BIN/harm15d.x $KINDAT  $NCON $DELT_M $LONGITUDE $FILEIN "$BEGINDATE" "$ENDDATE" 
      elif [ $max_time -lt 15 ]
      then
       echo 'time series is too short, do not perform HA'
      fi
      if [ -s cons.out ]
      then
         mv cons.out $stationname'_modeltides.std'
      fi
    fi   
    file_obs='../data/harmonic_con/'$stationname'.std'
    file_model=$stationname'_modeltides.std'
    if [ -s $file_obs -a -s $file_model ]
    then
      $BIN/table_Harmonic_C.x $KINDAT  $stationname "$longlabel" $DATA
    fi
done 3<&-

exit
