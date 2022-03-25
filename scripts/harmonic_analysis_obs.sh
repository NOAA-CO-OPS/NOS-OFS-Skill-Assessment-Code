#!/bin/sh
#Name:                harmonic_analysis_obs.sh 
#purpose:             perform harmonic analysis for the observed tidal time series.
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
source ${HOME1}/scripts/STEPS_SETUP.sh
cd $WRK_DIR
rm -f cons.out
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel 
do
    read LAT LONGITUDE XMAJOR  LAYER
    FILEIN=${DATA}/obs/$stationname'.obs'
    echo $FILEIN
    if [ -s $FILEIN ]
    then
      $BIN/gaps0.x $FILEIN  $KINDAT $DELT_O
     read max_time< tmp.dat
     echo  max_time=$max_time
     read TBEGIN < tmp1.dat
     echo 'BEGIN TIME=' $TBEGIN
     read TEND< tmp2.dat
     echo 'END TIME=' $TEND
      if [ $max_time -ge 180 ]
      then
        NCON=29
      elif [ $max_time -ge 30 -a $max_time -lt  180 ]
      then
         NCON=23
      else
         NCON=16
      fi	 
      	
      if [ $max_time -ge 40 ]
      then
#        $BIN/lsqha.notimelimits.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN 
        $BIN/lsqha.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN "$TBEGIN" "$TEND"
      elif [ $max_time -lt 40 -a $max_time -ge 29 ]
      then
#        $BIN/harm29d.notimelimits.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN   
        $BIN/harm29d.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN "$TBEGIN" "$TEND"  
      elif [ $max_time -lt 29 -a $max_time -ge 15 ]
      then
#        $BIN/harm15d.notimelimits.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN   
        $BIN/harm15d.x $KINDAT  $NCON $DELT_O $LONGITUDE $FILEIN "$TBEGIN" "$TEND"  
      elif [ $max_time -lt 15 -a $max_time -ge 5 ]
      then
#        $BIN/lsqha.notimelimits.x $KINDAT  8 $DELT_O $LONGITUDE $FILEIN 
        $BIN/lsqha.x $KINDAT  16 $DELT_O $LONGITUDE $FILEIN "$TBEGIN" "$TEND"
       echo 'time series is too short,  perform HA with less constituents lsqha.f'
      fi   
      if [ -s cons.out ]
      then
        mv cons.out ../data/harmonic_con/$stationname'.std'
      fi
      
    fi
done 3<&-

exit
