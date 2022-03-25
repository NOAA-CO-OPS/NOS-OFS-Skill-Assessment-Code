#!/bin/sh
# Name:                tide_prediction.sh 
# purpose:             make tidal predictions between BEGINDATE and ENDDATE
# Author:              Aijun Zhang
#   Coast Survey Development LaboratorySDL, NOS of NOAA
# Date:                11/20/2004  
# Language:            Korn Shell Script
# input parameters:    BEGINDATE,ENDDATE, DELT,KINDAT,STATIONDATA
# Programs Called:
#          Name               Location                            Description
#    reformat_ha.x            $BIN                 Fortran program to reformat harmonic constants 
#                                                  to a standard format
#      pred.x                 $BIN                 Fortran Program to make tidal predictions
# Modified by Lianyuan Zheng on 03/01/2017

 if [ $KINDAT -gt 2 ]; then
   echo KINDAT must be smaller than 3
   exit
 fi

 echo  " run pred.x from $BEGINDATE  to  $ENDDATE"
 BEGINDATE0=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
 ENDDATE0=`$BIN/datemath $ENDDATE + 0 0 2 0 0`
 bdate=`$BIN/dateformat $BEGINDATE0 "%Y%m%d"`
 edate=`$BIN/dateformat $ENDDATE0 "%Y%m%d"`

 exec 5<&0 <$STATIONDATA
 while read stnid stationname longlabel 
 do
   read lat longitude XMAJOR sdepth
   echo Station: $stationname
   if [ ! -s $CONSTANTS_DIR/$stationname'.std' ]; then
     echo $CONSTANTS_DIR/$stationname'.std' does not exist !!!
     if [ $KINDAT -eq 2 ]; then

##  Getting datum 
       sed -e s/VSTNIDV/$stnid/ \
          $CTL/request.template_datum_verified > request.GET_datum
       wget --no-check-certificate -o junk -O $stationname'.datum' -i request.GET_datum

##  Get Harmonic Constants from CO-OPS
       sed -e s/VSTNIDV/$stnid/ \
          $CTL/request.template_HC_verified > request.GET_ha
       wget --no-check-certificate -o junk -O $stationname'.ha' -i request.GET_ha
       $BIN/reformat_ha.x $stationname

       if [ -s $stationname'.std' ]; then
         mv $stationname'.std' $CONSTANTS_DIR
       fi
     fi
   fi

##  Run pred.f
  if [ -s $CONSTANTS_DIR/$stationname'.std' ]; then
    FILEIN=$CONSTANTS_DIR/$stationname'.std'
    FILEOUT=$stationname.prd
#    echo $FILEIN
    $BIN/pred.x "$BEGINDATE0" "$ENDDATE0" $KINDAT $DELT $FILEIN $FILEOUT
    mv $stationname.prd $PRD
  fi
 done 3<&-
 echo "Done of tide_prediction.sh!!!"
exit

