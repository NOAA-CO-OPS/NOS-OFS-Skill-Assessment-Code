#!/bin/sh
#Name:                harmonic_analysis_obs.sh 
#purpose:             perform harmonic analysis for the obseerved tidal time series.
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
source STEPS_SETUP.sh
cd $WRK_DIR
rm -f cons.out
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel 
do
    read LAT LONGITUDE XMAJOR  DEPTH1 DEPTH2
    FILEIN1=../data/obs/$stationname'.obs'
    FILEIN2=../data/obs/$stationname'B.obs'
    FILEOUT=../data/obs/$stationname'VG.obs'
    if [ -s $FILEIN1 -a -s $FILEIN2 ]
    then
        $BIN/gradient.x $FILEIN1 $FILEIN2 $FILEOUT $DEPTH1 $DEPTH2
    fi
    FILEIN1=../model_FVCOM/$stationname'_hindcast.dat'
    FILEIN2=../model_FVCOM/$stationname'B_hindcast.dat'
    FILEOUT=../model_FVCOM/$stationname'VG_hindcast.dat'
    if [ -s $FILEIN1 -a -s $FILEIN2 ]
    then
        $BIN/gradient.x $FILEIN1 $FILEIN2 $FILEOUT $DEPTH1 $DEPTH2 
    fi
    FILEIN1=../model_ROMS/$stationname'_hindcast.dat'
    FILEIN2=../model_ROMS/$stationname'B_hindcast.dat'
    FILEOUT=../model_ROMS/$stationname'VG_hindcast.dat'
    if [ -s $FILEIN1 -a -s $FILEIN2 ]
    then
        $BIN/gradient.x $FILEIN1 $FILEIN2 $FILEOUT $DEPTH1 $DEPTH2 
    fi
    FILEIN1=../model_POM/$stationname'_hindcast.dat'
    FILEIN2=../model_POM/$stationname'B_hindcast.dat'
    FILEOUT=../model_POM/$stationname'VG_hindcast.dat'
    if [ -s $FILEIN1 -a -s $FILEIN2 ]
    then
        $BIN/gradient.x $FILEIN1 $FILEIN2 $FILEOUT  $DEPTH1 $DEPTH2
    fi
    FILEIN1=../model_SELFE/$stationname'_hindcast.dat'
    FILEIN2=../model_SELFE/$stationname'B_hindcast.dat'
    FILEOUT=../model_SELFE/$stationname'VG_hindcast.dat'
    if [ -s $FILEIN1 -a -s $FILEIN2 ]
    then
        $BIN/gradient.x $FILEIN1 $FILEIN2 $FILEOUT $DEPTH1 $DEPTH2 
    fi

done 3<&-

exit
