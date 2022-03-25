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
    read LAT LONGITUDE XMAJOR  LAYER
    FILEIN=../data/obs/$stationname'.obs'
    FILEOUT=../data/obs/$stationname'_LP.obs'
    FILEIN=../work_SELFE/$stationname'_hindcast.dat'
    FILEOUT=../work_SELFE/$stationname'_LP_hindcast.dat'
    if [ -s $FILEIN ]
    then
        $BIN/filter.x $KINDAT  $DELT $CUTOFF $FILEIN $FILEOUT
    fi
done 3<&-

exit
