#!/bin/sh
#  Name:                 STEP4.sh
#  purpose:              Read tide simulation run 
#  Input:   MODELTIDES=
#  /comf/development/COMFmod_dev/ohms/CBOFS/worktide/meccastation.nc
#  A single file for a long tidal simulation
# 
#  Output:               work/$stationname_modeltides.dat		
#  Modified by Lianyuan Zheng on 03/01/2017

# ******************************************************************************
#     STEP 4    read model tidal simulation 
# ******************************************************************************

source ${HOME1}/scripts/STEPS_SETUP.sh
if [ ${IS[1]} -eq 1 ]; then 
  echo processing tidal model simulations
  if [ $NCYCLE_T -lt 1 ]; then 
    echo $STATIONDATA
    $BIN/read_netcdf_modeltides.x "$BEGINDATE" $MODELTIDES $STATIONDATA $KINDAT $OCEAN_MODEL $DELT_M
  else
    $SCRIPT_DIR/concatenate_hindcast.sh $ARCHIVE_DIR_T $NAME_TIDECAST $NCYCLE_T 
  fi 

  exec 5<&0 <$STATIONDATA
  while read stnid stationname longlabel
   do
     read Latlon
     echo StationNames $stnid ":" $stationname ":" $longlabel
     mv $stationname $stationname'_modeltides.dat'
   done 3<&-
   echo ' Tidal simulation data were read'
else
  echo "Warning: The IS(1) is not equal to 1, this process is terminated!"
fi

