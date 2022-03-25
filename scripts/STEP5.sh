#!/bin/sh
#Name:                 STEP5.sh
#purpose:              Read or Concatenate HindCast simulation runs 
#                      Output file is ?
#
# **********************************************************************************
#     STEP 5    read model hindcast simulation 
# **********************************************************************************
  source ${HOME1}/scripts/STEPS_SETUP.sh

  if [ ${IS[2]} -eq 1 ]; then 
   echo processing model hindcast simulations
     if [ $NCYCLE_H -lt 1 ]; then 
       $BIN/read_netcdf_modeltides.x "$BEGINDATE" $HINDCAST $STATIONDATA $KINDAT $OCEAN_MODEL $DELT_M
     else
       $SCRIPT_DIR/concatenate_hindcast.sh $ARCHIVE_DIR_H $NAME_HINDCAST $NCYCLE_H 
     fi 

     exec 5<&0 <$STATIONDATA
     while read stnid stationname longlabel
     do
        read Latlon
        mv $stationname $stationname'_hindcast.dat'
     done 3<&-
     echo ' Hindcast simulation data were read'
  else
    echo "Warning: The IS(2) is not equal to 1, this process is terminated!"
  fi

