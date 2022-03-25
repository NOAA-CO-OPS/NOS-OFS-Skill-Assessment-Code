#!/bin/sh
#Name:                 STEP7.sh
#purpose:              Concatenate ForeCast simulation runs 
#Output:	       work/$stationname".dat
#
# **********************************************************************************
#     STEP 7    concatenate model forecast simulation 
# **********************************************************************************

  source ${HOME1}/scripts/STEPS_SETUP.sh
  if [ ${IS[4]} -eq 1 ]; then 
     echo ' processing model forecast simulations'
     $SCRIPT_DIR/concatenate_forecast.sh

     exec 5<&0 <$STATIONDATA
     while read stnid stationname longlabel
     do
        read Latlon
        if [ -s $stationname ]; then
          mv $stationname $stationname'_forecast.dat'
        else
          echo "The station" $stationname "is not output in the model forecast."
        fi
     done 3<&-
     echo ' forecast data were concatenated'
  else
    echo " Warning: The IS(4) is not equal to 1, this process is terminated!"
  fi

