#!/bin/sh
#Name:                 STEP6.sh
#purpose:              Concatenate NowCast simulation runs 
#Output:	       work/$stationname".dat
#
# **********************************************************************************
#     STEP 6    concatenate model nowcast simulation 
# **********************************************************************************

  source ${HOME1}/scripts/STEPS_SETUP.sh
  if [ ${IS[3]} -eq 1 ]; then 
    echo processing model nowcast simulations
    $SCRIPT_DIR/concatenate_nowcast.sh $ARCHIVE_DIR $NAME_NOWCAST $NCYCLE_N 

    exec 5<&0 <$STATIONDATA
    while read stnid stationname longlabel
    do
        read Latlon
        if [ -s $stationname ]; then
          mv $stationname $stationname'_nowcast.dat'
        else
          echo "The station" $stationname "is not output in the model nowcast."
        fi
    done 3<&-
    echo ' nowcast data were concatenated'
  else
    echo "Warning: The IS(3) is not equal to 1, this process is terminated!"
  fi 

