#!/bin/sh
#Name:                 STEP8.sh
#purpose:              Make Persistence Forecast Cases 
#                     Output file is ?
#
# **********************************************************************************
#     STEP 8    make persistence forecasts 
# **********************************************************************************

  source ${HOME1}/scripts/STEPS_SETUP.sh
  if [ ${IS[5]} -eq 1 -a ${IS[4]} -eq 1 ]; then 
    echo makeing persistent forecasts 
    exec 5<&0 <$STATIONDATA
    while read stnid stationname longlabel 
    do
      read lat longitude XMAJOR
      echo $OBS/$stationname'.obs'> persistence.ctl
      echo $PRD/$stationname'.prd'>> persistence.ctl
      echo $stationname'_forecast.dat'>> persistence.ctl
      echo $stationname'_persistence.dat'>> persistence.ctl
      echo $KINDAT >> persistence.ctl
      echo $BEGINDATE >> persistence.ctl
      echo $ENDDATE >> persistence.ctl
      echo $DELT >> persistence.ctl
      echo $NCYCLE_F >> persistence.ctl
      echo $NFDURATION >> persistence.ctl
      $BIN/persistence.x < persistence.ctl > $LOG/skills.log
    done 3<&-
  fi

  if [ ${IS[5]} -eq 1 -a ${IS[4]} -ne 1 ]; then 
    echo makeing persistent forecasts 
    exec 5<&0 <$STATIONDATA
    while read stnid stationname longlabel 
    do
      read lat longitude XMAJOR
      echo $OBS/$stationname'.obs'> persistence.ctl
      echo $PRD/$stationname'.prd'>> persistence.ctl
      echo $stationname'_persistence.dat'>> persistence.ctl
      echo $KINDAT >> persistence.ctl
      echo $BEGINDATE >> persistence.ctl
      echo $ENDDATE >> persistence.ctl
      echo $DELT >> persistence.ctl
      echo $NCYCLE_F >> persistence.ctl
      echo $NFDURATION >> persistence.ctl
      $BIN/persistence1.x < persistence.ctl > $LOG/skills.log
    done 3<&-
  fi

