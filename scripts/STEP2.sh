#!/bin/sh
#Name:                STEP2.sh
#purpose:             Get Observations using station info in 
#                     $STATIONDATA   [e.g. $HOME1/control_files/cbofsstationdata.ctl]
#
#                     Output is ../data/obs/$stationname".obs" 
#                     $stationname is string read from $STATIONDATA
#
# Author:             Aijun Zhang
# Date:               11/20/2004  
#Language:            Korn Shell Script
# Modified:           Lianyuan Zheng on 03/01/2017

# **********************************************************************************
#     STEP 2     observation CO-OPS verified 6-minutes water level
# **********************************************************************************
  source ${HOME1}/scripts/STEPS_SETUP.sh

  if [ $KINDAT = 1 ]; then
    if [ $DBASE = "PORTS" -o  $DBASE = "NWLON" -o  $DBASE = "GLAKES" ]; then
      $SCRIPT_DIR/get_currents_COOPS.sh
    else 
      echo 'There is no data reader for' $DBASE 'database!'
      exit 
    fi

  elif [ $KINDAT = 2 ]; then
    if [ $DBASE = "PORTS" -o  $DBASE = "NWLON" -o  $DBASE = "GLAKES" ]; then
      $SCRIPT_DIR/get_WL_historic_COOPS.sh
    else
      echo 'Water level station is not in CO-OPS database!'
      exit 
    fi

  elif [ $KINDAT = 3 -o $KINDAT = 4 ]; then
    $SCRIPT_DIR/get_TS_USGS.sh
    $SCRIPT_DIR/get_TS_COOPS.sh

  else
    echo 'The value of parameter KINDA is not correct!'
    exit 
  fi

