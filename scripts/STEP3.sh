#!/bin/sh
#Name:                 STEP3.sh
#purpose:              Make Tide Predictions using station info in 
#                     $STATIONDATA   [e.g. $HOME1/control_files/cbofsstationdata.ctl]
#
#Output:              data/prediction/$stationname".prd"
#                      $stationname is 4 letter name from the @STATIONDATA
#
# **********************************************************************************
#     STEP 3     Make tidal predictions
# **********************************************************************************
source ${HOME1}/scripts/STEPS_SETUP.sh
$SCRIPT_DIR/tide_prediction.sh

