#!/bin/sh
#Name:                 STEP9.sh
#purpose:             Run the Skill Assessment Computation
#                     Step10, Run the Harmonic Analysis Also
#
#                     Output file is ?
#
# **********************************************************************************
#     STEP 9    statistics computation and generate skill assessment score tables 
# **********************************************************************************
set -x
   source ${HOME1}/scripts/STEPS_SETUP.sh
   echo $BEGINDATE > skill.ctl
   echo $ENDDATE   >> skill.ctl
   echo $NTYPE $DELT $DELT_O $DELT_T $DELT_M $NCYCLE_F $CUTOFF $FACTOR  >> skill.ctl
   echo ${IS[1]} ${IS[2]} ${IS[3]} ${IS[4]} ${IS[5]} ${IS[6]} $ISG >> skill.ctl
   echo $IGAPFILL $CRITERIA1 $CRITERIA2 $METHOD  >> skill.ctl 
   echo $IPRT         >> skill.ctl 
   echo $DELHR $DELAMP $DELPCT $IOPTA   >> skill.ctl  
   echo $X1 $X2 $X11  >> skill.ctl 
   echo $KINDAT   >> skill.ctl
   echo $STATIONDATA >> skill.ctl
   echo $NFDURATION >> skill.ctl
   echo $FINCLUDE0 >> skill.ctl
   echo $DATA  >> skill.ctl
   $BIN/skills.x < skill.ctl 
   
