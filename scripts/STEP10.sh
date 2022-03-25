#!/bin/sh
#Name:                 STEP9.sh
#purpose:             Run the Skill Assessment Computation
#                     Step10, Run the Harmonic Analysis Also
#
#                     Output file is ?
#

#source /disks/NASUSER/azhang/SKILL_OFS/scripts/STEPS_SETUP.sh
source ${HOME1}/scripts/STEPS_SETUP.sh
# *******************************************************************************************
#     STEP 10:  conduct harmonical constants comparison between the observed and the modeled
# *******************************************************************************************
if [ ${IS[1]} -eq 1 ]
then 
    $SCRIPT_DIR/harmonic_analysis.sh
fi
