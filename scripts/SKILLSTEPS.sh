#!/bin/sh
#Name:               SKILLSTEPS.sh
#purpose:            Runs the STEP*.sh scripts in sequence 
#                    
#Output:	
#
#Author:             Aijun Zhang
export HOME1=/opt/ofsdev/users/azhang/SKILL_OFS
# **********************************************************************************
#    Run all the Steps
# **********************************************************************************
# **********************************************************************************
#     STEP 2     observation CO-OPS verified 6-minutes water level
# **********************************************************************************
    $HOME1/scripts/STEP2.sh

# **********************************************************************************
#     STEP 3     Make tidal predictions
# **********************************************************************************
#     $HOME1/scripts/STEP3.sh

# **********************************************************************************
#     STEP 4    read model tidal simulation 
# **********************************************************************************
#    $HOME1/scripts/STEP4.sh
 
# **********************************************************************************
#     STEP 5    read model hindcast simulation 
# **********************************************************************************
#    $HOME1/scripts/STEP5.sh

# **********************************************************************************
#     STEP 6    concatenate model nowcast simulation 
# **********************************************************************************
    $HOME1/scripts/STEP6.sh

# **********************************************************************************
#     STEP 7    concatenate model forecast simulation 
# **********************************************************************************
    $HOME1/scripts/STEP7.sh

# **********************************************************************************
#     STEP 8    make persistence forecasts 
# **********************************************************************************
 #   $HOME1/scripts/STEP8.sh

# **********************************************************************************
#     STEP 9    statistics computation and generate skill assessment score tables 
# **********************************************************************************
    $HOME1/scripts/STEP9.sh

# *******************************************************************************************
#     STEP 10:  conduct harmonical constants comparison between the observed and the modeled
# *******************************************************************************************
#   $HOME1/scripts/STEP10.sh
