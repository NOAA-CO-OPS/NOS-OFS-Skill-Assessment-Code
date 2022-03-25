#!/bin/sh
#Name:                 STEPS_SETUP.sh
#purpose:              READ my_parameters.ctl and 
#                      set environment variables for STEP*.sh scripts
#
#*********************************************************************************************** 
#   STEP 1  modify two control files: my_parameters.ctl and station information control file
# **********************************************************************************
##. /disks/NASUSER/azhang/SKILL_OFS/control_files/my_parameters.ctl
. $HOME1/my_parameters.ctl
WRK_DIR=$HOME1/work
#WRK_DIR=` pwd `/work
SCRIPT_DIR=$HOME1/scripts
#DATA=` pwd `/data
DATA=$HOME1/data
CTL=$HOME1/control_files
OBS=$DATA/obs
PRD=$DATA/prediction
CONSTANTS_DIR=$DATA/harmonic_con
LOG=$HOME1/log
BIN=$HOME1/bin
if test ! -r $WRK_DIR 
then
   mkdir -p $WRK_DIR
fi
if test ! -r $LOG 
then
   mkdir -p $LOG
fi
if test ! -r $OBS
then
   mkdir -p $OBS
fi
if test ! -r $PRD
then
   mkdir -p $PRD
fi
if test ! -r $CONSTANTS_DIR
then
   mkdir -p $CONSTANTS_DIR
fi
export HOME1 WRK_DIR LOG SCRIPT_DIR BIN SORC_DIR CONSTANTS_DIR OBS PRD CTL DATA
export BEGINDATE ENDDATE NCYCLE_T NCYCLE_H NCYCLE_N NCYCLE_F DELT DELT_O DELT_M DELT_T
export CUTOFF IGAPFILL METHOD CRITERIA1 CRITERIA2 NFDURATION
export DBASE KINDAT NCON NTYPE FACTOR DATUM
export ARCHIVE_DIR_T ARCHIVE_DIR_H ARCHIVE_DIR
export STATIONDATA MODELTIDES HINDCAST NAME_TIDECAST NAME_HINDCAST NAME_NOWCAST NAME_FORECAST
export TZ_MODEL OCEAN_MODEL
export ofs variable
#BEGINDATE=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
#ENDDATE=`$BIN/dateformat $ENDDATE "%Y %m %d 00 00"`
cd $WRK_DIR
