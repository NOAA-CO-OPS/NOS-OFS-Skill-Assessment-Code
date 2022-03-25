#!/bin/sh
#---------------------------------------------------------------------------------------
#
# Script Name: get_TS_USGS.sh
#
# Abstract:
#           Gets data from the USGS Real-time temperature and salinity data from web page:
#  
#		 https://waterdata.usgs.gov/md/nwis/uv?01578310 
# 
#	    Request tab separated data and you will see the source file.
# 
#	    There is no choice about times on this web page, so this only
#	    gives you the last SEVEN days of data. 
#
#	    The script decodes these files to grab the different data
#	    types which might be available.  Not all stations have the
#	    same data (or in the same order.)  Possible choices are:
#
#  	    TEMP        TEMPERATURE, WATER (DEG. C)
# 	    COND        SPECIFIC CONDUCTANCE (MICROSIEMENS/CM AT 25 DEG. C)
#  	    DISCHARGE   DISCHARGE, CUBIC FEET PER SECOND
# 	    GAGE        GAGE HEIGHT, FEET
#
#	    The requested page is sent to  READUSGS.pl to parse out the data
#	    type requested.              
#		produces a ascii file
# 		Capable of returning any data variables from any river station
#	    Returns ascii like:
#	    2003 04 17 00 00  0  6.390000 63700.000000
#	    2003 04 17 00 30  0  6.380000 63600.000000
#
# Usage:  Interactively:    river_read_usgs.sh stationid listvar startdate enddate outputfilename
#			    river_read_usgs.sh  01570500 "GAGE DISCHARGE" "2003 04 01 00 00" \
#				"2003 04 26 12 00" riv.dat
#              Via cron:    Called by RIVERQCF.sh, SALTQCF.sh, TEMPQCF.sh, WLQCF.sh 
#
# Input Parameters:     stationid number  Ex.01570500
#                       list of variables  Ex.  "GAGE DISCHARGE"
#                       starting time  Ex  "2003 04 01 00 00"
#                       ending time  Ex  "2003 04 26 12 00"
#			output file name  Ex  riv.dat
#
# Language:  Bourne Shell Script
#
# Target Computer:    COMF computer, such as dsofs1.nos-tcn.noaa.gov
#
# Estimated Execution Time: 
#
# Scripts/Programs Called:
#    Name               Directory Location             Description
#   READUSGS.pl		/COMF/oqcs/scripts	Parses out the data from USGS river web page.
#   mktemp.c		/COMF/oqcs/sorc		Makes a temporary unique filename.
#
# Input Files:
#    Name               Directory  Location            Description
#
# Output Files:
#    Name               Directory Location             Description
#   riv.dat		user defined		2003 04 17 00 00  0  6.390 63700.00
#
# Error Conditions:
#
# Orginal coded by AJ
#
# Modified by Lianyuan Zheng on 03/01/2017

cd $WRK_DIR
echo get_TS_USGS.sh $BEGINDATE " to " $ENDDATE
BEGINDATE0=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
ENDDATE0=`$BIN/datemath $ENDDATE + 0 0 2 0 0`
WGETOUT=`mktemp -q wgetout.XXXXXX`
WGETLOG=`mktemp -q wgetlog.XXXXXX`
RIVSCRATCH=`mktemp -q river1.dat.XXXXXX`

BEGINDATE1=${BEGINDATE0%%??????}
ENDDATE1=${ENDDATE0%%??????}

BEGINDATE2=${BEGINDATE1// /-}  # replace all blank space with - 
ENDDATE2=${ENDDATE1// /-}

#echo $BEGINDATE2 $ENDDATE2
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel
do
  echo " " StationNames $stnid ":" $stationname ":" $longlabel
  read Latlon
# Per Susan Trapanese, USGS 7-21-2004, use the nwis. ..... url if the waterdata.usgs..... is down
#REQUESTGET="http://waterdata.usgs.gov/nwis/uv?format=rdb&period=31&site_no=$stnid"
#REQUESTGET="http://nwis.waterdata.usgs.gov/usa/nwis/uv/?&format=rdb&period=400&begin_date=2014-04-16&end_date=2014-04-23&site_no=08017118
#REQUESTGET="http://nwis.waterdata.usgs.gov/usa/nwis/uv/?&format=rdb&period=450&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
#REQUESTGET="http://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00095=on&cb_00010=on&cb_90860=on&format=rdb&period=450&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
#REQUESTGET="http://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00010=on&cb_00095=on&format=rdb&period=400&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
#REQUESTGET="http://nwis.waterdata.usgs.gov/usa/nwis/uv/?&format=rdb&period=450&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"

#  Zheng modified on 05/03/2018
#  REQUESTGET="https://nwis.waterdata.usgs.gov/usa/nwis/uv/?&cb_00010=on&cb_00095=on&format=rdb&period=450&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
  REQUESTGET="https://nwis.waterdata.usgs.gov/usa/nwis/uv/?&cb_00010=on&cb_00095=on&cb_00060=on&cb_00065=on&cb_90860=on&cb_00300=on&cb_63680=on&cb_99409=on&format=rdb&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
#  echo $REQUESTGET
  wget --no-check-certificate -o $WGETLOG -O $WGETOUT $REQUESTGET

# AJ Zhang, 11/20/2009, use READUSGS.f to replace READUSGS.pl
  $BIN/READUSGS.x $WGETOUT "COND TEMP" $RIVSCRATCH
#$HOME1/scripts/READUSGS.pl $WGETOUT "COND TEMP" $RIVSCRATCH

  rm -f *tmp*
  cat $RIVSCRATCH |   \
  awk '{ print $1" "$2" "$3" "$4" "$5"  " $6 "  "$7*(1+0.02*($8-25)) "  " $8  }' | sort -u > tmp1
  awk '{ print $1" "$2" "$3" "$4" "$5 "  " $8  }' tmp1 | sort -u > temperature.tmp

# Call the PERL script to convert conductivity and temperature to Salinity
  $HOME1/scripts/SALINITY.pl tmp1 tmp2
  awk '{printf("%4d %02d %02d %02d %02d   %10.5f \n",$1,$2,$3,$4,$5,$7)}' tmp2 > junk0

  if [ $KINDAT -eq 3 ]; then
    $BIN/reformat_USGS.x "$BEGINDATE0" "$ENDDATE0" $KINDAT temperature.tmp tmp.out
  elif [ $KINDAT -eq 4 ]; then
    rm -f tmp.salt salinity.tmp tmp.out
    $BIN/READUSGS.x $WGETOUT "SALT" tmp.salt
    if [ -s tmp.salt ]; then
      echo "Data from Salinity"
      awk '{ print $1" "$2" "$3" "$4" "$5 "  " $7 }'  tmp.salt | sort -u > salinity.tmp
      $BIN/reformat_USGS.x "$BEGINDATE0" "$ENDDATE0" $KINDAT  salinity.tmp tmp.out
    else 
      echo "Data from Conductivity"
      $BIN/reformat_USGS.x "$BEGINDATE0" "$ENDDATE0" $KINDAT junk0  tmp.out
    fi
  fi

  if [ -s tmp.out ]; then
    mv tmp.out $OBS/$stationname".obs"
  fi    
done
echo "done of get_TS_USGS.sh!!!"

rm junk* tmp*
rm $RIVSCRATCH &> /dev/null
rm $WGETOUT &> /dev/null
rm $WGETLOG &> /dev/null
exit

