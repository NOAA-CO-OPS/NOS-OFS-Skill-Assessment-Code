#!/bin/sh
#
# Script Name: get_TS_CDEC.sh
#
# Abstract:
#           Gets data from the USGS Real-time temperature and salinity data from web page:
#  
#		 http://waterdata.usgs.gov/md/nwis/uv?01578310 
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
# Modified:           Lianyuan Zheng on 03/01/2017 

echo get_TS_CDEC.sh $BEGINDATE " to " $ENDDATE
stnid=$1 
tstart=$2
tend=$3

WGETOUT1=`mktemp -q wgetout1.XXXXXX`
WGETLOG1=`mktemp -q wgetlog1.XXXXXX`
RIVSCRATCH1=`mktemp -q river1.dat.XXXXXX`

WGETOUT2=`mktemp -q wgetout2.XXXXXX`
WGETLOG2=`mktemp -q wgetlog2.XXXXXX`
RIVSCRATCH2=`mktemp -q river2.dat.XXXXXX`

BEGINDATE1=${tstart%%??????}
ENDDATE1=${tend%%??????}

BEGINDATE2=${BEGINDATE1// /-}  # replace all blank space with - 
ENDDATE2=${ENDDATE1// /-}

# Per Susan Trapanese, USGS 7-21-2004, use the nwis. ..... url if the waterdata.usgs..... is down
# REQUESTGET="http://waterdata.usgs.gov/nwis/uv?cb_00095=on&cb_00010=on&format=rdb&period=&begin_date=$BEGINDATE2&end_date=$ENDDATE2&site_no=$stnid"
#  For Conductivity, sensor_num=100; for Temperature, sensor_num=25
REQUESTGET="http://cdec.water.ca.gov/cgi-progs/selectQuery?station_id=$stnid&sensor_num=100&dur_code=H&start_date=$BEGINDATE2&end_date=$ENDDATE2"
wget --no-check-certificate -o $WGETLOG1 -O $WGETOUT1 $REQUESTGET

# AJ Zhang, 11/20/2009, use READUSGS.f to replace READUSGS.pl
# Machuan cut TEMP $BIN/READUSGS.x $WGETOUT "COND TEMP" $RIVSCRATCH
$BIN/READCDEC.x $WGETOUT1 "COND" $RIVSCRATCH1
cond_file=$RIVSCRATCH1

REQUESTGET="http://cdec.water.ca.gov/cgi-progs/selectQuery?station_id=$stnid&sensor_num=25&dur_code=H&start_date=$BEGINDATE2&end_date=$ENDDATE2"
wget --no-check-certificate -o $WGETLOG2 -O $WGETOUT2 $REQUESTGET
$BIN/READCDEC.x $WGETOUT2 "TEMP" $RIVSCRATCH2
temp_file=$RIVSCRATCH2

cat $RIVSCRATCH2 | awk '{ print ($7-32)*5.0/9.0}' > tmp11
paste  $RIVSCRATCH1  tmp11 > tmp12 

# https://www.solinst.com/products/dataloggers-and-telemetry/3001-levelogger-series/operating-instructions/user-guide/1-introduction/1-2-4-conductivity.php 
# to see the detail for conversion of electronic conductivity to specific conductivity
# Specific Conductance = Conductivity / (1 + 0.02 * (temp(C) - 25))
cat  tmp12 | awk '{ print $1" "$2" "$3" "$4" "$5"  " $6 "  "$7*(1+0.02*($8-25)) "  " $8  }'   > tmp1  
awk '{printf("%4d %02d %02d %02d %02d  %10.5f \n",$1,$2,$3,$4,$5,$8)}' tmp1  > temperature.tmp

$HOME1/scripts/SALINITY.pl tmp1 tmp2
awk '{printf("%4d %02d %02d %02d %02d  %10.5f \n",$1,$2,$3,$4,$5,$7)}' tmp2 > salinity.tmp

rm -f junk* tmp*
rm -f $RIVSCRATCH1 &> /dev/null
rm -f $RIVSCRATCH2 &> /dev/null
rm -f $WGETOUT1 &> /dev/null
rm -f $WGETOUT2 &> /dev/null
rm -f $WGETLOG1 &> /dev/null
rm -f $WGETLOG2 &> /dev/null
exit

