#!/bin/sh
#
# Documentation for Scripts get_data_ndbc_currents.sh 
#
#---------------------------------------------------------------------------------------
#
# Script Name: get_data_ndbc_currents.sh 
#
# Directory Location:   /COMF/oqcs/scripts  
# 
# Technical Contact(s): Name:  Aijun Zhang            Org:    NOS/CO-OPS
#                       Phone: 301-713-2890x127     E-Mail: aijun.zhang@noaa.gov
#
# Abstract:
#             
#		get_data_ndbc_currents.sh  gets data from NDBC and produces a ascii file
# 		Although this does sort for the date range, it does not go
# 		back in time more than one month.  ie it reads the realtime web site.

#		Output file has date, forecasthour, temp
       #2011 11 14 01 00 00    70.00    12.00
       #2011 11 14 01 30 00    70.00    11.00
       #2011 11 14 23 00 00    10.00    16.00

# 		uses dateformat, datemath and wget 
# 		Put them into the path with exports in maincalling scripts
# 		PATH=$PATH:.:/COMF/oqcs/bin:/COMF/oqcs/scripts
#
# Usage:   Interactively:  get_data_ndbc_currents.sh stationid startdate enddate outputfilename
#                          get_data_ndbc_currents.sh TPLM2 "2003 03 09 00 00" "2003 03 12 12 00" TPLM2.txt
#               Via cron:  Called by TEMPQCF.sh
#
# Input Parameters:     Station name  Ex. TPLM2
#                       starting time  Ex. "2003 03 09 00 00"
#                       ending time  Ex.  "2003 03 12 12 00"
#                       output file name  Ex.  TPLM2.txt
#
# Language:  Bourne Shell Script
#
# Target Computer:   Runs on COMF computers, such as dsofs1.nos-tcn.noaa.gov.  
#
# Estimated Execution Time: 
#
# Scripts/Programs Called:
#    Name                       Directory Location              Description
#   dateformat			/COMF/oqcs/sorc		        Flexible String builder using dates.
#   mktemp.c			/COMF/oqcs/sorc  		Makes a temporary unique filename.
#   wget			/COMF/oqcs/binlinux		Request data from WWW web.
#
# Input Files:
#    Name                       Directory  Location            Description
#    none
#
# Output Files:
#    Name                       Directory Location             Description
#    outputfilename		user defined		       year  m  d  h  m  spd  dir
#
# Error Conditions:
#
# Author Name:  Tom Gross       Creation Date:  
#
# Revisions:
#         Date          Author         Description
#	08-19-2004	H Lin 	     Standardized the documentation
#       11-19-2004      H Lin	     Add called scripts.
#	03-28-2005      H Lin        Update document.
#	11-03-2016      L Zheng      Modify script to download historical and real time data.
#
# Remarks:    
#	Gets data from the NDBC web site.  Really only tested on the
#  	CMAN stations, but ought to work for the buoys also.
#  
#	Their web page has fixed file names for the data files
#	so this works nicely:
#	wget http://www.ndbc.noaa.gov/data/realtime2/$stnid.txt -O $WGETOUT
#
#	But a limited amount of data is available.  
#	Only the most recent  45 days.  This has been extended to download the historical data modified by L. Zheng
#	
#	As a screen scraper this is susceptible to changes.  The format of
#	the downloaded files has changed in the past.  
#	
#	These files have 
#YY   MM DD hh mm  DEP01 DIR01 SPD01 DEP02 DIR02 SPD02 DEP03 DIR03 SPD03 DEP04 DIR04 SPD04 DEP05 DIR05 SPD05 DEP06 DIR06 SPD06 DEP07 DIR07 SPD07
#                    m   degT  cm/s     m  degT  cm/s     m  degT  cm/s     m  degT  cm/s     m  degT  cm/s
#2011 11 16 18 30    2    60    12
#	and are listed in upside down order.  Whatever.
#	They are converted to 
#	y m d h m fh temp
# -----------------------------------------------------------------------------------------------------------------
# Modified by Lianyuan Zheng on 03/01/2017

stnid=$1
begindate="$2"
enddate="$3"
fileout=$4
echo " " get_data_ndbc_currents.sh $begindate " to " $enddate

# temporary unique file names
WGETOUT=`mktemp -q wgetout.XXXXXX`
AWKED1=`mktemp -q awked1.XXXXXX`
AWKED=`mktemp -q awked.XXXXXX`
SCRATCH=`mktemp -q CURRSCRATCH.XXXXXX`

####  Begin of adding by Zheng
##  Get the years and months of begindate, enddate, and present date
Y_B=$(echo $begindate | cut -c1-4)
M_B=$(echo $begindate | cut -c6-7)
Y_E=$(echo $enddate | cut -c1-4)
M_E=$(echo $enddate | cut -c6-7)
export Y_P=$(date +%Y)
export M_P=$(date +%m)

#  convert stnid from uppercase to lowercase
stnidlow=$(echo $stnid | tr '[:upper:]' '[:lower:]')

#get the historical data
export WEBcur=${WEBcur:-https://www.ndbc.noaa.gov/data/historical/adcp}

#  Download the previous years data
for (( int_year=Y_B; int_year < Y_P; int_year=int_year+1)); do
  var1=`echo $int_year |  awk '{printf("%04i",$1)}'`

# download the adcp historical data  
  web_cur=${WEBcur}/${stnidlow}a${var1}.txt.gz
  wget -q ${web_cur} -O temp.gz
  actualsize=$(wc -c <"temp.gz")
  if [ $actualsize -gt 0 ]; then
    gunzip temp.gz
    cat temp >> $WGETOUT
    if [ -s temp ]; then
      rm -f temp
    fi
  else
    rm -f temp.gz
  fi
done

#  Download the present years data
#  Note: There are two adcp folders on ndbc website
#  If the requested end date is in the present year, do the following.  Otherwise skip it!

if [ $Y_E -eq $Y_P ]; then
  export WEBcur1=${WEBcur1:-https://www.ndbc.noaa.gov/data/adcp}
  export WEBcur2=${WEBcur2:-https://www.ndbc.noaa.gov/data/adcp2}
  for (( int_mon=1; int_mon <= 12; int_mon=int_mon+1)); do
    if [ $int_mon -eq 1 ]; then
      monname='Jan'
    elif [ $int_mon -eq 2 ]; then
      monname='Feb'
    elif [ $int_mon -eq 3 ]; then
      monname='Mar'
    elif [ $int_mon -eq 4 ]; then
      monname='Apr'
    elif [ $int_mon -eq 5 ]; then
      monname='May'
    elif [ $int_mon -eq 6 ]; then
      monname='Jun'
    elif [ $int_mon -eq 7 ]; then
      monname='Jul'
    elif [ $int_mon -eq 8 ]; then
      monname='Aug'
    elif [ $int_mon -eq 9 ]; then
      monname='Sep'
    elif [ $int_mon -eq 10 ]; then
      monname='Oct'
    elif [ $int_mon -eq 11 ]; then
      monname='Nov'
    else
      monname='Dec'
    fi

    web_cur1=${WEBcur1}/${monname}/${stnidlow}${int_mon}${Y_P}.txt.gz
    web_cur2=${WEBcur2}/${monname}/${stnidlow}${int_mon}${Y_P}.txt.gz

    wget -q ${web_cur1} -O temp.gz
    actualsize=$(wc -c <"temp.gz")
    if [ $actualsize -gt 0 ]; then
      gunzip temp.gz
      cat temp >> $WGETOUT
      if [ -s temp ]; then
        rm -f temp
      fi
    else
      rm -f temp.gz
    fi
    unset web_cur1

    wget -q ${web_cur2} -O temp.gz
    actualsize=$(wc -c <"temp.gz")
    if [ $actualsize -gt 0 ]; then
      gunzip temp.gz
      cat temp >> $WGETOUT
      if [ -s temp ]; then
        rm -f temp
      fi
    else
      rm -f temp.gz
    fi
    unset web_cur2

    web_cur1=${WEBcur1}/${monname}/${stnidlow}${int_mon}${Y_P}.txt
    web_cur2=${WEBcur2}/${monname}/${stnidlow}${int_mon}${Y_P}.txt

    wget -q ${web_cur1} -O temp
    actualsize=$(wc -c <"temp")
    if [ $actualsize -gt 0 ]; then
      cat temp >> $WGETOUT
      if [ -s temp ]; then
        rm -f temp
      fi
    else
      rm -f temp
    fi
    unset web_cur1

    wget -q ${web_cur2} -O temp
    actualsize=$(wc -c <"temp")
    if [ $actualsize -gt 0 ]; then
      cat temp >> $WGETOUT
      if [ -s temp ]; then
        rm -f temp
      fi
    else
      rm -f temp
    fi
    unset web_cur2
  done
####  End of adding by Zheng

#  THE WWW REQUEST  
  wget -q https://www.ndbc.noaa.gov/data/realtime2/$stnid.adcp -O temp

####  Begin of adding by Zheng
  actualsize=$(wc -c <"temp")
  if [ $actualsize -gt 0 ]; then
    tac temp > temp1
    cat temp1 >> $WGETOUT
    if [ -s temp ]; then
      rm -f temp temp1
    fi
  else
    rm -f temp
  fi

#else
#  echo It is not present year!
fi
####  End of adding by Zheng

DATEF=`$BIN/dateformat $begindate "%Y%m%d%H"`
DATEL=`$BIN/dateformat $enddate "%Y%m%d%H"`
awk -v datef=$DATEF '$1$2$3$4 >= datef {print }'  $WGETOUT > $AWKED1
awk -v datel=$DATEL '$1$2$3$4 <= datel {print }'  $AWKED1 > $AWKED

#  modified by Lianyuan Zheng  -- convert speed from cm/s to knot and unique increase sort
awk '$7 != "MM" { printf("%04i %02i %02i %02i %02i %8.2f %8.2f \n",$1,$2,$3,$4,$5,$8/51.444,$7)  }' \
    $AWKED | sort -u > $SCRATCH
#  end of modified by Lianyuan Zheng

if [ -s $SCRATCH ]; then
  mv $SCRATCH $fileout
fi  

rm -f $WGETOUT
rm -f $AWKED1
rm -f $AWKED
rm -f $SCRATCH
echo "  done of get_data_ndbc_current.sh!!!"
exit

