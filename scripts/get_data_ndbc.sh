#!/bin/sh
#
# Documentation for Scripts get_data_ndbc.sh 
#
#---------------------------------------------------------------------------------------
#
# Script Name: get_data_ndbc.sh 
#
# Directory Location:   /COMF/oqcs/scripts  
# 
# Technical Contact(s): Name:  Aijun zhang            Org:    NOS/CO-OPS
#                       Phone: 301-713-2890x127     E-Mail: aijun.zhang@noaa.gov
#
# Abstract:
#             
#		get_data_ndbc.sh  gets data from NDBC and produces a ascii file
# 		Although this does sort for the date range, it does not go
# 		back in time more than one month.  ie it reads the realtime web site.
#
#		Gets data from the NDBC web site.  
#		This actually gets any of the variables:
# 		AT   Air Temperature     Centigrade
# 		WT   Water Temperature   Centegrade
# 		AP   Air Pressure        mbars   1022.6                 
#
#		Output file has date, forecasthour, temp
#  		 (ndbc, nwlon, tide forecasthour ==0 )
#  		y    m  d  h  m   fh   temp
# 		2002 12 30 12 30   0   25.5
#
# 		uses dateformat, datemath and wget 
# 		Put them into the path with exports in maincalling scripts
# 		PATH=$PATH:.:/COMF/oqcs/bin:/COMF/oqcs/scripts
#
# Usage:  Interactively:   get_data_ndbc.sh stationid startdate enddate sensor outputfilename
#                          get_data_ndbc.sh  TPLM2 "2003 03 09 00 00" "2003 03 12 12 00" AT tempTPLM2.txt
#              Via cron:   Called by TEMPQCF.sh
#
# Input Parameters:     Station name  Ex. TPLM2
#                       starting time  Ex. "2003 03 09 00 00"
#                       ending time  Ex.  "2003 03 12 12 00"
#                       output file name  Ex.  tempTPLM2.txt
#			sensor  Ex.   AT, air temperature
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
#    outputfilename		user defined		       year  m  d  h  m   fh   temp
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
#
# Remarks:    
#	Gets data from the NDBC web site.  Really only tested on the
#  	CMAN stations, but ought to work for the buoys also.
#  
#	Their web page has fixed file names for the data files
#	so this works nicely:
#	wget https://www.ndbc.noaa.gov/data/realtime/$stnid.txt -O $WGETOUT -o $WGETLOG
#
#	But a limited amount of data is available.  
#	Only the most recent  45 days 
#	
#	As a screen scraper this is susceptible to changes.  The format of
#	the downloaded files has changed in the past.  
#	
#	These files have 
#	YYYY MM DD hh WD  WSPD  GST  WVHT   DPD   APD MWD  BARO   ATMP  WTMP  DEWP  VIS PTDY  TIDE
#	2003 03 27 18 190  2.1  2.6    MM    MM    MM  MM 1021.9  10.7   8.6   4.8   MM -0.5    MM
#	and are listed in upside down order.  Whatever.
#	They are converted to 
#	y m d h m fh temp
#	
# -----------------------------------------------------------------------------------------------------------------
# Modified by Lianyuan Zheng on 03/01/2017

stnid=$1
begindate="$2"
enddate="$3"
sensor=$4
fileout=$5
rm -f $fileout
echo get_data_ndbc.sh $sensor $begindate to $enddate

####  Temporary unique file names
WGETOUT=`mktemp -q wgetout.XXXXXX`
WGETOUT1=`mktemp -q wgetout1.XXXXXX`
WGETLOG=`mktemp -q wgetlog.XXXXXX`
AWKED1=`mktemp -q awked1.XXXXXX`
AWKED=`mktemp -q awked.XXXXXX`
TMPF=`mktemp -q tmpf.XXXXXX`

####  Begin of adding by Zheng
####  Get the years and months and months of begindate, enddate, and present date
Y_B=$(echo $begindate | cut -c1-4)
M_B=$(echo $begindate | cut -c6-7)
Y_E=$(echo $enddate | cut -c1-4)
M_E=$(echo $enddate | cut -c6-7)
export Y_P=$(date +%Y)
export M_P=$(date +%m)

####  Convert stnid from uppercase to lowercase
stnidlow=$(echo $stnid | tr '[:upper:]' '[:lower:]')

####  Get the historical data
export WEBmet=${WEBmet:-https://www.ndbc.noaa.gov/data/historical/stdmet}
export WEBoce=${WEBoce:-https://www.ndbc.noaa.gov/data/historical/ocean}

####  Download the previous years data
for (( int_year=Y_B; int_year < Y_P; int_year=int_year+1)); do
  var1=`echo $int_year |  awk '{printf("%04i",$1)}'`

####  Download the standard meteorological data  
  web_met=${WEBmet}/${stnidlow}h${var1}.txt.gz
  wget -q ${web_met} -O temp.gz
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

####  Download the oceanographical data  
  web_oce=${WEBoce}/${stnidlow}o${var1}.txt.gz
  wget -q ${web_oce} -O temp.gz
  actualsize=$(wc -c <"temp.gz")
  if [ $actualsize -gt 0 ]; then
    gunzip temp.gz
    cat temp >> $WGETOUT1
    if [ -s temp ]; then
      rm -f temp
    fi
  else
    rm -f temp.gz
  fi
done

####  Download the present years data
export WEBmet1=${WEBmet1:-https://www.ndbc.noaa.gov/data/stdmet}
export WEBoce1=${WEBoce1:-https://www.ndbc.noaa.gov/data/ocean}
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
  
  web_met1=${WEBmet1}/${monname}/${stnidlow}${int_mon}${Y_P}.txt.gz
  wget -q ${web_met1} -O temp.gz
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
  unset web_met1
  
  web_met1=${WEBmet1}/${monname}/${stnidlow}.txt
  wget -q ${web_met1} -O temp
  actualsize=$(wc -c <"temp")
  if [ $actualsize -gt 0 ]; then
    cat temp >> $WGETOUT
    if [ -s temp ]; then
      rm -f temp
    fi
  else
    rm -f temp
  fi
  unset web_met1
  
  web_oce1=${WEBoce1}/${monname}/${stnidlow}${int_mon}${Y_P}.txt.gz
  wget -q ${web_oce1} -O temp.gz
  actualsize=$(wc -c <"temp.gz")
  if [ $actualsize -gt 0 ]; then
    gunzip temp.gz
    cat temp >> $WGETOUT1
    if [ -s temp ]; then
      rm -f temp
    fi
  else
    rm -f temp.gz
  fi
  unset web_oce1
  
  web_oce1=${WEBoce1}/${monname}/${stnidlow}.txt
  wget -q ${web_oce1} -O temp
  actualsize=$(wc -c <"temp")
  if [ $actualsize -gt 0 ]; then
    cat temp >> $WGETOUT1
    if [ -s temp ]; then
      rm -f temp
    fi
  else
    rm -f temp
  fi
  unset web_oce1
done
####  End of adding by Zheng
 

####  Get recent data of present year
wget -q https://www.ndbc.noaa.gov/data/realtime2/$stnid.txt -O temp
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
####  End of adding by Zheng

#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft
#2011 11 17 13 30  50 10.0 13.0    MM    MM    MM  MM 1027.0    MM  20.5    MM   MM   MM    MM
#2011 11 17 13 00  50 10.0 14.0    MM    MM    MM  MM 1026.3    MM  20.5    MM   MM +4.2    MM
#2011 11 17 12 30  30 11.0 14.0    MM    MM    MM  MM 1025.6    MM  20.5    MM   MM   MM    MM

wget -q https://www.ndbc.noaa.gov/data/realtime2/$stnid.ocean -O temp
####  Begin of adding by Zheng
actualsize=$(wc -c < "temp")
if [ $actualsize -gt 0 ]; then
  tac temp > temp1
  cat temp1 >> $WGETOUT1
  if [ -s temp ]; then
    rm -f temp temp1
  fi
else
  rm -f temp
fi
####  End of adding by Zheng

# Now have a file with stuff like:
#YY  MM DD hh mm   DEPTH  OTMP   COND   SAL   O2% O2PPM  CLCON  TURB    PH    EH
#yr  mo dy hr mn       m  degC  mS/cm   psu     %   ppm   ug/l   FTU     -    mv
#2011 11 18 12 30     2.0 20.15     MM    MM    MM    MM     MM    MM    MM    MM
#2011 11 18 12 00     2.0 20.18     MM    MM    MM    MM     MM    MM    MM    MM
#2011 11 18 11 30     2.0 20.18     MM    MM    MM    MM     MM    MM    MM    MM

DATEF=`$HOME1/bin/dateformat $begindate "%Y%m%d%H"`
DATEL=`$HOME1/bin/dateformat $enddate "%Y%m%d%H"`
echo read_ndbc $sensor DATEF DATEL $DATEF $DATEL 
awk -v datef=$DATEF '$1$2$3$4 >= datef {print }'  $WGETOUT1 > $AWKED1
awk -v datel=$DATEL '$1$2$3$4 <= datel {print }'  $AWKED1 > $AWKED

#  For Temperature observation data
if [ $sensor == "WT" ]; then
   awk '$7 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$7)  }' \
   $AWKED | sort -u > $TMPF
   awk '$6 <= 90.0 { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$6)  }' \
   $TMPF | sort -u > $fileout
fi

#  For Conductivity observation data
if [ $sensor == "WC" ]; then
   awk '$8 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.2f \n",$1,$2,$3,$4,$5,$8)  }' \
   $AWKED | sort -u > $TMPF
   awk '$6 <= 999.0 { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$6)  }' \
   $TMPF | sort -u > $fileout
fi

#  For Salinity observation data
if [ $sensor == "SAL" ]; then
   awk '$9 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$9)  }' \
   $AWKED | sort -u > $TMPF
   awk '$6 <= 90.0 { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$6)  }' \
   $TMPF | sort -u > $fileout
fi

awk -v datef=$DATEF '$1$2$3$4 >= datef {print }'  $WGETOUT > $AWKED1
awk -v datel=$DATEL '$1$2$3$4 <= datel {print }'  $AWKED1 > $AWKED

# select sensor in $stnid.txt and $stnid.ocean   AT $14   WT $15   PRES $13 

if [ $sensor == "AT" ]; then
   awk '$14 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.2f \n",$1,$2,$3,$4,$5,$14)  }' \
   $AWKED | sort -u > $fileout
fi

if [ $sensor == "WT" ]; then
   awk '$15 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$15)  }' \
   $AWKED | sort -u > $TMPF
   awk '$6 <= 90.0 { printf("%04i %02i %02i %02i  %02i  %8.1f \n",$1,$2,$3,$4,$5,$6)  }' \
   $TMPF | sort -u >> $fileout
   sort -n $fileout > temp
   cat temp > $fileout
   rm -f temp
   sort -u $fileout > temp
   cat temp > $fileout
   rm -f temp
fi

if [ $sensor == "AP" ]; then
   awk '$13 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.2f \n",$1,$2,$3,$4,$5,$13)  }' \
   $AWKED | sort -u > $fileout
fi

if [ $sensor == "DEWP" ]; then
   awk '$16 != "MM" { printf("%04i %02i %02i %02i  %02i  %8.2f \n",$1,$2,$3,$4,$5,$16)  }' \
   $AWKED | sort -u > $fileout
fi

if [ -s $fileout ]; then
  head -1 $fileout
  tail -1 $fileout
fi

rm -f $WGETOUT
rm -f $WGETOUT1
rm -f $WGETLOG
rm -f $AWKED1
rm -f $AWKED
rm -f $TMPF
exit

