#!/bin/sh
#Name:                get_currents_COOPS.sh 
#purpose:             gets Historic 6-minutes  current, water Temp, and Salt at PORTS stations
#                          within BEGINDATE and ENDDATE
# Author:             Aijun Zhang
# Date:               11/20/2004  
#Language:            Korn Shell Script
#
#input parameters:    BEGINDATE,ENDDATE,STATIONDATA
# Programs Called:
#          Name       Location                            Description
#        refwl.x      $BIN              FORTRAN program to reformat to a standard format

#         <SELECT mtype>
#            mtype=7  Air Temperature
#            mtype=8  Barometric Pressure
#            mtype=5  Salinity/Gravity
#            mtype=4  Water Currents
#            mtype=9  Water Level
#            mtype=10 Water Temperature
#            mtype=6  Winds
#         </SELECT>
#  cd $WRK_DIR
#  set -x
# Modified:           Lianyuan Zheng on 03/01/2017

  cd $WRK_DIR
  NCOL=4
  echo get_currents_COOPS.sh $BEGINDATE " to " $ENDDATE
  BEGINDATE0=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
  ENDDATE0=`$BIN/datemath $ENDDATE + 0 0 2 0 0`
#  echo $BEGINDATE0 " to " $ENDDATE0
  TEMPLATE=$CTL/request.template_cu_verified_6min
 
  exec 5<&0 <$STATIONDATA
  while read stnid stationname longlabel
  do
    echo " " StationNames $stnid ":" $stationname ":" $longlabel
    read Latlon

# try NDBC first
    if [ -s tmp.out ]; then
      rm -f tmp.out
    fi

    $HOME1/scripts/get_data_ndbc_currents.sh  $stnid "$BEGINDATE0" "$ENDDATE0" tmp.out 
    if [ -s tmp.out ]; then
      mv tmp.out $stationname"_curr"
    fi

    tbegin=$BEGINDATE0
    tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
    WGETOUT=`mktemp -q wgetout.XXXXXX`
    yearb=`$BIN/dateformat $tbegin "%Y"`

#  Note: In the opendao.co-ops.nos.noaa.gov current 6-min data, the requested length of time 
#        must be less than 31 days.  Following is getting 30 days 6-min current data.
    while [ `$BIN/dateformat $tbeginp30 "%Y%m%d%H"` -lt `$BIN/dateformat $ENDDATE0 "%Y%m%d%H"` ]
    do
      if [ -s output.txt ]; then
        rm -f output.txt
      fi
      bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
      edate=`$BIN/dateformat $tbeginp30 "%Y%m%d"`
#        echo "getnwlon bdate edate "   $bdate $edate 
 
      sed -e s/VSTNIDV/$stnid/ \
          -e s/VBDATEV/$bdate/ \
          -e s/VEDATEV/$edate/  $TEMPLATE > request.GET

      wget --no-check-certificate -o junk -O output.txt -i request.GET
      cat output.txt >> $WGETOUT 
      tbegin=$tbeginp30
      tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
    done

    bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
    edate=`$BIN/dateformat $ENDDATE0 "%Y%m%d"`
#     echo "getnwlon bdate edate "   $bdate $edate 

    sed -e s/VSTNIDV/$stnid/ \
        -e s/VBDATEV/$bdate/ \
        -e s/VEDATEV/$edate/  $TEMPLATE > request.GET
    wget --no-check-certificate -o junk -O output.txt -i request.GET
    cat output.txt >> $WGETOUT

    if [ -s $WGETOUT ]; then
      awk -v year0=$yearb -v NCL=$NCOL ' substr($1,5,1) == "-" && substr($1,8,1) == "-" && NF >= NCL \
        {print substr($1,1,4)" "substr($1,6,2)" "substr($1,9,2)" " \
        substr($2,1,2)" "substr($2,4,2)" "$3 " "$4 }' $WGETOUT > tmp
      cp tmp $WGETOUT
    fi

    if [ -s $stationname"_curr" ]; then
      cp -p $stationname"_curr" junk0
    fi

    if [ -s $WGETOUT ]; then
      cat $WGETOUT >> junk0
    fi  

#  Note:  The current speed unit is knot, rather than m/s.
    $BIN/reformat_PORTS.x "$BEGINDATE0" "$ENDDATE0" $KINDAT junk0 $stationname"_msl.6min" 
    if [ -s $stationname"_msl.6min" ]; then
#        echo    $stationname"_msl.6min" 
      mv $stationname"_msl.6min" $OBS/$stationname".obs"
    else
      rm -f $stationname"_msl.6min"
    fi
    rm -f junk junk0 tmp output.txt 
    rm -f $WGETOUT
  done 3<&-
  echo "done of get_currents_COOPS.sh!!!"
exit

