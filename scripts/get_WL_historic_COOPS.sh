#!/bin/sh
#Name:                get_WL_historic.sh 
#purpose:             gets CO-OPS verified historic water levels at NWLON and Great Lakes Stations 
#                     within BEGINDATE and ENDDATE
# Author:             Aijun Zhang
# Date:               03/23/2006  
#Language:            Korn Shell Script
#input parameters:    BEGINDATE,ENDDATE,STATIONDATA
# modified from get_WL_verified.sh since CO-OPS switches to new web page
# Programs Called:
#          Name               Location                            Description
#        refwl.x               $BIN              FORTRAN program to reformat to a standard format
# source STEPS_SETUP.sh
#  modified by Aijun Zhang on Feb. 23, 2009
#  Data format on CO-OPS web site was changed
#  no tidal prediction is provided, so NCOL=4
#
# Modified by Lianyuan Zheng on 03/01/2017

  cd $WRK_DIR
  echo get_WL_historic.sh $BEGINDATE " to " $ENDDATE
  NCOL=5
  DATUM=0.0

#  Convert datum from IGLD 1985 to Chart Datum Low Water Datum  
  if [ $DBASE = "GLAKES" ]; then
    if [ $ofs = "leofs" ]; then
       DATUM=173.5
    elif [ $ofs = "lhofs" ]; then
       DATUM=176.0
    elif [ $ofs = "lmofs" ]; then
       DATUM=176.0
    elif [ $ofs = "loofs" ]; then
       DATUM=74.2
    elif [ $ofs = "lsofs" ]; then
       DATUM=183.2
    elif [ $ofs = "lmhofs" ]; then
       DATUM=176.0
    fi

    if [ $DELT_O = 60 -o $DELT_O = 60.0 ]; then
      TEMPLATE=$CTL/request.template_wl_verified_60min_GL
    elif [ $DELT_O = 6 -o $DELT_O = 6.0 ];  then
      TEMPLATE=$CTL/request.template_wl_verified_6min_GL
    fi
  else
    if [ $DELT_O = 60 -o $DELT_O = 60.0 ]; then
      TEMPLATE=$CTL/request.template_wl_verified_60min
    elif [ $DELT_O = 6 -o $DELT_O = 6.0 ]; then
      TEMPLATE=$CTL/request.template_wl_verified_6min
    fi
  fi

  WGETOUT=`mktemp -q wgetout.XXXXXX`
  BEGINDATE0=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
  ENDDATE0=`$BIN/datemath $ENDDATE + 0 0 2 0 0`

  exec 5<&0 <$STATIONDATA
  while read stnid stationname longlabel
  do
    rm -f $WGETOUT
    echo " " StationNames $stnid ":" $stationname ":" $longlabel
    read Latlon
    tbegin=$BEGINDATE0
    tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
    while [ `$BIN/dateformat $tbeginp30 "%Y%m%d%H"` -lt `$BIN/dateformat $ENDDATE0 "%Y%m%d%H"` ]
    do
      rm -f output.txt
        bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
        edate=`$BIN/dateformat $tbeginp30 "%Y%m%d"`
        sed -e s/VSTNIDV/$stnid/ \
            -e s/VBDATEV/$bdate/ \
            -e s/VEDATEV/$edate/ $TEMPLATE > request.GET

#  Begin modification by Zheng
#  Note: The month variables are not used here.  
#           -e s/BMONTH/$bmonth/ \
#           -e s/EMONTH/$emonth/ $TEMPLATE > request.GET
#  End of modification by Zheng

        wget --no-check-certificate -o junk -O output.txt -i request.GET
        cat output.txt >> $WGETOUT
        tbegin=$tbeginp30
        tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
     done   

     bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
     edate=`$BIN/dateformat $ENDDATE0 "%Y%m%d"`
     sed -e s/VSTNIDV/$stnid/ \
         -e s/VBDATEV/$bdate/ \
         -e s/VEDATEV/$edate/ $TEMPLATE > request.GET
     wget --no-check-certificate -o junk -O output.txt -i request.GET
     cat output.txt >> $WGETOUT

     if [ -s $WGETOUT ]; then
      awk -v sid=$stnid -v NCL=$NCOL -v datum=$DATUM ' $1 == sid && NF >= NCL \
      {print substr($2,1,4)" "substr($2,6,2)" "substr($2,9,2)" " \
       substr($3,1,2)" "substr($3,4,2)" "$4-datum }' $WGETOUT > tmp.out
      cp tmp.out $WGETOUT
     fi

     $BIN/refwl.x "$BEGINDATE0" "$ENDDATE0" $WGETOUT $stationname"_msl.6min" 
     if [ -s $stationname"_msl.6min" ]; then
       head -1 $stationname"_msl.6min"
       tail -1 $stationname"_msl.6min"
       mv $stationname"_msl.6min" $OBS/$stationname".obs"
     else
       rm -f $stationname"_msl.6min"
     fi
     rm -f junk $WGETOUT
  done 3<&-
  echo "done of get_WL_historic_COOPS.sh!!!"
exit

