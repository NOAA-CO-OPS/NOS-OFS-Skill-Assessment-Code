#!/bin/sh
#Name:                get_TS_COOPS.sh 
#purpose:             gets CO-OPS surface water temperature and conductivity 
#                     within BEGINDATE and ENDDATE
# Author:             Aijun Zhang
# Date:               03/23/2006  
#Language:            Korn Shell Script
#input parameters:    BEGINDATE,ENDDATE,STATIONDATA
# modified from get_WL_verified.sh since CO-OPS switches to new web page
# Programs Called:
#          Name               Location                            Description
#        refwl.x               $BIN              FORTRAN program to reformat to a standard format
#source STEPS_SETUP.sh
# Modified by Lianyuan Zheng on 03/01/2017

cd $WRK_DIR
echo get_TS_COOPS.sh $BEGINDATE " to " $ENDDATE
DATATYPET='Meteorological+Observations'
DATATYPES='Conductivity'

if [ $DELT_O = 60 -o $DELT_O = 60.0 ]; then
  DELTTYPE='h'     # hourly data
elif [ $DELT_O = 6 -o $DELT_O = 6.0 ]; then
  DELTTYPE=''      # 6-minute data
fi

# get water temperature data because it is needed to calculate the salnity from the conductivity
BEGINDATE0=`$BIN/dateformat $BEGINDATE "%Y %m %d 00 00"`
ENDDATE0=`$BIN/datemath $ENDDATE + 0 0 2 0 0`
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel
do
   rm -f temp.raw cond.raw salt.raw
   rm -f temp.out cond.out salt.out tandcond.out
   echo " " StationNames $stnid ":" $stationname ":" $longlabel
   read Latlon
   tbegin=$BEGINDATE0
   tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
   while [ `$BIN/dateformat $tbeginp30 "%Y%m%d%H"` -lt `$BIN/dateformat $ENDDATE0 "%Y%m%d%H"` ]
   do
      TEMPLATE=$CTL/request.template_wt_verified_6min
      rm -f output.txt request.GET
      bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
      edate=`$BIN/dateformat $tbeginp30 "%Y%m%d"`

#      echo "getnwlon bdate edate "   $bdate $edate
      sed -e s/VSTNIDV/$stnid/ \
          -e s/VBDATEV/$bdate/ \
          -e s/VEDATEV/$edate/ $TEMPLATE > request.GET

      wget -o junk -O output.txt -i request.GET
      cat output.txt >> temp.raw

#  get water conductivity
      if [ $KINDAT -eq 4 ]; then
         rm -f output.txt request.GET
	 TEMPLATE=$CTL/request.template_wc_verified_6min
	 sed -e s/VSTNIDV/$stnid/ \
             -e s/VBDATEV/$bdate/ \
            -e s/VEDATEV/$edate/  $TEMPLATE > request.GET

	 wget -o junk -O output.txt -i request.GET
	 cat output.txt >> cond.raw
      fi
      tbegin=$tbeginp30
      tbeginp30=`$BIN/datemath $tbegin + 0 0 30 0 0`
   done

   rm -f output.txt request.GET
   TEMPLATE=$CTL/request.template_wt_verified_6min
   bdate=`$BIN/dateformat $tbegin "%Y%m%d"`
   edate=`$BIN/dateformat $ENDDATE0 "%Y%m%d"`
#   echo "getnwlon bdate edate "   $bdate $edate
   sed -e s/VSTNIDV/$stnid/ \
       -e s/VBDATEV/$bdate/ \
       -e s/VEDATEV/$edate/  $TEMPLATE > request.GET

   wget -o junk -O output.txt -i request.GET
   cat output.txt >> temp.raw
   awk -v sid=$stnid -v NCL=$NCOL ' $1 == sid && NF >= NCL \
     {print substr($4,1,4)" "substr($4,6,2)" "substr($4,9,2)" " \
     substr($5,1,2)" "substr($5,4,2)" "$6 }' temp.raw  | sort -u > tmp1.out
   
   if [ -s tmp1.out  ]; then
     cp tmp1.out $stationname"_temp.raw"
     cat tmp1.out | $BIN/fillnan.x "$BEGINDATE0" "$ENDDATE0" $DELT_O " -999.99" > temp.out
   fi

   if [ $KINDAT -eq 4 ]; then
      rm -f output.txt request.GET
      TEMPLATE=$CTL/request.template_wc_verified_6min
      sed -e s/VSTNIDV/$stnid/ \
          -e s/VBDATEV/$bdate/ \
          -e s/VEDATEV/$edate/  $TEMPLATE > request.GET

      wget -o junk -O output.txt -i request.GET
      cat output.txt >> cond.raw
      awk -v sid=$stnid -v NCL=$NCOL ' $1 == sid && NF >= NCL \
	{print substr($4,1,4)" "substr($4,6,2)" "substr($4,9,2)" " \
	substr($5,1,2)" "substr($5,4,2)" "$6 }' cond.raw  | sort -u > tmp2.out

      if [ -s tmp2.out  ]; then
         cp tmp2.out $stationname"_cond.raw"
         cat tmp2.out | $BIN/fillnan.x "$BEGINDATE0" "$ENDDATE0" $DELT_O " -999.99" > cond.out
         cat temp.out | awk '{printf "%10.5f \n",$6}' > tmp     
         paste cond.out tmp > tmp1
         awk ' $6 >= -99 && $7 >= -99 {print $1" "$2" "$3" "$4" "$5" 00 "$6" "$7 }' tmp1  > tandcond.out

#echo "SALINITY.pl written for conductivity from NWLON, USGS gives temperature"
         $HOME1/scripts/SALINITY.pl tandcond.out tmp2
         awk '{printf("%4d %02d %02d %02d %02d   %10.5f \n",$1,$2,$3,$4,$5,$7)}' tmp2 >> salt.out
      fi 
   fi

   if [ $KINDAT -eq 3 ]; then
      FILEIN='temp.out'
      if [ -s $stationname"_temp.raw"  ]; then
        cp -p $stationname"_temp.raw" $FILEIN
      else	
        if [ ! -s $FILEIN  ]; then
	  $HOME1/scripts/get_data_ndbc.sh $stnid "$BEGINDATE0" "$ENDDATE0" WT DAT
        fi
        if [ -s DAT ]; then
          mv DAT $FILEIN
        fi  
      fi	

   elif [ $KINDAT -eq 4 ]; then
      FILEIN='salt.out'
      if [ ! -s $FILEIN  ]; then
	$HOME1/scripts/get_data_ndbc.sh $stnid "$BEGINDATE0" "$ENDDATE0" SAL DAT
      fi
      if [ -s DAT ]; then
        mv DAT $FILEIN
      fi  
   fi

# For SFBOFS, get real time observations from California Data Exchange Center for T & S
   rm -f temperature.tmp salinity.tmp
   if [ $stnid == "9415316" ]; then
     CDEC=RIV
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415064" ]; then
     CDEC=ANC
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415176" ]; then
     CDEC=CLL
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415112" ]; then
     CDEC=MAL
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415144" ]; then
     CDEC=PCT
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415102" ]; then
     CDEC=MRZ
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "9415265" ]; then
     CDEC=GOD
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   elif [ $stnid == "SFB1211" ]; then
     CDEC=374
     $HOME1/scripts/get_TS_CDEC.sh $CDEC "$BEGINDATE0" "$ENDDATE0"
   fi

   if [ -s temperature.tmp ]; then
     cp temperature.tmp temp.out
   fi

   if [ -s salinity.tmp ]; then
     cp salinity.tmp salt.out
   fi

   $BIN/refwl.x "$BEGINDATE0" "$ENDDATE0" $FILEIN $stationname"_msl.6min" 
   if [ -s $stationname"_msl.6min" ]; then
#     echo    $stationname"_msl.6min" 
#     head -1 $stationname"_msl.6min"
#     tail -1 $stationname"_msl.6min"
     mv $stationname"_msl.6min" $OBS/$stationname".obs"
   else
     rm -f   $stationname"_msl.6min"
   fi  

  done 3<&-
  echo "done of get_TS_COOPS.sh!!!"
exit

