#!/bin/sh
set -x
###  get datum and tidal constituents, convert tide constituents from MSL to MLLW by adding the difference 
cd /home/net/azhang/SKILL_OFS/work
STATIONDATA='../control_files/cbofs_wl_station.ctl'
BIN='../bin'
BEGINDATE="2003 01 02 00 00"
ENDDATE="2003 12 30 00 00"
bdate=`$BIN/dateformat $BEGINDATE "%Y%m%d"`
edate=`$BIN/dateformat $ENDDATE "%Y%m%d"`
echo "getnwlon bdate edate "   $bdate  $edate
# Loop on lines in $stationdata
exec 5<&0 <$STATIONDATA
while read stnid stationname longlabel 
do 
  echo $stnid $stationname 
  read lat longitude XMAJOR
  if [ ! -s ${stnid}.cons ]; then

##   getting datum 
   sed -e s/VSTNIDV/$stnid/ \
       -e s/VBDATEV/$bdate/ \
       -e s/VEDATEV/$edate/  \
     ../control_files/request.template_datum_verified > request.GET_datum
    wget --no-check-certificate -o junk -O $stnid'.datum' -i request.GET_datum

#    awk -v -v NCL=MSL ' substr($1,1,3) == 'MSL' {print substr($1,1,3)" "substr($1,6,2)" "substr($1,9,2)" " \
#    tail -2 junk0 | head -1> tmp0
#    awk '{printf("%10d  %6.3f %6.3f %6.3f \n",$1,$6-$8,$6,$8)}' tmp0 > $stnid'.datum'
#    read sidd dif msl mllw < $stnid'.datum'
#    echo $sidd $dif $msl $mllw

##   get harmonic constants from CO-OPS

   sed -e s/VSTNIDV/$stnid/ \
       -e s/VBDATEV/$bdate/ \
       -e s/VEDATEV/$edate/  \
     ../control_files/request.template_HC_verified > request.GET_ha

    wget --no-check-certificate -o junk -O $stnid'.ha' -i request.GET_ha
    $BIN/reformat_ha.x $stnid
  fi
done 3<&-
exit

