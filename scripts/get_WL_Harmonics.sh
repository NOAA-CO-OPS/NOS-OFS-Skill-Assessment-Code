#!/bin/sh
source ${HOME1}/scripts/STEPS_SETUP.sh
N=0
CTL='/home/net/azhang/SKILL_OFS/control_files'
#STATIONDATA='Station_Index_COOPS.txt'
bdate='20100101'
edate='20100707'
exec 5<&0 < $STATIONDATA
  while read stnid stationname longlabel 
  do
   echo 'stnid= ' $stnid
   read lat longitude XMAJOR
##  get Datum 
    sed -e s/VSTNIDV/$stnid/ \
        -e s/VBDATEV/$bdate/ \
        -e s/VEDATEV/$edate/  \
        $CTL/request.template_datum > request.GET_datum

    wget --no-check-certificate -o junk -O $stnid'.datum' -i request.GET_datum

##   get harmonic constants from CO-OPS
    sed -e s/VSTNIDV/$stnid/ \
        -e s/VBDATEV/$bdate/ \
        -e s/VEDATEV/$edate/  \
        $CTL/request.template_ha > request.GET_ha

    wget --no-check-certificate -o junk -O $stnid'.ha' -i request.GET_ha
    $BIN/reformat_ha.x  $stnid
    if [ -s $stnid'.std' ]
    then
          mv $stnid'.std' ${HOME1}/data/harmonic_con/$stnid'.std'
    fi
  done
done 3<&-
exit
