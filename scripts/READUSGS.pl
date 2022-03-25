#!/usr/bin/perl
#
# Documentation for Scripts READUSGS.pl
#
#---------------------------------------------------------------------------------------
#
# Script Name: READUSGS.pl
#
# Directory Location:    /COMF/oqcs/scripts 
# 
# Technical Contact(s): Name:  Tom Gross            Org:    NOS/CSDL
#                       Phone: 301-713-2809x139     E-Mail: tom.gross@noaa.gov
#
# Abstract:
#     		This Perl script parses the USGS web page 
#		(grabbed using river_read_usgs.sh) to isolate and return
#		one or more of these data variables:
#
#		 TEMP	       00010   - TEMPERATURE, WATER (DEG. C)
#		 COND	       00095   - SPECIFIC CONDUCTANCE (MICROSIEMENS/CM AT 25 DEG. C)
#		 DISCHARGE     00060   - DISCHARGE, CUBIC m PER SECOND
#		 GAGE	       00065   - GAGE HEIGHT, m
#		 Returns  -9999.00000 if the station doesn't have the data type                       
#		 Converstion from feet to MKS done here
#
# Usage:  Interactively:     READUSGS.pl $INPUTFILENAME "variable list" $OUTPUTFILENAME
#			     READUSGS.pl WGETOUT 'DISCHARGE GAGE' rivout.dat
#              Via cron:     Called by river_read_usgs.sh 
#
# Input Parameters:     $inputfilename : input file name 
#			The variable list should be something like:
#			"DISCHARGE"
#			"TEMP"  
#			"TEMP COND"
#			"TEMP DISCHARGE GAGE COND"
#			"DISCHARGE GAGE TEMP COND"
#			$outputfilename : output file name
#
# Language:  Bourne Shell Script
#
# Target Computer:  COMF computer, such as dsofs1.nos-tcn.noaa.gov   
#
# Estimated Execution Time: 
#
# Scripts/Programs Called:
#            Name               Directory Location      Description
# 
#
# Input Files:
#            Name               Directory  Location     Description
#	$INPUTFILENAME		User defined		a web page
#		
# Output Files:
#            Name               Directory Location      Description
#	$OUTPUTFILENAME		User defined		simple ascii file(with only as many v2,v3,v4 as requested)
#							y m d h m 00 v1 v2 v3 v4
#
# Error Conditions:
#
# Author Name:  Tom Gross                   Creation Date:  
#
#
# Revisions:
#         Date          Author          Description
#  	2005-02-04	H Lin		Learn Perl here and correct some error messages.
#
#
#
# Remarks: 
#	This is a real tour de force in hash table redirection sub-
#	scripting.  Study the line:
#	printf(OUTFILE " %f",$VAL[$POS{$VARDD{$VARIABLELIST[$i]}}]
#                            *$MKS{$VARIABLELIST[$i]});
#
# -----------------------------------------------------------------------------------------------------------------




$INPUTFILENAME=$ARGV[0];
@VARIABLELIST=split ' ', $ARGV[1];
$OUTPUTFILENAME=$ARGV[2]; 

open(OUTFILE, "> $OUTPUTFILENAME")
    || die "Cannot open".$OUTPUTFILENAME." for output\n";
    
open(INFILE, "< $INPUTFILENAME")
    || die "Cannot open".$INPUTFILENAME." for input\n";

# hash of possible variables by name and usgs parameter
# TEMP   00010   - TEMPERATURE, WATER (DEG. C)
# COND  00095   - SPECIFIC CONDUCTANCE (MICROSIEMENS/CM AT 25 DEG. C)
# DISCHARGE   00060   - DISCHARGE, CUBIC FEET PER SECOND
# GAGE   00065   - GAGE HEIGHT, FEET

%VARDD=('DISCHARGE','00060','GAGE','00065','COND','00095','TEMP','00010');
%POS=('00060',0,'00065',0,'00095',0,'00010',0);

# conversion factors  ft3/s>m3/s  ft>m  MICROSIEMENS/CM > mS/cm ,  DEG. C
%MKS=('DISCHARGE',0.028317,'GAGE',0.3048,'COND',0.001,'TEMP',1.0);

#print "variablelist @TEST>",@VARIABLELIST,"\n" ;
#print "variablelist \$TEST>",$VARIABLELIST[0],"\n" ;
#print "HASH TEST>",$VARDD{GAGE},"\n" ;
#print "HASH TEST2>",$VARDD{$VARIABLELIST[0]},"\n" ;

$nv=scalar(@VARIABLELIST);
print "length VARIABLELIST>",$nv,"\n";

    
while (<INFILE>) {
#@VAL=split(/\s+/);
@VAL=split(/[\s+_]/);

if ($VAL[0] eq agency )
{
print ;
$num=split(/[\s+_]/);
print $num,"=num\n" ;
#printf(OUTFILE " if  %s %s %s\n",$VAL[2],$VAL[6],$VAL[8]);
#printf( " if  %s %s %s\n",$VAL[2],$VAL[6],$VAL[8]);
# process the line to find what variables are in what columns
# n points to positions in this line, getting the data descriptors.
# agency_cd	site_no	datetime	02_00065	01_00060
# num objects returned above.  The data are in 
#  n= num:-2:7     -1   (0:num-1)    6,8,10,12,14,num-1
#  i points to data position in data line
# USGS	01570500	2003-04-17 00:00	6.39	63,700
#  i=7:8    (counts start at 0,  i++ before usage)
$i=6;
for ( $n=6; $n<=$num ; $n=$n+2 ) 
  { $i=$i+1;
    $POS{$VAL[$n]}=$i;
#print "i,n,val ",$i,' ',$n,' ',$VAL[$n],"\n";
   }
#print "Hash POS@>",%POS,"\n" ;


}
elsif ($VAL[0] eq USGS)
 { 
 s/,//;
 $_=~s/(\d)[-:]/$1 /g;
  #printf(STDOUT "%s\n", $_);

#@VAL=split(/[\s-:]/);
@VAL=split(/\s/);

#print "USGS Equal\n" ; 
#print ;
#print $VAL[7];
if ($VAL[7]>0.000) {
  printf(OUTFILE "%d %02d %02d %02d %02d  0 ",
  $VAL[2],$VAL[3],$VAL[4],$VAL[5],$VAL[6]);
  $VAL[0]=-9999.;
  for ($i=0; $i<$nv; $i++) {
     printf(OUTFILE " %f",$VAL[$POS{$VARDD{$VARIABLELIST[$i]}}]
                            *$MKS{$VARIABLELIST[$i]});
    #printf(" V%i",$POS{$hash{$VARIABLELIST[$i]}});
   }
  printf(OUTFILE " \n");
}

}
}



