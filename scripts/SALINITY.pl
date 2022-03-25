#!/usr/bin/perl
#
# Documentation for Scripts SALINITY.pl  filein  fileout
#
#---------------------------------------------------------------------------------------
#
# Script Name: SALINITY.pl
#
# Directory Location:   /COMF/oqcs/scripts  
# 
# Technical Contact(s): Name:  Tom Gross            Org:    NOS/CSDL
#                       Phone: 301-713-2809x139     E-Mail: tom.gross@noaa.gov
#
# Abstract:
#                        
#  	  A perl script used to convert a file containing water conductivity
#	  and temperature into just Salinity.  
# 	  Converts the Conductivity and Temperature in input file into Salinity
# 	  Assumes the Pressure = 0.00
#
# 	  Uses snippets of code from the matlab toolbox "seawater"
#	  http://www.marine.csiro.au/~morgan/seawater/
#	  It is based on the Unesco
#	  equations:  http://www.ices.dk/ocean/procedures/standard_seawater.htm
#
#	  Check against this handy web based salinity calculator
#	  http://ioc.unesco.org/oceanteacher/resourcekit/M3/Converters/SeaWaterEquationOfState/
#		Sea%20Water%20Equation%20of%20State%20Calculator.htm
#
#	  It expects the conductivity to be in milli Siemens/cm.  Those
#	  are the units which are about 1/10 the ppt.  ie the conductivity
#	  at 35ppt, 15C = 4.2914 milli Siemens/cm
#
#	  Pressure is assumed to be 0, Sea Surface data only.  However the pressure
#	  correction is slight until you get past 100m or more.  
#	  
#	  The conductivity is temperature corrected by this routine.  So it
#         requires the water temperture.  
#
#	  This is incompatible with the USGS conductivities which are specific,
#	  meaning that they are converted to conductivity at 15C. 
#
#
# Usage:  Interactively:     SALINITY.pl "$CDAT" "$CDAT2"
#              Via cron:     Called by NetCDFgetstation_nwlon_fast.sh
#                            Called by SALTQCF.sh
#
# Input Parameters:     
#          "$CDAT" : Input file, has           
#		     y m d h m fh  conductivity  temperature
#          "$CDAT2" : Output file, has 
#		     y m d h m fh  salinity
#
# Language:  Bourne Shell Script
#
# Target Computer:   Runs on COMF computers, such as dsofs1.nos-tcn.noaa.gov. 
#
# Estimated Execution Time: 
#
# Scripts/Programs Called:
#       Name            Directory Location                  Description
#       none
#
# Input Files:
#       Name            Directory  Location                 Description
#      $CDAT		user defined			y m d h m fh  conductivity  temperature
#
# Output Files:
#       Name            Directory Location                  Description
#      "$CDAT2"		user defined			y m d h m fh  salinity
#
#
# Error Conditions:
#
# Author Name:  Tom Gross                   Creation Date:  
#
#
# Revisions:
#         Date                  Author             Description
#	08-19-2004		H Lin 		Standardized the documentation
#	03-28-2005    		H Lin		Update document.
#
# Remarks: 
#
# -----------------------------------------------------------------------------------------------------------------





{

# conductivity 35,15   milli Siemens/cm  
$c3515 = 4.2914;

# the other bad units, but what nwlon is stuck on....
$c3515 = 42.914;

# Eqn (3) p.7 Unesco.
$c0 =  0.6766097;
$c1 =  2.00564e-2;
$c2 =  1.104259e-4;
$c3 = -6.9698e-7;
$c4 =  1.0031e-9;

$a0 =  0.0080;
$a1 = -0.1692;
$a2 = 25.3851;
$a3 = 14.0941;
$a4 = -7.0261;
$a5 =  2.7081;

$b0 =  0.0005;
$b1 = -0.0056;
$b2 = -0.0066;
$b3 = -0.0375;
$b4 =  0.0636;
$b5 = -0.0144;

$k  =  0.0162;


$INPUTFILENAME=$ARGV[0];
$OUTPUTFILENAME=$ARGV[1];

open(OUTFILE, ">".$OUTPUTFILENAME) || die "Cannot open ".$OUTPUTFILENAME." for input \n";
open(INFILE, "<".$INPUTFILENAME) || die "Cannot open ".$INPUTFILENAME." for input \n";

while (<INFILE>) {

	@VAL=split(/\s+/);
	$COND=$VAL[6];
	$T=$VAL[7];

	if ($COND<0 || $T<0) {
		$SALT=-99999.000000
	} 
	else {

		# cndr =  C(S,T,P)/C(35,15,0)
		$cndr=$COND/$c3515;

		# rt = rt(T) = C(35,T,0)/C(35,15,0)
		# Eqn (3) p.7 Unesco.       (checks out with web page calculator)
		$rt = $c0 + ($c1 + ($c2 + ($c3 + $c4*$T)*$T)*$T)*$T;

		# Rt = Rt(S,T) = C(S,T,0)/C(35,T,0)     (checks out with web page calculator)
		# Rt= C(S,T,P)/C(35,15,0)  / ( C(35,T,0)/C(35,15,0) )
		$Rt= $cndr/$rt;

		# salinity = sw_sals(Rt,T) 
		$Rtx   = sqrt($Rt);
		$del_T = $T - 15.;
		$del_S = ($del_T / (1.+$k*$del_T) ) * 
        		( $b0 + ($b1 + ($b2+ ($b3 + ($b4 + $b5*$Rtx)*$Rtx)*$Rtx)*$Rtx)*$Rtx);
	
		$S = $a0 + ($a1 + ($a2 + ($a3 + ($a4 + $a5*$Rtx)*$Rtx)*$Rtx)*$Rtx)*$Rtx;

		$SALT = $S + $del_S;

		#printf(OUTFILE "%d %d %d %d %d %d s=%f c=%f t=%f rt=%f Rtx=%f delS=%f\n",
		#  $VAL[0],$VAL[1],$VAL[2],$VAL[3],$VAL[4],$VAL[5],$SALT ,$COND, $T, $rt, $Rtx,$del_S); 
	
	}# end of if negative 

	printf(OUTFILE "%d %d %d %d %d %d %f \n",
 	 $VAL[0],$VAL[1],$VAL[2],$VAL[3],$VAL[4],$VAL[5],$SALT ); 

   }
}

