// Documentation for C language program : dateformat.c
//
//-------------------------------------------------------------------------------------------------------------------
//
// C Program Name: dateformat.c
//
// Directory Location:   /COMF/oqcs/sorc  
//
// Technical Contact(s):   Name: Tom Gross		Org:  NOS/CSDL
//			   Phone: 301-713-2809x139  	E-Mail: tom.gross@noaa.gov
//
// Abstract:
//		To create arbitrary strings using date components such
//        	as year, month, day, julian day, etc...	   
//		dateformat.c is compiled to the executable dateformat.
//		Input year month day hour min format and it outputs date string
//
// Usage:  dateformat year month day hour min format
//         year month day hour min are strings for the date.
//         format is a format string using the same syntax as
//         the format strings of the UNIX command  date +format
//         A string output is created with the specified date and format
//
//	   Example Usage
//		Command line usage example:
//		>> dateformat 1998 3 13 13 30 "tdl%Y%m/cbbt/%Y%j%H%M.out" 
//		>> tdl199803/cbbt/19980721330.out
//
//	Usage inside a script:
//        	str=`dateformat 1998 3 13 13 30 "tdl%Y%m/cbbt/%Y%j%H%M.out"`
//        	echo the string is:$str
//
//        	It can be useful to put the dateformat, datemath routines into
//       	your path.  Add this line near the beginning of your script:
//        	put dateformat datemath in the path:
//        	PATH=$PATH:/COMF/oqcs/sorc/binlinux
//
//	Notes:
//	     	Years must be specified as four digit years  1998, 2001 not 98 or 01
//
//      	So called julian year days can be used as input if they are
//        	assumed to be days of the month of January.  That is Mar. 4
//        	is yearday 63 and can be specified as either:
//        	>>dateformat 1998 3  4 12 30 "%Y %m %d %j "
//		>>1998 03 04 063
//       	>>dateformat 1998 1 63 12 30 "%Y %m %d %j "
//		>>1998 03 04 063
//
//		Time starts with Jan 01. When you want to know 90th of the day,  
//		Using dateformatr 2005 01 90.43 0 0 "%Y %m %d %H %M %S"
//		NOT dateformatr 2005 01 89.43 0 0 "%Y %m %d %H %M %S"
//
//	Bugs:   Dates are only valid in years following 1900.
//        	Therefore do not use
//        	dateformat 0 0 3 12 30 " %d %H %m "
//        	Instead put in a dummy year and month:
//       	dateformat 1998 1 3 12 30 " %d %H %m "
//
//	General Purpose of Script Usage Example:
//
//	A directory and file system has been created with date information.
//	We are looking for the file whose name specifies a date and time at least
//	twelve hours before the present time but not more than 50 hours before
//	(bailout puts the 50 hour limit.  Without a bailout variable this can
//	become an endless loop.):
//
//	# Read the realtime UNIX clock
//	nowdate=`date +"%Y %m %d %H %M"`
//	# Start with 12 hour old date by decrementing the date by 12 hours
//	nowdate=`datemath $nowdate - 0 0 0 12 0`
//	# Construct a file name e.g.  /home/datedir/199807/cbbt/199820512.out
//	filename=`dateformat $nowdate "/home/datedir/%Y%m/cbbt/%Y%j%H.out"`
//	bailout=0
//	echo count $bailout FILE= $filename
//	# Loop until the filename exists or bailout exceeds 48
//	while [ ! -s "$filename" ] && [ $bailout -lt 48 ]
//	 do
//	# Decrement the date by one hour
// 	 nowdate=`eval datemath $nowdate - 0 0 0 1 0`
//	  filename=`eval dateformat $nowdate "/home/datedir/%Y%m/cbbt/%Y%j%H.out"`
//	  bailout=`expr $bailout  + 1`
// 	 echo "test" $bailout "exist ? FILE=" $filename
//	 done
//	echo $filename
//
//
//	Files:  dateformat, datemath  the executables
//        	dateformat.c datemath.c  the C source code (cc -lm ....)
//        	dateformat.txt  Text files describing the usage
//        	testdates.sh   A script to test the usage.
//
// Input Parameters:  	year		must be specified as four digit years  1998, 2001 not 98 or 01
//			month		digits from 1 to 12
//			day		digits from 1 to 31
//			hour		digits from 1 to 24
//			min		digits from 0 to 59
//			format 		 	
//
// Language:	C	  
//
// Compiling/Linking Syntax: cc dateformat.c -o dateformat.x -lm
//
// Target Computer:  Runs on COMF computers, such as dsofs1.nos-tcn.noaa.gov  
//
// Estimated Execution Time: 
//
// Suboutines/Functions Called:
//	Name			Directory Location	       Description
// 
//
// Input Files:
//  Unit No.	   Name 		Directory Location	        Description
//	none
//
// Output Files:
//  Unit No.	    Name		Directory Location		Description
//	none
//
// Libraries Used:     none
//
// Error Conditions:
//
// Author Name: Tom Gross		Creation Date:  Jan.  1998
//
// Revisions:
//	   Date 		Author 		Description
//	08-23-2004		H Lin		Standardized the documentation.
//	03-07-2005		H Lin		Test and clean the code. 
//
//
// Remarks: 
//
//
// -------------------------------------------------------------------------------------------------------------


#include <stdlib.h>
#include <time.h>
#include <stdio.h>

int main(argc,argv)
int argc;
char **argv;
{
	struct tm now;
	time_t clock;
	char str[120];
        char format[120];
	
        clock = time(NULL);

	/*  All dates are internally GMT and displayed as GMT
    	This avoids conversions for daylight savings time*/
	putenv("TZ=GMT0");

	if(argc==7){

		now.tm_year = atoi(argv[1]) - 1900;
  		if(now.tm_year <0) now.tm_year+=1900;
  		if(atoi(argv[1])==0) now.tm_year = 50;
		now.tm_mon  = atoi(argv[2]) - 1;
  		if(atoi(argv[2])==0) now.tm_mon = 0;
		now.tm_mday = atoi(argv[3]);
		now.tm_hour = atoi(argv[4]);  
		now.tm_min  = atoi(argv[5]);
		now.tm_sec  = 0;
		/* now.tm_isdst = 0; */
		clock = mktime(&now);
		strncpy(format,argv[6],120);
		
        	strftime(str,120,format,&now);
        	printf(str);
	}
	else{ 
		printf("Need_6_arguments_year_month_day_hr_min_format\n");
	}
}
	


