// Documentation for C program datemath.c
//
//---------------------------------------------------------------------------------------
//
// C program Name: datemath.c
//
// Directory Location:  /COMF/oqcs/sorc   
// 
// Technical Contact(s): Name:  Tom Gross            Org:    NOS/CSDL
//                       Phone: 301-713-2809x139     E-Mail: tom.gross@noaa.gov
//
// Abstract:
//              Perform math on dates.  Days hours can be added or subtracted
//        from a date.  The program correctly keeps track of the roll over
//        of hours to days, days to month etc...
//        Differences of two full dates will give number of days and hours min
//        separating them.           
//  		Input two dates separated by "+" or "-" 
// 	 The dates are converted to julian, the math worked and the result
//  	 is returned as a string of year month day hour min 
//  	 Usually you want to add day hour min to a full date (not two dates)
//  	 No error checking is done for adding two dates.
//
// Usage: 
//	 datemath y1 mon1 d1 h1 min1  +[-] y2 mon2 d2 h2 min2
//        y mon d h min are strings for the date.
//        + or - is specified
//        Output is a string of
//        y mon d h m which is designed to be given to dateformat
//        or individual items extracted with
//
//	Example Usage
//	Command Line Usage Example:
//        	>> datemath 1998 3 2 12 30 - 0 0 3 20 0
//        	1998 02 26 16 30
//        	>> dateformat `datemath 1998 3 2 12 30 - 0 0 3 20 0` "%Y %h %d %H:%M"
//        	1998 Feb 26 16:30
//          	(Note the clever use of back quotes to make that work!)
//
//	Script Usage Example:
//		str1="1998 3 2 12 30"
//		str2="0 0 3 20 0"
//		strp=`datemath $str1 - $str2`
//		strpdate=dateformat $strp "%Y %h %d %H:%M"
//		echo "$str1  -  $str2  =  $strpdate"
//
//	When using "yeardays"  one may think of them as days in the
//	month of Jan.  
//	Thus year day 180 can be written as 1998 1 180 12 30  
//
// Compilation:
//	    cc datemath.c -o datemath -lm
//
// Input Parameters:     
//          See above Example Usage.         
//                       
// Language:  C language.
//
// Target Computer:   Runs on COMF computers, such as dsofs1.nos-tcn.noaa.gov  
//
// Estimated Execution Time: 
//
// Scripts/Programs Called:
//      Name                    Directory Location                  Description
// 	dateformat   		/COMF/oqcs/binlinux		the executables
//      dateformat.c		/COMF/oqcs/sorc			the C source code (cc -lm ....)
//
// Input Files:
//      Name               	Directory  Location                  Description
//	none
//
// Output Files:
//      Name               	Directory Location                   Description
//	None
//
// Error Conditions:
//
// Author Name:  Tom Gross                   Creation Date:  July. 1998
//
// Revisions:
//         Date                 Author              Description
//	July.  1998 		Tom Gross	    Original C coding
//	Oct. 22, 1998		Tom Gross	    Bug fix for multi-year date differences
//	July 2002		Tom Gross	    Linux, sgi tested
//	July 2003		Tom Gross	    round off of 1 sec added to j1 and jout
//	08-23-2004		H Lin		    Standard documentation.
//	03-07-2005		H Lin		    Test and clean the code.
//
// Remarks: 
//        It can be useful to put the dateformat, datemath routines into
//        your path.  Add this line near the beginning of your script:
//        # put dateformat datemath in the path:
//        PATH=$PATH:/COMF/oqcs/binliunx
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
//	   # Decrement the date by one hour
//	   nowdate=`eval datemath $nowdate - 0 0 0 1 0`
//	   filename=`eval dateformat $nowdate "/home/datedir/%Y%m/cbbt/%Y%j%H.out"`
//	   bailout=`expr $bailout  + 1`
//	   echo "test" $bailout "exist ? FILE=" $filename
//	 done
//	echo $filename
//
// -----------------------------------------------------------------------------------------------------------------





#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


void julian2tm(double *x,struct tm *date);
void tm2julian(struct tm *date, double *x);

int main(argc,argv)
int argc;
char **argv;
{
	struct tm date1,dateout;
	struct tm date2;
        double j1,j2,jout;
	char str[80];
        char plusminus[1] ;
	int day,hour,min ;

	if(argc==12){
		date1.tm_year = atoi(argv[1]) - 1900;
		date1.tm_mon  = atoi(argv[2]) - 1;
		date1.tm_mday = atoi(argv[3]);
		date1.tm_hour = atoi(argv[4]);
		date1.tm_min  = atoi(argv[5]);
		date1.tm_sec  = 0;
		date1.tm_isdst = 0;
		strcpy(plusminus,argv[6]);
		date2.tm_year = atoi(argv[7]) - 1900;
		date2.tm_mon  = atoi(argv[8]) - 1;
		date2.tm_mday = atoi(argv[9]);
		date2.tm_hour = atoi(argv[10]);
		date2.tm_min  = atoi(argv[11]);
		date2.tm_sec  = 0;
		date2.tm_isdst = 0;
	}
	else { 
		printf("datemath_needs_11_arguments__y_m_d_h_min_+-_y_m_d_h_min_\n");
		return;
	}

 
	tm2julian(&date1,&j1);
	j1=j1+0.000001;

	tm2julian(&date2,&j2);


	if (0==strcmp( plusminus,"+")) { 
		jout = j1 + j2; 
	} 
	else { 
		jout = j1 - j2; 
	}
	jout=jout+0.000001;


	julian2tm(&jout,&dateout);

/*  	debugging tools...
	printf(" date1 %s  <%i>  \n",argv[1],date1.tm_year);
	printf(" date1 %s  <%s>  \n",argv[6],plusminus);
	printf(" date1 %s  <%s>  \n",argv[6],plusminus);
	printf(" date2 %s \n",argv[7]);
	printf(" date2 %s  <%i>  \n",argv[7],date2.tm_year);
	printf(" through %s  <%s>  %s \n",argv[1],argv[6],argv[7]);
        strftime(str,100,"%Y %m %d %H %M",&date1);
        printf("\nfirst date1 %18.8f , %s \n",j1,str);
        printf("tm_year  %d \n",date1.tm_year);
        printf("tm_mon   %d \n",date1.tm_mon);
        printf("tm_mday  %d \n",date1.tm_mday);
        printf("tm_hour  %d \n",date1.tm_hour);
        printf("tm_min   %d \n",date1.tm_min);
         
        strftime(str,100,"%Y %m %d %H %M",&date2);
        printf(" second date2  %18.8f  , %s\n",j2,str);  
        printf("tm_year  %d \n",date2.tm_year);
        printf("tm_mon   %d \n",date2.tm_mon);
        printf("tm_mday  %d \n",date2.tm_mday);
        printf("tm_hour  %d \n",date2.tm_hour);
        printf("tm_min   %d \n",date2.tm_min);
         
        strftime(str,100,"%Y %m %d %H %M",&dateout);
        printf(" dateout final julian %18.8f , %s \n",jout,str);  
        printf("tm_year  %d   \n",dateout.tm_year);
        printf("tm_mon   %d   \n",dateout.tm_mon);
        printf("tm_mday  %d   \n",dateout.tm_mday);
        printf("tm_hour  %d   \n",dateout.tm_hour);
        printf("tm_min   %d   \n",dateout.tm_min);
*/

	putenv("TZ=GMT0");
	strftime(str,100,"%Y %m %d %H %M",&dateout);

	/*  If julian difference produced a date with small year, then we were going
    	for days between two dates.  So be more explicit than using the calender
    	generation routines  */
    	if ( dateout.tm_year < 1)  { 
		day = jout;
		hour = jout*24. - day*24;
		min = jout*24.*60. - day*24*60 - hour*60;
	 	sprintf(str,"0000 00 %5d %3d %3d"
	                       ,day
                               ,hour
	                       ,min);
        }
        printf(str);

}
	
 
void tm2julian(struct tm *datep, double *jdayp)
{
	char str[80];
	long y,m;
	double jday;
	struct tm date;

	putenv("TZ=GMT0");

	date = *datep;

	if ( date.tm_year == -1900 ) { 	
		y = 0; 
		m = 0;
		jday = date.tm_mday+( date.tm_hour + date.tm_min/60.)/24.;
 	}
	else {     
		if (1+date.tm_mon <= 2) { 
			y = date.tm_year -1; 
			m = date.tm_mon +1 +12;
		}
      		else { 
			y = date.tm_year; 
			m = date.tm_mon +1; 
		}
      		y = 365.25 * (y+1900) ;
      		m = 30.6001*(m+1);
     		jday = y + m + date.tm_mday + 1720981.5 + 
         		( date.tm_hour + date.tm_min/60.)/24.;
	}

	strftime(str,100,"tm2julian %Y %m %d %H %M",&date);
	*jdayp = jday;

}

void julian2tm(double *jdayp, struct tm *datep)
{
	char str[80];
	long a,b,c,d,e,y,m;
	long Y,M,D,H,Min,Sec;
	double fracday,fracrnd;
	double jday;
	struct tm date;

	putenv("TZ=GMT0");

	jday = *jdayp;

	fracday = jday - floor(jday);
	fracrnd = jday+.5 - floor(jday+.5);

 	if ( jday < 500.) { /* just days hours min */ 
		Y = 0;
		M = 0;
		D = jday;
		H = fracday*24;
		Min = fracday*24*60 -H*60 + .5/60;
		Sec = fracday*24*60*60 -(H*60 +Min)*60;
	}
	else { /* full year and month jday */
		a = jday + 0.5;
		b = a + 1537;
		c = (b-122.1)/365.25;
		d = c* 365.25;
		e = (b-d)/30.6001;
		/*printf("a=%d  b=%d  c=%d  d=%d  e=%d \n",a,b,c,d,e);*/
		D = b-d-floor(30.6001*e)+fracday;
		M = e-1-12*floor(e/14);
		Y = c - 4715 -floor((7+M)/10);
		H = fracrnd*24;
		Min = fracrnd*24*60 -H*60 + .5/60;   /*  correct for round off errors in sec*/
		Sec = fracrnd*24*60*60 -(H*60 +Min)*60;
		if (Sec==60) { 
			Min=Min+1; 
			Sec=0; 
		}
		if (Min>=60) { 
			Min=Min-60; 
			H = H+1;
		}
	} 

	date.tm_year = Y - 1900;
	date.tm_mon  = M - 1;
	date.tm_mday = D;
	date.tm_hour = H;
	date.tm_min  = Min;
	date.tm_sec  = Sec;
	date.tm_isdst = 0;


	strftime(str,100,"%Y %m %d %H %M",&date);

	*datep = date;
}

