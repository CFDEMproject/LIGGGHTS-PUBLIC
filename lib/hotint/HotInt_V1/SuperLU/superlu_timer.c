//#include "stdafx.h"
/* 
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */


#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#else

//#include <sys/types.h>
//#include <sys/times.h>

/*#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <memory.h>*/

#include <time.h>
#include <sys/timeb.h>
#include <math.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double SuperLU_timer_()
{
	struct timeb tb;
	tb.time = 0;
	tb.millitm = 0;
	ftime(&tb);
	return tb.time+(double)tb.millitm*0.001;
/*
	struct tms use;
    double tmp;
    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / CLK_TCK;
		*/
}

#endif

