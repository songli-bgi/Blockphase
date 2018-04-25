#ifndef __TIME_COST_H__

#define __TIME_COST_H__

#include <time.h>
#include <stdio.h>

class _Time_Cost
{
	time_t time_start;
	time_t time_end;

public:

	_Time_Cost()
	{
		time( &time_start );
		time( &time_end );
	}

	void start()
	{
		time( &time_start );
	}
	
	void end()
	{
		time( &time_end );
	}

	void show( FILE * fp )
	{
		double cost = difftime( time_end, time_start );
		int c = int(cost);
		fprintf( fp, "  * Second : %d\n", c % 60 );
		c /= 60;
		fprintf( fp, "  * Minute : %d\n", c % 60 );
		c /= 60;
		fprintf( fp, "  * Hour   : %d\n", c % 24 );
		c /= 24;
		fprintf( fp, "  * Day    : %d\n", c );
	}
};

#endif
