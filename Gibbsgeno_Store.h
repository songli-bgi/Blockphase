#ifndef __GIBBSGENO_STORE_H__

#define __GIBBSGENO_STORE_H__

#include <stdint.h>
#include "Init_Type.h"

struct _Gibbsgeno_Store
{
	int snp_number;
	int individuals;
	int circles;
	int length;

	uint8_t *** gibbs_genotypes;
	
	_Gibbsgeno_Store( int _individuals, int _snp_number, int _circles )
	{
		snp_number = _snp_number;
		individuals = _individuals;
		circles = _circles;
		
		length = snp_number / 8 + 1;
		
		// alloc memory for gibbs genotypes.
		gibbs_genotypes = new uint8_t ** [ circles ];
		for ( int i = 0; i < circles; i++ ) {
			gibbs_genotypes[ i ] = new uint8_t * [ individuals * 2 ];
			for ( int j = 0; j < individuals * 2; j++ ) {
				gibbs_genotypes[ i ][ j ] = new uint8_t [ length ];
			}
		}
	}

	void set( int marker, int flag, int cycle, int phase, int individuals_ID );
	int  get( int marker, int cycle, int phase, int individuals_ID );
	
	void load( Init_Type & engine, int cycle );
	void load( Haplotype * haplotypes, int individuals_ID, int cycle );
	void binary2uint8_t( uint8_t *** _gibbs_genotypes );

	~_Gibbsgeno_Store()
	{
		for ( int i = 0; i < circles; i++ ) {
			for ( int j = 0; j < individuals * 2; j++ ) {
				delete [] gibbs_genotypes[ i ][ j ];
			}
			delete [] gibbs_genotypes[ i ];
		}
		delete [] gibbs_genotypes;
	}
};

#endif
