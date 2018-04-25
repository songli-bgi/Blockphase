#ifndef __GENOHAP_INIT_H__
#define __GENOHAP_INIT_H__

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include "Haplotype.h"
#include "freader.h"
#include "Arg_Wizzard.h"
#include "Run_Wizzard.h"

using namespace std;

struct hap_t
{
	int * first;
	int * second;
	int length;

	hap_t()
	{
		first = NULL;
		second = NULL;
		length = 0;
	}

	hap_t( int len )
	{
		length = len;
		first = new int [ length ];
		second = new int [ length ];
		for ( int i = 0; i < length; i++ ) {
			first[ i ] = -1;
			second[ i ] = -1;
		}
	}

	hap_t( const hap_t & h )
	{
		length = h.length;
		first = new int [ length ];
		second = new int [ length ];
		memcpy( first, h.first, 4 * h.length );
		memcpy( second, h.second, 4 * h.length );
	}
	void operator = ( const hap_t & h )
	{
		this->length = h.length;
		this->first = new int [ length ];
		this->second = new int [ length ];
		memcpy( this->first, h.first, 4 * h.length );
		memcpy( this->second, h.second, 4 * h.length );
	}

	~hap_t()
	{
		delete [] first;
		delete [] second;
	}

};

struct _GenoHap
{
	int individuals;
	int snp_number;
	hap_t * haplotypes;
	vector< vector< int > > genotypes;

	float * glf; 			 ///< glf of all samples, all snps and all genotypes

	Haplotype * whole_reference;
	int whole_n_reference;

	_GenoHap()
	{
		individuals = 0;
		snp_number  = 0;
		haplotypes = NULL;
		glf = NULL;

		whole_reference = NULL;
		whole_n_reference = 0;
	}

	~_GenoHap()
	{
		this->destroyHap();
	}

	/// Initialize by GLF info.
//	void selfinit_glf( string & infile, vector< int > & snpids, vector< string > & samples, char ** & ref_Alt_bases );	
	void selfinit_glf( string & infile, _Run_Wizzard & RW );
	/// Create genotypes by GLF info.
	void random_createGeno( int start, int stop );
	void em_createGeno( int start, int stop );
	void load_reference( const char * ref_file, vector< int > & snpids, char ** ref_Alt_bases );
	void show_variables( FILE * fp, _Program_Variables & pv );

	void createHap( int _snp_number );
	void randomHap( int _snp_number );
	void free_hap();
	void destroyHap();
};

#endif
