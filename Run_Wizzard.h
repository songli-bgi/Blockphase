#ifndef __RUN_WIZZARD_H__

#define __RUN_WIZZARD_H__

#include <vector>
#include <string>
#include <assert.h>
#include "Arg_Wizzard.h"

using namespace std;

struct _Run_Wizzard
{
	int individuals; // total number of individuals.
	int snp_number;  // total number of snps.

	char ** ref_Alt_bases;
	vector< pair<int, int> > Run_Pieces;
	vector< string > samples;
	vector< int > snpids;

	_Run_Wizzard()
	{
		individuals = 0;
		snp_number = 0;
		ref_Alt_bases = NULL;
	}

	void init( int _individuals, int _snp_number )
	{
		individuals = _individuals;
		snp_number = _snp_number;
	}

	void split( _Program_Variables & pv );
	int  get_index( int pos );
	void show_variables( _Program_Variables & pv, FILE * fp );

	~_Run_Wizzard()
	{	
		if ( ref_Alt_bases != NULL ) {
			for ( int i = 0; i < snp_number; i++ ) {
				delete [] ref_Alt_bases[ i ];
			}
			delete [] ref_Alt_bases;
		}
		ref_Alt_bases = NULL;
	}
};

#endif
