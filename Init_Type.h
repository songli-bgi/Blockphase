#ifndef __INIT_TYPE_H__

#define __INIT_TYPE_H__

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Haplotype_graph.h"
#include <string>
#include "Error.h"

using std::string;

struct Init_Type
{
	int snp_number;
	int individuals;
	int n_reference;
	Haplotype * haplotypes;
	Haplotype * reference;
	float * mutation_rate; // size is snp_number;
	float * theta_array;   // size is snp_number - 1;
	Error * error_model;
	Theta * crossover;

	Init_Type( int _individuals , int _snp_number )
	{
		snp_number = _snp_number;
		individuals = _individuals;
				
		haplotypes = new Haplotype[ individuals * 2 ];
		Haplotype hap( snp_number );
		for ( int i = 0; i < individuals; i++ ){
			haplotypes[ i * 2 + 0 ] = hap;
			haplotypes[ i * 2 + 1 ] = hap;
		}
		mutation_rate = new float [ snp_number ];
		theta_array = new float [ snp_number - 1 ];

		error_model = new Error[ snp_number ];
		crossover = new Theta[ snp_number - 1 ];

		n_reference = 0;
		reference = NULL;
	}
	void init( double _mutation_rate, double _recombination_rate );
	void update_mutation();
	void update_theta();
	void clear();
	void show_mutation();
	void show_theta();
	void load_reference( int start, int stop, Haplotype * whole_reference, int whole_n_reference );

	~Init_Type()
	{
		delete [] haplotypes;
		delete [] mutation_rate;
		delete [] theta_array;
		delete [] error_model;
		delete [] crossover;
		if ( reference != NULL ) {
			delete [] reference;
			reference = NULL;
		}
	}
};

#endif
