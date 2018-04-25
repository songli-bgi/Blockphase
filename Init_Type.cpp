#include "Init_Type.h"
#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

void Init_Type::init( double _mutation_rate, double _recombination_rate )
{
	for ( int i = 0; i < snp_number; i++ ) {
		mutation_rate[ i ] = _mutation_rate;
	}
	for ( int i = 0; i < snp_number - 1; i++ ) {
		theta_array[ i ] = _recombination_rate;
	}
}

void Init_Type::update_mutation()
{
	Error baseModel;
	for ( int i = 0; i < snp_number; i++ ) {
		if ( error_model[ i ].mismatch <= 2 ){
			baseModel += error_model[ i ];
		}else{
			mutation_rate[ i ] = error_model[ i ].update();
		}
	}
	baseModel.update();

	for ( int i = 0; i < snp_number; i++ ) {
		if ( error_model[ i ].mismatch <= 2 ){
			mutation_rate[ i ] = baseModel.rate;
		}
	}
}

void Init_Type::update_theta()
{
	Theta baseModel;
	for ( int i = 0; i < snp_number - 1; i++ ) {
		if ( crossover[ i ].cross <= 1 ) {
			baseModel += crossover[ i ];
		} else {
			theta_array[ i ] = crossover[ i ].update();
		}
	}
	float rate = baseModel.update();

	for ( int i = 0; i < snp_number - 1; i++ ) {
		if ( crossover[ i ].cross <= 1 ) {
			theta_array[ i ] = rate;
		}
	}
	
#ifdef DEBUG
	for ( int i = 0; i < snp_number - 1; i++ ) {
		cerr << "marker " << i << " recombination rate is " << theta_array[ i ] << endl;
	}
#endif
}

void Init_Type::clear()
{
	for ( int i = 0; i < snp_number; i++ ) {
		error_model[ i ].reset();
	}
	for ( int i = 0; i < snp_number - 1; i++ ) {
		crossover[ i ].reset();
	}
}

void Init_Type::show_mutation()
{
	cerr << "mutation      rate :";
	for ( int i = 0; i < snp_number; i++ ) {
		cerr.setf( ios_base::showpoint );
		cerr.setf( ios_base::fixed, ios_base::floatfield );
		cerr.precision( 6 );
		cerr << "\t\t" << "<" << i << ">" << mutation_rate[ i ];
		cerr.setf( ios_base::right, ios_base::adjustfield );
		cerr.width( 5 );
		cerr << "(" << error_model[ i ].mismatch << ")";
	}
	cerr << "\n";
}

void Init_Type::show_theta()
{
	cerr << "recombination rate:";
	for ( int i = 0; i < snp_number - 1; i++ ) {
		cerr.setf( ios_base::showpoint );
		cerr.setf( ios_base::fixed, ios_base::floatfield );
		cerr.precision( 6 );
		cerr << "\t\t" << "<" << i << ">" << theta_array[ i ];
		cerr.setf( ios_base::right, ios_base::adjustfield );
		cerr.width( 5 );
		cerr << "(" << crossover[ i ].cross << ")";
	}
	cerr << "\n";
}

void Init_Type::load_reference( int start, int stop, Haplotype * whole_reference, int whole_n_reference )
{
	if ( NULL == whole_reference ) {
		reference = NULL;
		n_reference = 0;
		return;
	}

	int pseudo_snp_number = stop - start;
	if ( pseudo_snp_number != snp_number ) {
		cerr << "[Init_Type::load_reference] fail to load this piece of reference from whole reference set.\n";
		abort();
	}

	n_reference = whole_n_reference;
	reference = new Haplotype[ n_reference ];
	Haplotype hap( pseudo_snp_number );
	for ( int i = 0; i < n_reference; i++ ) {
		reference[ i ] = hap;
		memcpy( reference[ i ].haplotype, whole_reference[ i ].haplotype + start, pseudo_snp_number );
	}
}
