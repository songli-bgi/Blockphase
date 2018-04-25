#ifndef __HAP_INIT_H__

#define __HAP_INIT_H__

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include "GenoHap_Init.h"
#include "Haplotype.h"
#include <ostream>
#include <boost/random.hpp>
#include "Random_Wizzard.h"

using namespace std;
using namespace boost;

struct node_t
{
	int * piece;
	int start;
	int end;
	int length;
	int count;
	vector< int > indivs;
	node_t( int len, int _start, int _end )
	{
		length = len;
		start = _start;
		end = _end;
		piece = new int [ length ];
		for ( int i = 0; i < length; i++ ) {
			piece[ i ] = -1;
		}
		count = 0;
	}

	node_t( const node_t & _n )
	{
		length = _n.length;
		start = _n.start;
		end = _n.end;
		count = _n.count;
		piece = new int [ length ];
		memcpy( piece, _n.piece, 4 * length );
		for ( size_t i = 0; i < _n.indivs.size(); i++ ) {
			indivs.push_back( _n.indivs[ i ] );
		}
	}

	bool compare( const node_t & _n )
	{
		if ( _n.start < start || _n.end > end ) {
			return false;
		}
		for ( int i = 0; i < _n.length; i++ ) {
			if ( _n.piece[ i ] != piece[ _n.start - start + i ] ) {
				return false;
			}
		}
		return true;
	}

	bool operator == ( const node_t & _n )
	{
		if ( _n.start >= end || _n.end <= start ) {
			return false;
		}
		int pesudo_start = _n.start > start ? _n.start : start;
		int pesudo_end   = _n.end   < end   ? _n.end   : end;
		for ( int i = pesudo_start; i <= pesudo_end; i++ ) {
			if ( _n.piece[ i - _n.start ] != piece[ i - start ] ) {
				return false;
			}
		}
		return true;
	}

	bool conjunct( const node_t & _n );
	bool is_conjunct( const node_t & _n )
	{
		if ( start <= _n.start && end >= _n.end ) {
			return false;
		} else {
			return true;
		}
	}

	friend ostream & operator << ( ostream & os, const node_t & _n )
	{
		os << "node show\t|";
		for ( int i = 0; i < _n.start; i++ ) {
			os << " ";
		}
		for ( int i = 0; i < _n.length; i++ ) {
			os << _n.piece[ i ];
		}
	#ifdef DEBUG
		os << "\tcount: " << _n.count;
	#endif
		os << "\n";
		return os;
	}
	
	~node_t()
	{
		delete [] piece;
		indivs.clear();
	}
};

struct _Init_Arg
{
	int begin;
	int end;
	int snp_number;
	int individuals;
	int marker_start;
	int gap;
	int16_t *** matrix;
	bool start_from_zero;

	vector< vector< int > > * p_genotypes;
	hap_t * haplotypes;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform;
	
	_Init_Arg()
	{
		begin = 0;
		end = 0;
		snp_number = 0;
		individuals = 0;
		marker_start = 0;
		gap = 0;
		matrix = NULL;
		start_from_zero = true;
		
		p_genotypes = NULL;
		haplotypes = NULL;
		p_uniform = NULL;
	}
	
	void init_args( int _begin, int _end, int _snp_number, int _individuals, int _marker_start, int _gap, bool _start_from_zero )
	{
		begin = _begin;
		end = _end;
		snp_number = _snp_number;
		individuals = _individuals;
		marker_start = _marker_start;
		gap = _gap;
		start_from_zero = _start_from_zero;
	}
	
	void init_pointers( int16_t *** _matrix, vector< vector< int > > * _p_genotypes, hap_t * _haplotypes, \
			    variate_generator< kreutzer1986, uniform_real<> > * _p_uniform )
	{
		matrix = _matrix;
		p_genotypes = _p_genotypes;
		haplotypes = _haplotypes;
		p_uniform = _p_uniform;
	}
	~_Init_Arg()
	{
		;
	}
};

struct _Hap_Init
{
	int individuals;
	int snp_number;
	int threshold;
	_Init_Arg * args;
	
	_Hap_Init( int _threshold )
	{
		individuals = 0;
		snp_number  = 0;
		threshold = _threshold;
		args = NULL;
	}

	void haplotype_init( vector< vector< int > > & genotypes, hap_t * haplotypes );
	void haplotype_init_piece( vector< vector< int > > & genotypes, hap_t * haplotypes, int marker_start, int marker_end, int gap, bool start_from_zero );
	void split_works( int _marker_start, int _gap, int _start_from_zero, int16_t *** _matrix, vector< vector< int > > * _p_genotypes, \
			  hap_t * _haplotypes, _Real_Uniforms & RU );
	void free_args();
	void haplotype_init( vector< vector< int > > & genotypes, hap_t * haplotypes, _Real_Uniforms & RU );
	void haplotype_init_piece( vector< vector< int > > & genotypes, hap_t * haplotypes, int marker_start, int marker_end, int gap, \
				   bool start_from_zero, _Real_Uniforms & RU );
	void hap_int2uint8_t( Haplotype * uint8_haps, hap_t * int_haps );


	/// Generate haplotypes randomly
	void haplotype_init_rand(vector<vector<int> > & genotypes, hap_t * haplotypes, variate_generator< kreutzer1986 &, uniform_real<> > & uniform);

	~_Hap_Init();
};

void * hap_init_function( void * );

#endif
