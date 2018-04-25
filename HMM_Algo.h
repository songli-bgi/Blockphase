#ifndef __HMM_ALGORITHM_H__

#define __HMM_ALGORITHM_H__

#include <stdio.h>
#include "Haplotype_graph.h"
#include <vector>
#include <string>
#include <boost/random.hpp>
#include "Hap_Init.h"
#include <map>
#include "Error.h"

using namespace std;
using namespace boost;

struct _pair_compare_
{
	bool operator()( const pair< int, int > & a, const pair< int, int > & b ) const
	{
		return a.second < b.second;
	}
};

struct _reverse_pair_compare_
{
	bool operator()( const pair< int, int > & a, const pair< int, int > & b ) const
	{
		return a.second > b.second;
	}
};

struct _XK
{
	int start;
	int end;
	int id;
	int length;
	
	_XK()
	{
		start = 0;
		end = 0;
		id = -1;
		length = 0;
	}

	_XK( int _start, int _end, int _id )
	{
		start = _start;
		end = _end;
		id = _id;
		length = end - start;
	}
	
	_XK( const _XK & x )
	{
		this->start = x.start;
		this->end = x.end;
		this->id = x.id;
		this->length = x.length;
	}
};

struct _xk_compare_
{
	bool operator()( const _XK & a, const _XK & b ) const
	{
		if ( a.start < b.start ) {
			return true;
		} else if ( a.start > b.start ) {
			return false;
		} else {
			if ( a.length > b.length ) {
				return true;
			} else {
				return false;
			}
		}
	}
};

struct _xk_sort_
{
	bool operator()( const _XK & a, const _XK & b ) const
	{
		return a.length > b.length;
	}
};

struct _xk_reverse_compare_
{
	bool operator()( const _XK & a, const _XK & b ) const
	{
		if ( a.end > b.end ) {
			return true;
		} else if ( a.end < b.end ) {
			return false;
		} else {
			if ( a.length > b.length ) {
				return true;
			} else {
				return false;
			}
		}
	}
};

struct _HMM_Algo
{
	int      _individuals;
	int      _snp_number;
	int	 _start;
	int 	 _marginal_len;
	int      _indiv_ost;
	int    * _start_array;
	int    * _end_array;
	int    * _block_IDs;
	int    * _block_offsets;
	bool   * _block_flags;
	float    _haplotype_number;
	float ** _LeftMatrix;
	Haplo_graph _hg;
	float *** _errors; // mutation rate, no mutation for 0, one mutation for 1, two mutation for 2
	float * _marginal_i;
	float * _marginal_j;

	uint8_t ** prev_piece;
	uint8_t ** next_piece;
	uint32_t * prev_weigth;
	uint32_t * next_weigth;
	// vector<vector<pair<int, double> > > freqs;

	_HMM_Algo( int individuals, int snp_number, int start, int indiv_ost )
	{
		_snp_number  = snp_number;
		_individuals = individuals;
		_start = start;
		_marginal_len = 2 * individuals;
		_hg.snp_number  = _snp_number;
		_hg.individuals = _individuals;
		_indiv_ost = indiv_ost;
		_haplotype_number = (float)( (_individuals - 1) * 2 );
		_block_IDs = new int [ _snp_number ];
		for ( int i = 0; i < _snp_number; i++ ) {
			_block_IDs[ i ] = 0;
		}
		_block_offsets = new int [ _snp_number ];
		for ( int i = 0; i < _snp_number; i++ ) {
			_block_offsets[ i ] = 0;
		}
		_block_flags = new bool [ _snp_number ];
		for ( int i = 0; i < _snp_number; i++ ) {
			_block_flags[ i ] = false;
		}
		_errors = new float ** [ _snp_number ];
		for ( int i = 0; i < _snp_number; i++ ) {
			_errors[ i ] = new float * [ 3 ];
			for ( int j = 0; j < 3; j++ ) {
				_errors[ i ][ j ] = new float [ 3 ];
			}
		}

		prev_piece = NULL;
		next_piece = NULL;
		prev_weigth = NULL;
		next_weigth = NULL;
	}
	
	void transmit_state();
	void set_hapNumber( int n_hidden_states );
	void calculate_errors( float * mutation_rate );
	void alloc_leftmatrix();
	void show_leftmatrix();
	void free_leftmatrix();
	void check_haplotypes( Haplotype * haplotypes );

	void clear_matrix()
	{
		memset( _marginal_i, 0, sizeof( float ) * _marginal_len );
		memset( _marginal_j, 0, sizeof( float ) * _marginal_len );
	}

	void update_haplotype( Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
		 Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, float * mutation_rate );
	void update_marker( int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
		 variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, float * mutation_rate );

	/// Update haplotypes by maximum probability
	void update_haplotype_max(Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
		 Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, \
		 float * mutation_rate );
	void update_marker_max(int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
		 variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, float * mutation_rate );

	void prepare_pointors( int prev_marker, int next_marker );
	void destory_pointors( int prev_marker, int next_marker );

	void prior(); // initial states' probabilities
	void match_data( int marker, int individuals_ID, _GenoHap & genohap ); // emission probabilies
	void transpose( int prev_marker, int next_marker, float * theta_array, _GenoHap & genohap, int individuals_ID );
	void scale( int marker );
	void process( float * theta_array, _GenoHap & genohap, int individuals_ID );

	/**
	 * Design codes for huge sample size, all functions will be named with prefix "prune";
	*/
	vector< vector< int > > prune_states;

	void shared_haplotypes( ele & hap, ele & state, vector< pair< int, int > > & blks, int & ibd_distance );
	void shared_haplotypes( ele & hap, ele & state, vector< _XK > & _XKs, int id, bool output );
	void prune_states_find( int individuals_ID, Haplotype * haplotypes, int maximum );
	void prune_states_find( int individuals_ID, int maximum, Haplotype * haplotypes );

	void prune_states_find( Haplotype * haplotypes, int individuals_ID, int maximum );
	void prune_states_find( Haplotype * haplotypes, int individuals_ID, map< int, int > & freqs, int minimun, int maximum );
	void prune_alloc_leftmatrix();
	void prune_free_leftmatrix();
	
	void prune_prior();
	void prune_match_data( int marker, int individuals_ID, _GenoHap & genohap );
	void prune_prepare_pointors( int prev_marker, int next_marker );
	void prune_destory_pointors();
	void prune_transpose( int prev_marker, int next_marker, float * theta_array, _GenoHap & genohap, int individuals_ID );
	void prune_scale( int marker );
	void prune_process( float * theta_array, _GenoHap & genohap, int individuals_ID );
	void prune_update_haplotype( Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
			Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, float * mutation_rate );
	void prune_update_marker( int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
				variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, float * mutation_rate );
	
	~_HMM_Algo();
};

#endif
