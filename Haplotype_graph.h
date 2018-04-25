#ifndef __HAPLOTYPE_GRAPH_H__

#define __HAPLOTYPE_GRAPH_H__

#include <list>
#include <map>
#include <vector>
#include <assert.h>
#include "Haplotype.h"
#include <boost/random.hpp>

using namespace std;
using namespace boost;

struct ele
{
	uint8_t * p;
	int length;
	ele()
	{
		p = NULL;
		length = 0;
	}

	ele ( int n )
	{
		p = new uint8_t [ n ];
		memset( p, 0, n );
		length = n;
	}

	ele ( const ele & e )
	{
		p = new uint8_t [ e.length ];
		length = e.length;
		memcpy( p, e.p, e.length );
	}
	
	void Alloc( int n )
	{
		p = new uint8_t [ n ];
		length = n;
		memset( p, 0, length );
	}
	
	~ele()
	{
		if ( p != NULL ) {
			delete [] p;
			p = NULL;
		}
	}
};

struct Segment
{
	ele * piece;
	uint32_t * weigth;
	int size;
	int length;
	Segment( int n )
	{
		piece = new ele[ n ];
		weigth = new uint32_t[ n ];
		memset( weigth , 0 , 4 * n );
		size = n;
		length = 0;
	}

	~Segment()
	{
		delete [] piece;
		delete [] weigth;
	}
};

struct Comp
{
	bool operator() ( const ele & a, const ele & b ) const
	{
		assert( a.length == b.length );
                for ( int i = 0; i < a.length; i++ ){
                        if ( a.p[i] < b.p[i] ){
                                return true;
                        } else if ( a.p[i] > b.p[i] ){
                                return false;
                        }
                }
                return false;
        }
};

typedef map< ele, pair< int, int >, Comp > _KJD_TYPE;
typedef map< ele, pair< int, int >, Comp >::iterator _KJD_TYPE_ITERATOR;

struct _KJD
{
private:
	int id;
public:
	_KJD_TYPE haps;
	map< int, map< int, int > >  prev_edges;
	int length;
	_KJD()
	{
		id = 0;
		length = 0;
	}

	void add( const ele & _hap )
	{
		haps[ _hap ] = pair< int, int >( id++, 1 );
	}

	int getid( const ele & _hap )
	{
		return haps[ _hap ].first;
	}
	
	void setid( int _id )
	{
                id = _id;
        }

        int curr_id()
        {
                return id;
        }
	
	void clear()
	{
		haps.clear();
		map< int, map< int, int > >::iterator iter1 = prev_edges.begin();
		map< int, map< int, int > >::iterator iter2 = prev_edges.end();
		
		for ( ; iter1 != iter2; iter1++ ) {
			iter1->second.clear();
		}
		
		prev_edges.clear();
		id = 0;
	}
};

struct Haplo_graph
{
	int snp_number;
	int individuals;
	int n_reference;

	_KJD * _kjds;
	int n_kjds;
	vector< Segment * > _kjd_array;
	list< int > _nodes;

	

	Haplo_graph()
	{
		snp_number  = 0;
		individuals = 0;
		n_reference = 0;
		_kjds = NULL;
		n_kjds = 0;
	}
	
	Haplo_graph( int _individuals, int _snp_number )
	{
		snp_number = _snp_number;
		individuals = _individuals;
		n_reference = 0;
		_kjds = NULL;
		n_kjds = 0;
	}

	void split_blocks( Haplotype * haplotypes , int state_k, Haplotype * reference, int _n_reference, \
			   variate_generator< kreutzer1986, uniform_real<> > & uniform, bool reference_only );
	void split_blocks( vector< int > & blocks, int size, variate_generator< kreutzer1986, uniform_real<> > & uniform );
	void list2vector( vector< int > & blocks );
	void construct_hg( Haplotype * haplotypes, int individuals_ID, vector< int > & blocks, Haplotype * reference, int _n_reference, bool reference_only );
	void construct_hg( Haplotype * haplotypes, map< int, int > & exclude_samples, vector< int > & blocks, Haplotype * reference, int _n_reference, \
			   bool reference_only );
	void modify_hg( Haplotype * haplotypes, int prev_individuals_ID, int next_individuals_ID, vector< int > & blocks );
	void modify_hg( Haplotype * haplotypes, vector< int > & prev_samples, vector< int > & next_samples, vector< int > & blocks );
	void modify_KJDS();
	void free_kjds();
	void show_kjds();
	void load_kjd_array();
	void free_kjd_array();
	void show_kjd_array();

	void clear();
	~Haplo_graph()
	{
		this->free_kjd_array();
		this->free_kjds();
	}
};

#endif
