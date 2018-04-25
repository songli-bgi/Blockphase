#include "Thread_Pool.h"
#include <vector>

using namespace std;

void _Thread_Arg::init_iters()
{
        int nsamples_perChunk = (end - begin) / subChunk;
        int n_remains  = (end - begin) % subChunk;
        int locate = begin;
        for ( int i = 0; i < subChunk; i++ ) {
                if ( n_remains > 0 ) {
                        vec_iters.push_back( pair< int, int >( locate, locate + nsamples_perChunk + 1 ) );
                        n_remains--;
                        locate += nsamples_perChunk + 1;
                } else {
                        vec_iters.push_back( pair< int, int >( locate, locate + nsamples_perChunk ) );
                        locate += nsamples_perChunk;
                }
        }
        assert( locate == end );
}

void _Thread_Pool::init_args( Init_Type * _engine, _GenoHap * _genohap, _Real_Uniforms & RU, _Thread_Syn * _p_TS, _Sample_Orders * _p_so )
{
	int n_piece = individuals / n_threads;
	int n_remain = individuals % n_threads;
	vector< int > segs;
	for ( int i = 0; i < n_threads; i++ ) {
		if ( n_remain > 0 ) {
			segs.push_back( n_piece + 1 );
			n_remain--;
		} else {
			segs.push_back( n_piece );
		}
	}
	
	int check = 0;
	for ( int i = 0; i < segs.size(); i++ ) {
		check += segs[ i ];
	}

	if ( check != individuals ) {
		cerr << "[_Thread_Pool::init_args] fail to split individuals to " << n_threads << " .\n";
		abort();
	}

	args = new _Thread_Arg[ n_threads ];
	int _start( 0 ), _end( 0 ), locate( 0 );
	for ( int i = 0; i < n_threads; i++ ) {
		_start = locate;
		_end   = locate + segs[ i ];
		args[ i ].init( individuals, snp_number, states, _start, _end, start, i, burnin, reference_only, maximum, indiv_ost, sample_time, gl, subChunk, blk_size, odd );
		args[ i ].engine = _engine;
		args[ i ].genohap = _genohap;
		args[ i ].p_uniform = &(RU.uniforms[ i ]);
		args[ i ].p_TS = _p_TS;
		args[ i ].p_so = _p_so;
		args[ i ].init_iters();
		locate += segs[ i ];
	}
}

void _Thread_Pool::destory_args()
{
	for ( int i = 0; i < n_threads; i++ ) {
		args[ i ].destory();
	}
	delete [] args;
}

void _Thread_Pool::threads_start( start_routine_t thread_start_function )
{
	for ( int i = 0; i < n_threads; i++ ) {
		if ( 0 != pthread_create( &(thread_ids[ i ]), NULL, thread_start_function, ( void * )(&(args[ i ])) ) ) {
			cerr << "[_Thread_Pool::threads_start] fail to create " << i << " thread.\n";
			abort();
		}
	}
}

void _Thread_Pool::threads_end()
{
	for ( int i = 0; i < n_threads; i++ ) {
		if ( 0 != pthread_join( thread_ids[ i ], NULL ) ) {
			cerr << "[_Thread_Pool::threads_end] fail to wait thread " << i << " exit.\n";
			abort();
		}
	}
}

void _Thread_Pool::integrate_paramaters( Init_Type * _engine )
{
	for ( int i = 0; i < snp_number; i++ ) {
		Error t;
		for ( int j = 0; j < n_threads; j++ ) {
			t += args[ j ].error_model[ i ];
		}
		_engine->error_model[ i ] += t;
	}
	
	for ( int i = 0; i < snp_number - 1; i++ ) {
		Theta t;
		for ( int j = 0; j < n_threads; j++ ) {
			t += args[ j ].crossover[ i ];
		}
		_engine->crossover[ i ] += t;
	}
}
