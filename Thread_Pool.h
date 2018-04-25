#ifndef __THREAD_POOL_H__

#define __THREAD_POOL_H__

#include "Error.h"
#include "GenoHap_Init.h"
#include "Init_Type.h"
#include <pthread.h>
#include "Random_Wizzard.h"
#include "Thread_Synchronization.h"
#include "Sample_Orders.h"

struct _Thread_Arg
{
	int individuals;
	int snp_number;
	int begin;
	int end;
	int states;
	int maximum;
	int start;
	int thread_id;
	int blk_size;
	int indiv_ost;
	int sample_time;
	int subChunk;
	bool burnin;
	bool reference_only;
	bool gl;
	bool odd;
	Error * error_model;
	Theta * crossover;
	Init_Type * engine;
	_GenoHap  * genohap;
	_Thread_Syn * p_TS;
	_Sample_Orders * p_so;

	variate_generator< kreutzer1986, uniform_real<> > * p_uniform;
	vector< pair< int, int > > vec_iters;

	_Thread_Arg()
	{
		individuals = 0;
		snp_number = 0;
		states = 0;
		start = 0;
		
		begin = 0;
		end = 0;
		thread_id = 0;
		sample_time = 1;
		subChunk = 4;
		
		burnin = false;
		reference_only = false;
		gl = false;
	
		error_model = NULL;
		crossover = NULL;

		engine = NULL;
		genohap = NULL;
		p_uniform = NULL;
		p_TS = NULL;
		p_so = NULL;
	}

	void init( int _individuals, int _snp_number, int _states, int _begin, int _end, int _start, int _thread_id, bool _burnin, \
		   bool _reference_only, int _maximum, int _indiv_ost, int _sample_time, bool _gl, int _subChunk, int _blk_size, bool _odd )
	{
		individuals = _individuals;
		snp_number = _snp_number;
		states = _states;
		start = _start;
		
		begin = _begin;
		end = _end;
		thread_id = _thread_id;
		sample_time = _sample_time;
		subChunk = _subChunk;
		burnin = _burnin;
		reference_only = _reference_only;
		maximum = _maximum;
		blk_size = _blk_size;
		indiv_ost = _indiv_ost;
		gl = _gl;
		odd = _odd;

		error_model = new Error [ _snp_number ];
		crossover = new Theta [ _snp_number - 1 ];
	}

	void init_iters();

	void destory()
	{
		if ( error_model != NULL ) {
			delete [] error_model;
			error_model = NULL;
		}
		if ( crossover != NULL ) {
			delete [] crossover;
			crossover = NULL;
		}
	}
};

typedef void * ( * start_routine_t ) ( void * );

struct _Thread_Pool
{
	int individuals;
	int snp_number;
	int states;
	int n_threads;
	int start;
	int maximum;
	int blk_size;
	int indiv_ost;
	int sample_time;
	int subChunk;

	bool burnin;
	bool reference_only;
	bool gl;
	bool odd;

	_Thread_Arg * args;

	pthread_t * thread_ids;
	
	_Thread_Pool( int _individuals, int _snp_number, int _states, int _n_threads, int _start, bool _burnin, int _maximum, \
		      bool _reference_only, int _indiv_ost, int _sample_time, bool _gl, int _subChunk, int _blk_size, bool _odd )
	{
		individuals = _individuals;
		snp_number = _snp_number;
		states = _states;
		n_threads = _n_threads;
		start = _start;
		burnin = _burnin;
		maximum = _maximum;
		sample_time = _sample_time;
		subChunk = _subChunk;
		blk_size = _blk_size;
		reference_only = _reference_only;
		gl = _gl;
		odd = _odd;
		indiv_ost = _indiv_ost;

		thread_ids = ( pthread_t * )malloc( sizeof( pthread_t ) * n_threads );
	}

	void init_args( Init_Type * _engine, _GenoHap * _genohap, _Real_Uniforms & RU, _Thread_Syn * _p_TS, _Sample_Orders * _p_so );
	void destory_args();

	void threads_start( start_routine_t thread_start_function );
	void threads_end();

	void integrate_paramaters( Init_Type * _engine );

	~_Thread_Pool()
	{
		free( thread_ids );
		destory_args();
	}
};

#endif
