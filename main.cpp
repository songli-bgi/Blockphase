#include "GenoHap_Init.h"
#include "Hap_Init.h"
#include "Haplotype_graph.h"
#include "Init_Type.h"
#include <boost/random.hpp>
#include "HMM_Algo.h"
#include "Thread_Pool.h"
#include "Arg_Wizzard.h"
#include "Run_Wizzard.h"
#include "Random_Wizzard.h"
#include "Thread_Synchronization.h"
#include "Output_Wizzard.h"
#include "Time_Cost.h"
#include "Gibbsgeno_Store.h"
#include "options.h"
#include "EM_hfs.h"

using namespace std;
using namespace boost;

void check_init_haplotypes(int individuals, int snp_number, _GenoHap & genohap)
{
	for ( int i = 0; i < individuals; i++ ) {
		for ( int j = 0; j < snp_number; j++ ) {
			if ( genohap.haplotypes[ i ].first[ j ] == -1 ) {
				cerr << "[_Hap_Init] fail to init haplotypes.\n";
				abort();
			}
		}
		for ( int j = 0; j < snp_number; j++ ) {
			if ( genohap.haplotypes[ i ].second[ j ] == -1 ) {
				cerr << "[_Hap_Init] fail to init haplotypes.\n";
                                abort();
                        }
		}
	}
}

void * normal_start_function( void * arg );
void * prune_start_function( void * arg );
void * hybrid_start_function( void * arg );
void * prune_normal_start_function( void * arg );

int main( int argc, char * argv[] )
{
	time_t time1, time2;
	struct tm * tm2;
	time( &time1 );

	_Options OPH;
	_Program_Variables pv;
        anyarg opt;
        _Arg_Wizzard::init( opt );
	if ( argc == 1 ) {
                OPH.help();
                exit( 0 );
        }
        if ( ! opt.parse_command_line( argc, argv ) ) {
                abort();
        }
        _Arg_Wizzard::process( opt, pv, OPH );

	pv.show_variables( stderr );

	_Run_Wizzard RW;
	
	if ( pv.in.empty() ) {
		cerr << "[main] Your must provide an input file\n";
		abort();
	}
	string logfile;
	if ( pv.out.empty() ) {
		logfile = "phase_haplotypes.log";
	} else {
		logfile = pv.out + ".log";
	}
	FILE * logp = fopen( logfile.c_str(), "w" );
	pv.show_variables( logp );
	string infile = pv.in;
	_GenoHap genohap;
	genohap.selfinit_glf( infile, RW );
	int original_individuals = genohap.individuals;
	int original_snp_number  = genohap.snp_number;

	_Sample_Orders SO( original_individuals );

	if ( ! pv.reference.empty() ) {
		genohap.load_reference( pv.reference.c_str(), RW.snpids, RW.ref_Alt_bases );
	}

	RW.init( original_individuals, original_snp_number );
	RW.split( pv );
	RW.show_variables( pv, stderr );
	genohap.show_variables( stderr, pv );
	
	RW.show_variables( pv, logp );
	genohap.show_variables( logp, pv );
	fflush( logp );

	cerr << "\n\033[0mSplit Data to ";
	cerr << "\033[0;31m" << RW.Run_Pieces.size();
	cerr << "\033[0m Chunks\n";

	/* Random Genorator */
	_Real_Uniforms RU( pv.n_threads, pv.seed );
	kreutzer1986 rh_engine;
	uniform_real<> rh_real( 0, 1 );
	variate_generator< kreutzer1986&, uniform_real<> > rh_uniform( rh_engine, rh_real );

	/* Threads synchronization */
	_Thread_Syn TS( pv.n_threads );

	/* Haplotypes Conjuncting */
	_Hap_Conjunct HC( original_individuals, pv.overlap );
	_Time_Cost TC;

	for ( int piece_id = 0; piece_id < RW.Run_Pieces.size(); piece_id++ ) {
		TC.start();
		int start = RW.Run_Pieces[ piece_id ].first;
		int stop  = RW.Run_Pieces[ piece_id ].second;
		int pseudo_snp_number = stop - start;
		int individuals = original_individuals;
		int snp_number  = pseudo_snp_number;
		
		EM_hfs em_engine;
		em_engine.enumerate( pv.eml );

		if ( pv.random ) {
			genohap.random_createGeno( start, stop );
			genohap.createHap( pseudo_snp_number );
			genohap.randomHap( pseudo_snp_number );
		} else {
			genohap.em_createGeno( start, stop );
			genohap.createHap( pseudo_snp_number );
			em_engine.initialize( genohap, genohap.haplotypes, start, pseudo_snp_number );
		}

		cerr << "\n\033[0mChunk ";
		cerr << "\033[0;36m" << piece_id + 1 << "/" << RW.Run_Pieces.size();
		cerr << "\033[0m Contains markers in region [ " << RW.snpids[ start ] << ", " \
		     << RW.snpids[ stop>=original_snp_number?original_snp_number - 1:stop ] << " )" << endl;
		fprintf( logp, "\nChunk %d Contains markers in region [ %d, %d )\n", piece_id + 1, RW.snpids[ start ], \
				RW.snpids[ stop>=original_snp_number?original_snp_number - 1:stop ] );
		
		Init_Type engine( individuals, snp_number );
		_Hap_Init hapinit( pv.threshold );
	
		if ( pv.ibd ) {
			hapinit.haplotype_init( genohap.genotypes, genohap.haplotypes );
			hapinit.snp_number = snp_number;
			hapinit.individuals = individuals;
			hapinit.hap_int2uint8_t( engine.haplotypes, genohap.haplotypes );
		} else {
			hapinit.snp_number = snp_number;
			hapinit.individuals = individuals;
			hapinit.hap_int2uint8_t( engine.haplotypes, genohap.haplotypes );
		}
		check_init_haplotypes( individuals, snp_number, genohap );
		genohap.free_hap();
		engine.init( pv.mutation_rate, pv.recombination_rate );
	
		if ( ! pv.reference.empty() ) {
			engine.load_reference( start, stop, genohap.whole_reference, genohap.whole_n_reference );
		}	
		int states = pv.states;
		int BURNIN = pv.burnin;
		int CIRCLE = pv.circle;
		int maximum = pv.maximum;
		int sample_time = pv.sample_time;
		int subChunk = pv.subChunk;
		bool reference_only = pv.reference_only;
		bool gl = pv.gl;
		bool normal = pv.normal;
		bool prune  = pv.prune;
		bool hybrid = pv.hybrid;
		
		int n_threads = pv.n_threads;
		_Gibbsgeno_Store GI( individuals, snp_number, BURNIN );
		for ( int burnin = 0; burnin < BURNIN; burnin++ ) {
			cerr << "\033[0mBurnin Iterating ";
			cerr << "\033[0;36m" << burnin + 1 << "/" << BURNIN;
			cerr << "\033[0m \033[A" << endl;
			TS.reset();
			_Thread_Pool t_pool( individuals, snp_number, states, n_threads, start, true, maximum, reference_only, original_snp_number * 3, \
					     sample_time, gl, subChunk, 0, 0 );
			SO.randOrder();
			t_pool.init_args( &engine, &genohap, RU, &TS, &SO );
			if ( normal && !prune )
				t_pool.threads_start( normal_start_function );
			else if ( prune && !normal )
				t_pool.threads_start( prune_start_function );
			else if ( hybrid )
				t_pool.threads_start( hybrid_start_function );
			else
				t_pool.threads_start( prune_normal_start_function );
			t_pool.threads_end();
			t_pool.integrate_paramaters( &engine );
			engine.update_mutation();
			engine.update_theta();
			engine.clear();
			GI.load( engine, burnin );
		}
		_Output_Wizzard OWI;
		OWI.decide_all( GI, engine, BURNIN );
		cerr << "Burnin Iterating Finish                \n";

		bool odd;
		_Gibbsgeno_Store GS( individuals, snp_number, CIRCLE );
		for ( int circle = 0; circle < CIRCLE; circle++ ) {
			int blk_size = 5 + 3 * circle;
			cerr << "\033[0mMain   Iterating ";
			cerr << "\033[0;36m" << circle + 1 << "/" << CIRCLE;
			cerr << "\033[0m \033[A" << endl;
			odd = circle >=  (2 * CIRCLE / 3) ? true : false;
			TS.reset();
			_Thread_Pool t_pool( individuals, snp_number, states, n_threads, start, false, maximum, false, original_snp_number * 3, \
					     sample_time, gl, subChunk, blk_size, odd );
			SO.randOrder();
			t_pool.init_args( &engine, &genohap, RU, &TS, &SO );
			if ( normal && !prune )
                                t_pool.threads_start( normal_start_function );
                        else if ( prune && !normal )
                                t_pool.threads_start( prune_start_function );
                        else if ( hybrid )
                                t_pool.threads_start( hybrid_start_function );
                        else
				t_pool.threads_start( prune_normal_start_function );
			t_pool.threads_end();
			t_pool.integrate_paramaters( &engine );
			engine.update_mutation();
			engine.update_theta();
			engine.clear();
			GS.load( engine, circle );
		}
		cerr << "Main   Iterating Finish                \n";
		cerr << "\nPlease wait, Refine iterating will prepare output\n";
		cerr << "\n\033[0mRefine Iterating Finish, Outputting Chunk ";
		cerr << "\033[0;36m" << piece_id + 1;
		cerr << "\033[0m ...\n";
		
		// output
		string outfile;
		if ( pv.out.empty() ) {
			outfile = "blockphase";
		} else {
			outfile = pv.out;
		}

		bool add = false;
		int out_begin( 0 ), out_end( 0 );
		if ( RW.Run_Pieces.size() == 1 ) {
			out_begin = 0;
			out_end   = snp_number;
			add = false;
		} else {
			if ( piece_id == 0 ) {
				out_begin = 0;
				out_end   = snp_number - pv.overlap / 2;
				add = false;
			} else if ( piece_id == RW.Run_Pieces.size() - 1 ) {
				out_begin = pv.overlap / 2;
				out_end   = snp_number;
				add = true;
			} else {
				out_begin = pv.overlap / 2;
				out_end   = snp_number - pv.overlap / 2;
				add = true;
			}
		}
		
		_Output_Wizzard OW( outfile, original_individuals, &HC );
		OW.process( GS, engine, out_begin, out_end, add, RW, start, CIRCLE );
		
		TC.end();
		fprintf( logp, "\nChunk %d Using time :\n", piece_id + 1 );
		TC.show( logp );
		fflush( logp );
	}

	time( &time2 );
	tm2 = localtime( &time2 );
	int using_time = int(difftime( time2, time1 ));
	fprintf( stderr, "\nEnd time : %0.2d/%0.2d/%d", tm2->tm_mday, ( 1 + tm2->tm_mon), (1900 + tm2->tm_year) );
	fprintf( stderr, " %0.2d:%0.2d:%0.2d\n\n", tm2->tm_hour, tm2->tm_min, tm2->tm_sec );
	fprintf( logp, "\nEnd time : %0.2d/%0.2d/%d", tm2->tm_mday, ( 1 + tm2->tm_mon), (1900 + tm2->tm_year) );
        fprintf( logp, " %0.2d:%0.2d:%0.2d\n\n", tm2->tm_hour, tm2->tm_min, tm2->tm_sec );
	fprintf( stderr, "Using time :\n" );
	fprintf( stderr, "  * Second : %d\n", using_time % 60 );
	fprintf( logp, "Using time :\n" );
        fprintf( logp, "  * Second : %d\n", using_time % 60 );
	using_time /= 60;
	fprintf( stderr, "  * Minute : %d\n", using_time % 60 );
	fprintf( logp, "  * Minute : %d\n", using_time % 60 );
	using_time /= 60;
	fprintf( stderr, "  * Hour   : %d\n", using_time % 24 );
	fprintf( logp, "  * Hour   : %d\n", using_time % 24 );
	using_time /= 24;
	fprintf( stderr, "  * Day    : %d\n", using_time );
	fprintf( logp, "  * Day    : %d\n", using_time );
	
	fclose( logp );
	return 0;
}

void * normal_start_function( void * arg )
{
	_Thread_Arg * t_arg = ( _Thread_Arg * ) arg;
	int individuals = t_arg->individuals;
	int snp_number  = t_arg->snp_number;
	int states      = t_arg->states;
	int begin	= t_arg->begin;
	int end		= t_arg->end;
	int start	= t_arg->start;
	int thread_id   = t_arg->thread_id;
	int maximum     = t_arg->maximum;
	int blk_size    = t_arg->blk_size;
	int indiv_ost   = t_arg->indiv_ost;
	int sample_time = t_arg->sample_time;
	bool burnin     = t_arg->burnin;
	bool gl         = t_arg->gl;
	bool reference_only = t_arg->reference_only;
	Error * error_model = t_arg->error_model;
	Theta * crossover   = t_arg->crossover;
	Init_Type * p_engine = t_arg->engine;
	_GenoHap * p_genohap = t_arg->genohap;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform = t_arg->p_uniform;
	_Thread_Syn * p_TS = t_arg->p_TS;
	_Sample_Orders * p_so = t_arg->p_so;

	_HMM_Algo hmm_algorithm( individuals, snp_number, start, indiv_ost );

	vector< pair< int, int > > & iters = t_arg->vec_iters;
	int sample_beg = iters[ 0 ].first;
	int sample_end = iters[ 0 ].second;
	map< int, int > exclude_samples;
	vector< int > prev_samples, next_samples;
	for ( int s = sample_beg; s < sample_end; s++ ) {
		exclude_samples[ p_so->orders[ s ] ] = 0;
		prev_samples.push_back( p_so->orders[ s ] );
		next_samples.push_back( p_so->orders[ s ] );
	}

	vector< int > blocks;
	hmm_algorithm._hg.split_blocks( p_engine->haplotypes, states, p_engine->reference, p_engine->n_reference, *p_uniform, reference_only );
	hmm_algorithm._hg.list2vector( blocks );
	hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	
	hmm_algorithm._hg.load_kjd_array();
	p_TS->set( thread_id );

	while ( ! p_TS->is_Can_go() ) {
		usleep( 10 );
	}

	int n_hidden_states;

	hmm_algorithm.transmit_state();
	hmm_algorithm.calculate_errors( p_engine->mutation_rate );

	for ( size_t i = 0; i < iters.size(); i++ ) {
		int s_beg = iters[ i ].first;
		int s_end = iters[ i ].second;
		exclude_samples.clear();
		for ( int s = s_beg; s < s_end; s++ ) {
			exclude_samples[ p_so->orders[ s ] ] = 0;
		}
		if ( i != 0 ) {
			prev_samples.swap( next_samples );
			next_samples.clear();
			for ( int s = s_beg; s < s_end; s++ ) {
				next_samples.push_back( p_so->orders[ s ] );
			}
			hmm_algorithm._hg.free_kjd_array();
			hmm_algorithm._hg.modify_hg( p_engine->haplotypes, prev_samples, next_samples, blocks );
			hmm_algorithm._hg.load_kjd_array();
		}

		n_hidden_states = 2 * ( individuals - exclude_samples.size() );
		n_hidden_states = reference_only ? p_engine->n_reference : n_hidden_states + p_engine->n_reference;
		hmm_algorithm.set_hapNumber( n_hidden_states );	
		hmm_algorithm.alloc_leftmatrix();
		for ( int s = s_beg; s < s_end; s++ ) {
			hmm_algorithm.process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
			_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
			for ( int nc = 0; nc < sample_time; nc++ ) {
				hmm_algorithm.update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
								error_model, *p_uniform, crossover, p_engine->mutation_rate );
				GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
			}
			_Output_Wizzard OMW;
			OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
		}
		hmm_algorithm.free_leftmatrix();
	}
}

void * prune_start_function( void * arg )
{
	_Thread_Arg * t_arg = ( _Thread_Arg * ) arg;
	int individuals = t_arg->individuals;
	int snp_number  = t_arg->snp_number;
	int states      = t_arg->states;
	int begin	= t_arg->begin;
	int end		= t_arg->end;
	int start	= t_arg->start;
	int thread_id   = t_arg->thread_id;
	int maximum     = t_arg->maximum;
	int blk_size    = t_arg->blk_size;
	int indiv_ost   = t_arg->indiv_ost;
	int sample_time = t_arg->sample_time;
	bool burnin     = t_arg->burnin;
	bool gl         = t_arg->gl;
	bool reference_only = t_arg->reference_only;
	Error * error_model = t_arg->error_model;
	Theta * crossover   = t_arg->crossover;
	Init_Type * p_engine = t_arg->engine;
	_GenoHap * p_genohap = t_arg->genohap;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform = t_arg->p_uniform;
	_Thread_Syn * p_TS = t_arg->p_TS;
	_Sample_Orders * p_so = t_arg->p_so;

	_HMM_Algo hmm_algorithm( individuals, snp_number, start, indiv_ost );

	vector< pair< int, int > > & iters = t_arg->vec_iters;
	int sample_beg = iters[ 0 ].first;
	int sample_end = iters[ 0 ].second;
	map< int, int > exclude_samples;
	vector< int > prev_samples, next_samples;
	for ( int s = sample_beg; s < sample_end; s++ ) {
		exclude_samples[ p_so->orders[ s ] ] = 0;
		prev_samples.push_back( p_so->orders[ s ] );
		next_samples.push_back( p_so->orders[ s ] );
	}

	vector< int > blocks;
	if ( burnin ) {
		hmm_algorithm._hg.split_blocks( p_engine->haplotypes, states, p_engine->reference, p_engine->n_reference, *p_uniform, reference_only );
		hmm_algorithm._hg.list2vector( blocks );
		hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	} else {
		hmm_algorithm._hg.split_blocks( blocks, blk_size, *p_uniform );
		hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	}
	hmm_algorithm._hg.load_kjd_array();
	p_TS->set( thread_id );

	while ( ! p_TS->is_Can_go() ) {
		usleep( 10 );
	}

	int n_hidden_states;

	hmm_algorithm.transmit_state();
	hmm_algorithm.calculate_errors( p_engine->mutation_rate );

	for ( size_t i = 0; i < iters.size(); i++ ) {
		int s_beg = iters[ i ].first;
		int s_end = iters[ i ].second;
		exclude_samples.clear();
		for ( int s = s_beg; s < s_end; s++ ) {
			exclude_samples[ p_so->orders[ s ] ] = 0;
		}
		if ( i != 0 ) {
			prev_samples.swap( next_samples );
			next_samples.clear();
			for ( int s = s_beg; s < s_end; s++ ) {
				next_samples.push_back( p_so->orders[ s ] );
			}
			hmm_algorithm._hg.free_kjd_array();
			hmm_algorithm._hg.modify_hg( p_engine->haplotypes, prev_samples, next_samples, blocks );
			hmm_algorithm._hg.load_kjd_array();
		}
		
		if ( burnin ) {
			n_hidden_states = 2 * ( individuals - exclude_samples.size() );
			n_hidden_states = reference_only ? p_engine->n_reference : n_hidden_states + p_engine->n_reference;
			hmm_algorithm.set_hapNumber( n_hidden_states );
			hmm_algorithm.alloc_leftmatrix();
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
			}
			hmm_algorithm.free_leftmatrix();
		} else {
			n_hidden_states = p_engine->n_reference + 2 * ( individuals - exclude_samples.size() );
			hmm_algorithm.set_hapNumber( n_hidden_states );
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.prune_states_find( p_so->orders[ s ], p_engine->haplotypes, maximum );
				hmm_algorithm.prune_alloc_leftmatrix();
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				hmm_algorithm.prune_process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.prune_update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
				hmm_algorithm.prune_free_leftmatrix();
			}
		}
	}
}

void * hybrid_start_function( void * arg )
{
	_Thread_Arg * t_arg = ( _Thread_Arg * ) arg;
	int individuals = t_arg->individuals;
	int snp_number  = t_arg->snp_number;
	int states      = t_arg->states;
	int begin	= t_arg->begin;
	int end		= t_arg->end;
	int start	= t_arg->start;
	int thread_id   = t_arg->thread_id;
	int maximum     = t_arg->maximum;
	int blk_size    = t_arg->blk_size;
	int indiv_ost   = t_arg->indiv_ost;
	int sample_time = t_arg->sample_time;
	bool burnin     = t_arg->burnin;
	bool gl         = t_arg->gl;
	bool reference_only = t_arg->reference_only;
	bool odd        = t_arg->odd;
	Error * error_model = t_arg->error_model;
	Theta * crossover   = t_arg->crossover;
	Init_Type * p_engine = t_arg->engine;
	_GenoHap * p_genohap = t_arg->genohap;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform = t_arg->p_uniform;
	_Thread_Syn * p_TS = t_arg->p_TS;
	_Sample_Orders * p_so = t_arg->p_so;

	_HMM_Algo hmm_algorithm( individuals, snp_number, start, indiv_ost );

	vector< pair< int, int > > & iters = t_arg->vec_iters;
	int sample_beg = iters[ 0 ].first;
	int sample_end = iters[ 0 ].second;
	map< int, int > exclude_samples;
	vector< int > prev_samples, next_samples;
	for ( int s = sample_beg; s < sample_end; s++ ) {
		exclude_samples[ p_so->orders[ s ] ] = 0;
		prev_samples.push_back( p_so->orders[ s ] );
		next_samples.push_back( p_so->orders[ s ] );
	}

	vector< int > blocks;
	if ( burnin || odd ) {
		hmm_algorithm._hg.split_blocks( p_engine->haplotypes, states, p_engine->reference, p_engine->n_reference, *p_uniform, reference_only );
		hmm_algorithm._hg.list2vector( blocks );
		hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	} else {
		hmm_algorithm._hg.split_blocks( blocks, blk_size, *p_uniform );
		hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	}
	hmm_algorithm._hg.load_kjd_array();
	p_TS->set( thread_id );

	while ( ! p_TS->is_Can_go() ) {
		usleep( 10 );
	}

	int n_hidden_states;

	hmm_algorithm.transmit_state();
	hmm_algorithm.calculate_errors( p_engine->mutation_rate );

	for ( size_t i = 0; i < iters.size(); i++ ) {
		int s_beg = iters[ i ].first;
		int s_end = iters[ i ].second;
		exclude_samples.clear();
		for ( int s = s_beg; s < s_end; s++ ) {
			exclude_samples[ p_so->orders[ s ] ] = 0;
		}
		if ( i != 0 ) {
			prev_samples.swap( next_samples );
			next_samples.clear();
			for ( int s = s_beg; s < s_end; s++ ) {
				next_samples.push_back( p_so->orders[ s ] );
			}
			hmm_algorithm._hg.free_kjd_array();
			hmm_algorithm._hg.modify_hg( p_engine->haplotypes, prev_samples, next_samples, blocks );
			hmm_algorithm._hg.load_kjd_array();
		}
		
		if ( burnin || odd ) {
			n_hidden_states = 2 * ( individuals - exclude_samples.size() );
			n_hidden_states = reference_only ? p_engine->n_reference : n_hidden_states + p_engine->n_reference;
			hmm_algorithm.set_hapNumber( n_hidden_states );
			hmm_algorithm.alloc_leftmatrix();
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
			}
			hmm_algorithm.free_leftmatrix();
		} else {
			n_hidden_states = p_engine->n_reference + 2 * ( individuals - exclude_samples.size() );
			hmm_algorithm.set_hapNumber( n_hidden_states );
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.prune_states_find( p_so->orders[ s ], p_engine->haplotypes, maximum );
				hmm_algorithm.prune_alloc_leftmatrix();
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				hmm_algorithm.prune_process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.prune_update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
				hmm_algorithm.prune_free_leftmatrix();
			}
		}
	}
}

void * prune_normal_start_function( void * arg )
{
	_Thread_Arg * t_arg = ( _Thread_Arg * ) arg;
	int individuals = t_arg->individuals;
	int snp_number  = t_arg->snp_number;
	int states      = t_arg->states;
	int begin	= t_arg->begin;
	int end		= t_arg->end;
	int start	= t_arg->start;
	int thread_id   = t_arg->thread_id;
	int maximum     = t_arg->maximum;
	int blk_size    = t_arg->blk_size;
	int indiv_ost   = t_arg->indiv_ost;
	int sample_time = t_arg->sample_time;
	bool burnin     = t_arg->burnin;
	bool gl         = t_arg->gl;
	bool reference_only = t_arg->reference_only;
	Error * error_model = t_arg->error_model;
	Theta * crossover   = t_arg->crossover;
	Init_Type * p_engine = t_arg->engine;
	_GenoHap * p_genohap = t_arg->genohap;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform = t_arg->p_uniform;
	_Thread_Syn * p_TS = t_arg->p_TS;
	_Sample_Orders * p_so = t_arg->p_so;

	_HMM_Algo hmm_algorithm( individuals, snp_number, start, indiv_ost );

	vector< pair< int, int > > & iters = t_arg->vec_iters;
	int sample_beg = iters[ 0 ].first;
	int sample_end = iters[ 0 ].second;
	map< int, int > exclude_samples;
	vector< int > prev_samples, next_samples;
	for ( int s = sample_beg; s < sample_end; s++ ) {
		exclude_samples[ p_so->orders[ s ] ] = 0;
		prev_samples.push_back( p_so->orders[ s ] );
		next_samples.push_back( p_so->orders[ s ] );
	}

	vector< int > blocks;
	hmm_algorithm._hg.split_blocks( p_engine->haplotypes, states, p_engine->reference, p_engine->n_reference, *p_uniform, reference_only );
	hmm_algorithm._hg.list2vector( blocks );
	hmm_algorithm._hg.construct_hg( p_engine->haplotypes, exclude_samples, blocks, p_engine->reference, p_engine->n_reference, reference_only );
	hmm_algorithm._hg.load_kjd_array();
	p_TS->set( thread_id );

	while ( ! p_TS->is_Can_go() ) {
		usleep( 10 );
	}

	int n_hidden_states;

	hmm_algorithm.transmit_state();
	hmm_algorithm.calculate_errors( p_engine->mutation_rate );

	for ( size_t i = 0; i < iters.size(); i++ ) {
		int s_beg = iters[ i ].first;
		int s_end = iters[ i ].second;
		exclude_samples.clear();
		for ( int s = s_beg; s < s_end; s++ ) {
			exclude_samples[ p_so->orders[ s ] ] = 0;
		}
		if ( i != 0 ) {
			prev_samples.swap( next_samples );
			next_samples.clear();
			for ( int s = s_beg; s < s_end; s++ ) {
				next_samples.push_back( p_so->orders[ s ] );
			}
			hmm_algorithm._hg.free_kjd_array();
			hmm_algorithm._hg.modify_hg( p_engine->haplotypes, prev_samples, next_samples, blocks );
			hmm_algorithm._hg.load_kjd_array();
		}
		
		if ( burnin ) {
			n_hidden_states = 2 * ( individuals - exclude_samples.size() );
			n_hidden_states = reference_only ? p_engine->n_reference : n_hidden_states + p_engine->n_reference;
			hmm_algorithm.set_hapNumber( n_hidden_states );
			hmm_algorithm.alloc_leftmatrix();
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
			}
			hmm_algorithm.free_leftmatrix();
		} else {
			n_hidden_states = p_engine->n_reference + 2 * ( individuals - exclude_samples.size() );
			hmm_algorithm.set_hapNumber( n_hidden_states );
			for ( int s = s_beg; s < s_end; s++ ) {
				hmm_algorithm.prune_states_find( p_so->orders[ s ], p_engine->haplotypes, maximum );
				hmm_algorithm.prune_alloc_leftmatrix();
				_Gibbsgeno_Store GMI( 1, snp_number, sample_time );
				hmm_algorithm.prune_process( p_engine->theta_array, *p_genohap, p_so->orders[ s ] );
				for ( int nc = 0; nc < sample_time; nc++ ) {
					hmm_algorithm.prune_update_haplotype( p_engine->haplotypes, p_so->orders[ s ], *p_genohap, p_engine->theta_array, \
									error_model, *p_uniform, crossover, p_engine->mutation_rate );
					GMI.load( p_engine->haplotypes, p_so->orders[ s ], nc );
				}
				_Output_Wizzard OMW;
				OMW.decide_one( GMI, *p_engine, sample_time, p_so->orders[ s ], *p_uniform );
				hmm_algorithm.prune_free_leftmatrix();
			}
		}
	}
}
