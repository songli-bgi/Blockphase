#include "Arg_Wizzard.h"
#include <assert.h>
#include <time.h>

void _Program_Variables::show_variables( FILE * fp )
{
	fprintf( fp, "Block Phase :\n");
	fprintf( fp, "  * Authors : Song Li, Jin Wei, Li Qibin\n" );
	fprintf( fp, "  * Contact : songli3@bgitechsolutions.org.cn\n" );
	fprintf( fp, "  * Version : 1.0.raw\n" );
	time_t timep;
	struct tm *p;
	time( &timep );
	p = localtime( &timep );
	fprintf( fp, "  * Date    : ");
	fprintf( fp, "%0.2d/%0.2d/%d ", p->tm_mday, ( 1 + p->tm_mon), (1900 + p->tm_year) );
	fprintf( fp, "%0.2d:%0.2d:%0.2d\n\n", p->tm_hour, p->tm_min, p->tm_sec );

	fprintf( fp, "MODE :\n");
	fprintf( fp, "  * Autosome( chr1 ... chr22 ) & chrX\n" );
	fprintf( fp, "  * Window phasing and genotypes inferencing\n" );
	fprintf( fp, "  * Gibbs Sampling\n\n" );

	fprintf( fp, "Parameters :\n" );
	fprintf( fp, "  * Seed : %ld\n", seed );
	fprintf( fp, "  * Parallelisation : %d thread(s)\n", n_threads );
	fprintf( fp, "  * MCMC : %d burnin, %d circle, %d maximum\n", burnin, circle, maximum );
	fprintf( fp, "  * Model : %d states, %d markers in each chunk, %d shared, sample %d time(s)\n", states, piece, overlap, sample_time );
	fprintf( fp, "  * Initailize : random : %s, ibd %s, em %s, reference_only %s\n\n", random?"true":"false", ibd?"true":"false", em?"true":"false",\
			  reference_only?"true":"false" );

	fprintf( fp, "Modules :\n" );
	fprintf( fp, "  * normal : %s, prune : %s, hybrid : %s\n", normal?"true":"false", prune?"true":"false", hybrid?"true":"false" );
	fprintf( fp, "  * prune model : %s\n\n", gl?"GL":"IBD" );
}

void _Arg_Wizzard::init( anyarg & opt )
{
	opt.set_option( "n_threads", 't', "INT", "number of threads, default(1)" );
	opt.set_option( "begin", "INT", "set program begin marker, default( the first marker )" );
	opt.set_option( "end", "INT", "set program end marker, default( the last marker )" );
	opt.set_option( "burnin", 'b', "INT", "set burnin iterators, default(10)" );
	opt.set_option( "circle", 'c', "INT", "set main iterators, default(30)" );
	opt.set_option( "states", 's', "INT", "set conditioning states, default(30)" );
	opt.set_option( "eml", "INT", "set the length of local haplotype for EM initailizing" );

	opt.set_option( "chunk", "INT", "set number of snps for echo piece at low memory module, default(10000)" );
	opt.set_option( "share", "INT", "set number of sites at neighbor pieces, default(400)" );
	opt.set_option( "maximum", "INT", "Maximum hidden states in prune step, default(60)" );
	opt.set_option( "threshold", "INT", "set maximun number of homozygote segments for initializing, default(100)" );

	opt.set_option( "ibd", "INT", "IBD haplotypes, non-zero for using, default(non-zero)" );
	opt.set_option( "sample_time", "INT", "Sampling time for backward algorithm" );
	opt.set_option( "subChunk", "INT", "number of sub chunk for each thread" );
	
	opt.set_option( "in", 'i', "STRING", "input file, can be gz or txt" );
	opt.set_option( "out", 'o', "STRING", "output file, be gz" );
	opt.set_option( "reference", 'r', "STRING", "provide reference panel" );
	
	opt.set_option( "mutation_rate", 'm', "FLOAT", "initialise mutation rate, default(0.01)" );
	opt.set_option( "recombination_rate", 'R', "FLOAT", "initialise sites recombination rate, default(0.01)" );

	opt.set_option( "seed", "ULONG", "initialise random seed, default(9999)" );

	// set flag
	opt.set_flag( "random", "random initialize haplotype" );
	opt.set_flag( "em", "EM initialize haplotype" );

	opt.set_flag( "reference_only", "only use reference panel to initialize " );
	opt.set_flag( "prune", 'P', "prune module flag" );
	opt.set_flag( "hybrid", 'H', "hybrid module flag" );
	opt.set_flag( "normal", 'N', "normal module flag" );

	opt.set_flag( "gl", "choose hidden states by gl" );
	opt.set_flag( "hap", "choose hidden states by ibd" );

	opt.set_flag( "help", 'h', "show this help" );
	opt.set_flag( "version", 'v', "show version" );
}

void _Arg_Wizzard::process( anyarg & opt, _Program_Variables & pv, _Options & OPH )
{
	const char *s;
	if ( NULL != (s = opt.get_value("n_threads")) ) pv.n_threads = atoi( s );
	if ( NULL != (s = opt.get_value("begin")) ) pv.begin = atoi( s );
	if ( NULL != (s = opt.get_value("end")) ) pv.end = atoi( s );
	if ( NULL != (s = opt.get_value("burnin")) ) pv.burnin = atoi( s );
	if ( NULL != (s = opt.get_value("circle")) ) pv.circle = atoi( s );
	if ( NULL != (s = opt.get_value("states")) ) pv.states = atoi( s );
	if ( NULL != (s = opt.get_value("chunk")) ) pv.piece = atoi( s );
	if ( NULL != (s = opt.get_value("share")) ) pv.overlap = atoi( s );
	if ( NULL != (s = opt.get_value("maximum")) ) pv.maximum = atoi( s );
	if ( NULL != (s = opt.get_value("threshold")) ) pv.threshold = atoi( s );
	if ( NULL != (s = opt.get_value("sample_time")) ) pv.sample_time = atoi( s );
	if ( NULL != (s = opt.get_value("eml")) ) pv.eml = atoi( s );
	if ( NULL != (s = opt.get_value("subChunk")) ) pv.subChunk = atoi( s );

	assert( pv.n_threads > 0 );

	if ( NULL != (s = opt.get_value("in")) ) {
		pv.in = s;
	}

	if ( NULL != (s = opt.get_value("out")) ) {
		pv.out = s;
	}

	if ( NULL != (s = opt.get_value("reference")) ) {
                pv.reference = s;
	}

	if ( pv.reference.empty() ) {
		pv.reference_only = false;
	}

	if ( NULL != (s = opt.get_value("mutation_rate")) ) pv.mutation_rate = atof( s );
	if ( NULL != (s = opt.get_value("recombination_rate")) ) pv.recombination_rate = atof( s );
	
	if ( NULL != (s = opt.get_value("seed")) ) pv.seed = (unsigned long)(atoi( s ));

	if ( NULL != (s = opt.get_value("ibd")) ) {
		if ( 0 == atoi( s ) ) {
			pv.ibd = false;
		} else {
			pv.ibd = true;
		}
	}

	// get flag
	if ( opt.get_flag( "random" ) ) {
		pv.random = true;
		pv.em = false;
	}
	
	if ( opt.get_flag( "em" ) ) {
		pv.em = true;
		pv.random = false;
	}

	if ( opt.get_flag( "reference_only" ) ) {
		if ( pv.reference.empty() ) {
			fprintf( stderr, "Warning : No reference panel, reference_only is set to false.\n" );
			pv.reference_only = false;
		} else {
			pv.reference_only = true;
		}
	}
	
	if ( opt.get_flag( "normal" ) ) {
		pv.normal = true;
		pv.hybrid = false;
	}
	
	if ( opt.get_flag( "prune" ) ) {
		pv.prune = true;
		pv.hybrid = false;
	}
	
	if ( opt.get_flag( "hybrid" ) ) {
		pv.hybrid = true;
		pv.normal = false;
	}

	if ( opt.get_flag( "gl" ) ) {
		pv.gl = true;
		pv.hap = false;
	}

	if ( opt.get_flag( "hap" ) ) {
		pv.hap = true;
		pv.gl = false;
	}

	if ( opt.get_flag("version") ) {
		fprintf( stderr, "blockphase version 1.0\n" );
		exit( 0 );
	}

	if ( opt.get_flag("help") ) {
		OPH.help();
		exit( 0 );
	}
}
