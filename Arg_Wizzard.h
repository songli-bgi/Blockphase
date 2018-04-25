#ifndef __ARG_WIZZARD_H__

#define __ARG_WIZZARD_H__

#include "anyarg.h"
#include "options.h"

struct _Program_Variables
{
	int n_threads;
	int begin;
	int end;
	int burnin;
	int circle;
	int states;

	int piece;
	int overlap;
	int eml;
	int maximum;
	int threshold;
	int sample_time;
	int subChunk;

	float mutation_rate;
	float recombination_rate;

	bool random;
	bool ibd;
	bool em;
	bool normal;
	bool prune;
	bool hybrid;
	
	bool gl;
	bool hap;
	bool reference_only;
	
	string in;
	string out;
	string reference;

	unsigned long seed;

	_Program_Variables()
	{
		n_threads = 1;
		begin = 0; // later determination
		end   = 0; // later determination
		burnin = 10;
		circle = 30;
		states = 30;
		eml = 3;
	
		piece = 10000;
		overlap = 500;
		maximum = 60;
		threshold = 100;
		sample_time = 1;
		subChunk = 4;

		mutation_rate = 0.01;
		recombination_rate = 0.01;

		random = true;	
		ibd = true;
		em = false;
		gl = false;
		hap = true;
		normal = false;
		prune = false;
		hybrid = true;
		reference_only = false;
		
		seed = 9999;
	}
	
	void show_variables( FILE * fp );
};

struct _Arg_Wizzard
{
	static void init( anyarg & opt );
	static void process( anyarg & opt, _Program_Variables & pv, _Options & OPH );
};

#endif
