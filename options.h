#ifndef __OPTIONS_H__

#define __OPTIONS_H__

#include <stdio.h>

struct _Options
{
	void help()
	{
		fprintf( stderr, "Block Phase :\n" );
		fprintf( stderr, "  * Authors : Song Li, Jin Wei\n" );
		fprintf( stderr, "  * Contact : \033[0;36mli.song@bgi.com\033[0m\n" );
		fprintf( stderr, "  * Version : 1.0\n\n" );
		fprintf( stderr, "Description :\n" );
		fprintf( stderr, "  * Chromosomes : Autosomes and X chromosome\n" );
		fprintf( stderr, "  * Features    : Gibbs Sampling\n" );
		fprintf( stderr, "  * Functions   : Imputation and Phasing for Unrelated Individuals\n\n" );
		fprintf( stderr, "Options :\n" );
		fprintf( stderr, "  * Parameters  :\n" );
		fprintf( stderr, "    \033[0;36m-b   --burnin\033[0m                   Burnin iteratings [10]\n" );
		fprintf( stderr, "    \033[0;36m-c   --circle\033[0m                   Main iteratings [30]\n" );
		fprintf( stderr, "    \033[0;36m-s   --states\033[0m                   Hidden states in burnin step [30]\n" );
		fprintf( stderr, "    \033[0;36m-t   --n_threads\033[0m                Number of threads [1]\n" );
		fprintf( stderr, "         \033[;36m--begin\033[0m                    Physical marker to begin with [first]\n" );
		fprintf( stderr, "         \033[;36m--end\033[0m                      Physical marker to end with [last]\n" );
		fprintf( stderr, "         \033[;36m--chunk\033[0m                    Number of markers in one chunk [10000]\n" );
		fprintf( stderr, "         \033[;36m--share\033[0m                    Shared markers between adjacent chunks [500]\n" );
		fprintf( stderr, "         \033[;36m--maximum\033[0m                  Maximum hidden states in main step [60]\n" );
		fprintf( stderr, "         \033[;36m--sample_time\033[0m              Number of backward sampling [1]\n" );
		fprintf( stderr, "         \033[;36m--subChunk\033[0m                 Number of sub chunks in each thread [4]\n" );
		fprintf( stderr, "  * Initialize  :\n" );
		fprintf( stderr, "    \033[0;36m-R   --recombination_rate\033[0m       Initialize recombination rate [0.01]\n" );
		fprintf( stderr, "    \033[0;36m-m   --mutation_rate\033[0m            Initialize mutation rate [0.01]\n" );
		fprintf( stderr, "         \033[;36m--random\033[0m                   Initialize haplotypes in random [true]\n");
		fprintf( stderr, "         \033[;36m--em\033[0m                       Initialize haplotypes in EM strategy [false]\n" );
		fprintf( stderr, "         \033[;36m--eml\033[0m                      Length of local haplotypes in EM strategy [3]\n");
		fprintf( stderr, "         \033[;36m--ibd\033[0m                      Initialize haplotypes using IBD information [true]\n" );
		fprintf( stderr, "         \033[;36m--threshold\033[0m                Maximum number of homozygote segments for IBD \\\n                                    initializing, it goes with \033[;36m--ibd\033[0m argument [100]\n" );
		fprintf( stderr, "         \033[;36m--seed\033[0m                     Random seed [9999]\n\n" );
		fprintf( stderr, "Reference :\n" );
		fprintf( stderr, "    \033[0;36m-r   --reference\033[0m                Reference panel file, compressed or not\n" );
		fprintf( stderr, "         \033[;36m--reference_only\033[0m           Only use reference panel to construct haplotype \\\n                                    graph in burnin step, it goes with \033[0;36m--reference\033[0m \\\n                                    argument [false]\n\n" );
		fprintf( stderr, "Input/Output :\n" );
		fprintf( stderr, "    \033[0;36m-i   --in\033[0m                       Input file, compressed or not\n" );
		fprintf( stderr, "    \033[0;36m-o   --out\033[0m                      Output file, compressed files\n\n" );
		fprintf( stderr, "Modules :\n" );
		fprintf( stderr, "    \033[0;36m-N   --normal\033[0m                   Use step out strategy to split chromosome [false]\n" );
		fprintf( stderr, "    \033[0;36m-P   --prune\033[0m                    Use increasing length strategy to split chromosome [false]\n" );
		fprintf( stderr, "    \033[0;36m-H   --hybrid\033[0m                   Use step out and increasing length strategy alternately [true]\n\n" );
		fprintf( stderr, "         \033[0;36m--gl\033[0m                       Choose hidden states by GL [false]\n" );
		fprintf( stderr, "         \033[0;36m--hap\033[0m                      Choose hidden states by IBD [true]\n\n" );
		fprintf( stderr, "Others :\n" );
		fprintf( stderr, "    \033[0;36m-h   --help\033[0m                     Help information\n" );
		fprintf( stderr, "    \033[0;36m-v   --version\033[0m                  Version\n\n" );
	}
};

#endif
