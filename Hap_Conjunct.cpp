#include "Hap_Conjunct.h"

/* it will begin with snp_number - overlap / 2 */
void _Hap_Conjunct::load_prev( Init_Type & engine, int begin )
{
	if ( begin + overlap / 2 > engine.snp_number ) {
		return;
	}
	for ( int i = 0; i < individuals; i++ ) {
		if ( orders[ i ] ) {
			memcpy( prev_haplotypes[ i * 2 + 0 ].haplotype, engine.haplotypes[ i * 2 + 0 ].haplotype + begin, overlap / 2 );
			memcpy( prev_haplotypes[ i * 2 + 1 ].haplotype, engine.haplotypes[ i * 2 + 1 ].haplotype + begin, overlap / 2 );
		} else {
			memcpy( prev_haplotypes[ i * 2 + 0 ].haplotype, engine.haplotypes[ i * 2 + 1 ].haplotype + begin, overlap / 2 );
			memcpy( prev_haplotypes[ i * 2 + 1 ].haplotype, engine.haplotypes[ i * 2 + 0 ].haplotype + begin, overlap / 2 );
		}
	}
}

/* it will begin with overlap / 2 */
void _Hap_Conjunct::load_next( Init_Type & engine, int begin )
{
	for ( int i = 0; i < individuals; i++ ) {
		memcpy( next_haplotypes[ i * 2 + 0 ].haplotype, engine.haplotypes[ i * 2 + 0 ].haplotype + begin, overlap / 2 );
		memcpy( next_haplotypes[ i * 2 + 1 ].haplotype, engine.haplotypes[ i * 2 + 1 ].haplotype + begin, overlap / 2 );
	}
}

void _Hap_Conjunct::show()
{
	for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
		cerr << "prev indiv " << indiv_id << " haplotype 0 |";
		for ( int snp_id = 0; snp_id < overlap / 2; snp_id++ ) {
			cerr << int(prev_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ]);
		}
		cerr << "\n";
		cerr << "prev indiv " << indiv_id << " haplotype 1 |";
		for ( int snp_id = 0; snp_id < overlap / 2; snp_id++ ) {
                        cerr << int(prev_haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ]);
                }
                cerr << "\n";
		cerr << "next indiv " << indiv_id << " haplotype 0 |";
                for ( int snp_id = 0; snp_id < overlap / 2; snp_id++ ) {
                        cerr << int(next_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ]);
                }
                cerr << "\n";
                cerr << "next indiv " << indiv_id << " haplotype 1 |";
                for ( int snp_id = 0; snp_id < overlap / 2; snp_id++ ) {
                        cerr << int(next_haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ]);
                }
                cerr << "\n";
	}
}

void _Hap_Conjunct::update_orders()
{
	for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
		for ( int snp_id = 0; snp_id < overlap / 2; snp_id++ ) {
			if ( prev_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ] + prev_haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ] != \
			     next_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ] + next_haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ] )
			{
				continue;
			} else {
				if ( 1 == prev_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ] + prev_haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ] ) {
					if ( prev_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ] == next_haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ] ) {
						orders[ indiv_id ] = true;
						break;
					} else {
						orders[ indiv_id ] = false;
						break;
					}
				}
			}
		}
	}
}
