#include "Output_Wizzard.h"

void _Output_Wizzard::phased_out( _Gibbsgeno_Store & GS, Init_Type & engine, int out_begin, int out_end, bool add, \
				  _Run_Wizzard & RW, int start, int cycle )
{
	gzFile fp;
	string phased_file = outfile + ".phased.gz";
	if ( add ) {
		fp = gzopen( phased_file.c_str(), "a" );
	} else {
		fp = gzopen( phased_file.c_str(), "w" );
		gzprintf( fp, "position\talleleA\talleleB" );
		for ( int i = 0; i < RW.samples.size(); i++ ) {
			gzprintf( fp, "\t%s\t%s", RW.samples[ i ].c_str(), RW.samples[ i ].c_str() );
		}
		gzprintf( fp, "\n" );
		
	}
	this->decide_all(GS, engine, cycle);
	if ( add ) {
		_HC->load_next( engine, out_begin );
		_HC->update_orders();
		_HC->load_prev( engine, out_end );
	} else {
		_HC->load_prev( engine, out_end );
	}
	for ( int snp_id = out_begin; snp_id < out_end; snp_id++ ) {
		gzprintf( fp, "%d\t%c\t%c", RW.snpids[ start + snp_id ], RW.ref_Alt_bases[ start + snp_id ][ 0 ], RW.ref_Alt_bases[ start + snp_id ][ 1 ] );
		for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
			int base1( -1 ), base2( -1 );
			bool spin = true;
			if ( add ) {
				if( !_HC->orders[ indiv_id ] ) {
					spin = false;
				}
			}
			if ( spin ) {
				base1 = engine.haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ];
				base2 = engine.haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ];
			} else {
				base1 = engine.haplotypes[ indiv_id * 2 + 1 ].haplotype[ snp_id ];
				base2 = engine.haplotypes[ indiv_id * 2 + 0 ].haplotype[ snp_id ];
			}
			gzprintf( fp, "\t%c\t%c", RW.ref_Alt_bases[ start + snp_id ][ base1 ], RW.ref_Alt_bases[ start + snp_id ][ base2 ] );
		}
		gzprintf( fp, "\n" );
	}
	gzclose( fp );
}

int _Output_Wizzard::max_id( uint8_t * g )
{
	int ret = -1;
        int t = 0;
        for ( int i = 0; i < 3; i++ ) {
                if ( g[ i ] > t ) {
                        t = g[ i ];
                        ret = i;
                }
        }
        return ret;
}

void _Output_Wizzard::gprobs_out( _Gibbsgeno_Store & GS, int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start )
{
	gzFile fp;
	string gprobs_file = outfile + ".gprobs.gz";
	if ( add ) {
		fp = gzopen( gprobs_file.c_str(), "a" );
	} else {
		fp = gzopen( gprobs_file.c_str(), "w" );
                gzprintf( fp, "position" );
                for ( int i = 0; i < RW.samples.size(); i++ ) {
                        gzprintf( fp, "\t%s\t%s\t%s", RW.samples[ i ].c_str(), RW.samples[ i ].c_str(), RW.samples[ i ].c_str() );
                }
                gzprintf( fp, "\n" );
	}

	int gibbs_genotypes[ 3 ];
	int base1, base2;

	for ( int snp_id = out_begin; snp_id < out_end; snp_id++ ) {
		gzprintf( fp, "%d", RW.snpids[ start + snp_id ] );
		for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
			memset( gibbs_genotypes, 0, sizeof( int ) * 3 );
			for ( int c = 0; c < GS.circles; c++ ) {
				base1 = GS.get( snp_id, c, 0, indiv_id );
				base2 = GS.get( snp_id, c, 1, indiv_id );
				gibbs_genotypes[ base1 + base2 ]++;
			}
			double sum = 0;
			for ( int i = 0; i < 3; i++ ) {
				sum += double( gibbs_genotypes[ i ] );
			}
			for ( int i = 0; i < 3; i++ ) {
				gzprintf( fp, "\t%f", double( gibbs_genotypes[ i ] ) / sum );
			}
		}
		gzprintf( fp, "\n" );
	}
	gzclose( fp );
}

void _Output_Wizzard::dosage_out( _Gibbsgeno_Store & GS, \
				  int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start )
{
	sq_Pearson_cc sPc;
	gzFile fp;
	string dosage_file = outfile + ".dosage.gz";
	if ( add ) {
		fp = gzopen( dosage_file.c_str(), "a" );
	} else {
		fp = gzopen( dosage_file.c_str(), "w" );
		gzprintf( fp, "position" );
		for ( int i = 0; i < RW.samples.size(); i++ ) {
			gzprintf( fp, "\t%s", RW.samples[ i ].c_str() );
		}
		gzprintf( fp, "\n" );
	}

	int gibbs_genotypes[ 3 ];
	int base1, base2;

	for ( int snp_id = out_begin; snp_id < out_end; snp_id++ ) {
		gzprintf( fp, "%d", RW.snpids[ start + snp_id ] );
		for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
			memset( gibbs_genotypes, 0, sizeof( int ) * 3 );
			for ( int c = 0; c < GS.circles; c++ ) {
				base1 = GS.get( snp_id, c, 0, indiv_id );
				base2 = GS.get( snp_id, c, 0, indiv_id );
				gibbs_genotypes[ base1 + base2 ]++;
			}
			double g_dosage = sPc.cal_dosage( gibbs_genotypes, indiv_id, snp_id );
			gzprintf( fp, "\t%f", g_dosage );
		}
		gzprintf( fp, "\n" );
	}
	gzclose( fp );
}

void _Output_Wizzard::r2_out( _Gibbsgeno_Store & GS, Init_Type & engine, \
			      int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start )
{
	sq_Pearson_cc sPc;
	gzFile fp;
	string r2_file = outfile + ".r2.gz";
	if ( add ) {
		fp = gzopen( r2_file.c_str(), "a" );
	} else {
		fp = gzopen( r2_file.c_str(), "w" );
		gzprintf( fp, "position\tgibbs.r2\n" );
	}

	int ** gibbs_genotypes = new int * [ individuals ];
	for ( int i = 0; i < individuals; i++ ) {
		gibbs_genotypes[ i ] = new int [ 3 ];
	}

	for ( int snp_id = out_begin; snp_id < out_end; snp_id++ ) {
		gzprintf( fp, "%d", RW.snpids[ start + snp_id ] );
		int base1, base2;
		for ( int indiv_id = 0; indiv_id < individuals; indiv_id++ ) {
			memset( gibbs_genotypes[ indiv_id ], 0, sizeof( int ) * 3 );
			for ( int c = 0; c < GS.circles; c++ ) {
				base1 = GS.get( snp_id, c, 0, indiv_id );
				base2 = GS.get( snp_id, c, 1, indiv_id );
				gibbs_genotypes[ indiv_id ][ base1 + base2 ]++;
			}
		}
		double g_r2 = sPc.cal_r_squared( engine, gibbs_genotypes, snp_id );
		gzprintf( fp, "\t%f\n", g_r2 );
	}
	gzclose( fp );

	for ( int i = 0; i < individuals; i++ ) {
		delete [] gibbs_genotypes[ i ];
	}
	delete [] gibbs_genotypes;
}

void _Output_Wizzard::process( _Gibbsgeno_Store & GS, Init_Type & engine, \
			       int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start, int cycle )
{
	phased_out( GS, engine, out_begin, out_end, add, RW, start, cycle );
	gprobs_out( GS, out_begin, out_end, add, RW, start );
	dosage_out( GS, out_begin, out_end, add, RW, start );
	r2_out( GS, engine, out_begin, out_end, add, RW, start );
}

void _Output_Wizzard::decide_all(_Gibbsgeno_Store & GS, Init_Type & engine, int cycle)
{
	vector<int> hete_idx;
	for (int i = 0; i != engine.individuals; ++i) {
		this->decide_homo(i, GS, engine, cycle, hete_idx);
		this->decide_hete(i, GS, engine, cycle, hete_idx);
	}
}

void _Output_Wizzard::decide_one( _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, int individuals_ID, \
				  variate_generator< kreutzer1986, uniform_real<> > & uniform )
{
	vector<int> hete_idx;
	decide_homo(individuals_ID, GS, engine, cycle, hete_idx, 0, uniform);
	decide_hete(individuals_ID, GS, engine, cycle, hete_idx, 0);
}

void _Output_Wizzard::decide_homo(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx, int cur_id, \
				  variate_generator< kreutzer1986, uniform_real<> > & uniform)
{
	hete_idx.clear();
        map<int, int> best_gt;
        int geno, appe; // diplotype & appearance
        for (int j = 0; j != engine.snp_number; ++j) {
                best_gt.clear();
                geno = appe = 0;
                for (int k = 0; k != cycle; ++k) {
			int base1 = GS.get( j, k, 0, cur_id );
			int base2 = GS.get( j, k, 1, cur_id );
                        best_gt[base1 + base2] += 1;
                }
                for (map<int, int>::iterator it = best_gt.begin(); it != best_gt.end(); ++it) {
                        if (it->second == appe) {
                                int rand_n = uniform() > 0.5 ? 0 : 1;
                                if (rand_n) {
                                        geno = it->first;
                                        appe = it->second;
                                }
                        } else if (it->second > appe) {
                                geno = it->first;
                                appe = it->second;
                        } else { continue; }
                }
                if (geno == 0 || geno == 2)
                        engine.haplotypes[ i * 2 + 0 ].haplotype[ j ] = engine.haplotypes[ i * 2 + 1 ].haplotype[ j ] = geno / 2;
                else
                        hete_idx.push_back(j);
        }
}

void _Output_Wizzard::decide_homo(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx)
{
	hete_idx.clear();
        map<int, int> best_gt;
        int geno, appe; // diplotype & appearance
        for (int j = 0; j != engine.snp_number; ++j) {
                best_gt.clear();
                geno = appe = 0;
                for (int k = 0; k != cycle; ++k) {
			int base1 = GS.get( j, k, 0, i );
			int base2 = GS.get( j, k, 1, i );
                        best_gt[base1 + base2] += 1;
                }
                for (map<int, int>::iterator it = best_gt.begin(); it != best_gt.end(); ++it) {
                        if (it->second == appe) {
                                int rand_n = rand() % 2;
                                if (rand_n) {
                                        geno = it->first;
                                        appe = it->second;
                                }
                        } else if (it->second > appe) {
                                geno = it->first;
                                appe = it->second;
                        } else { continue; }
                }
                if (geno == 0 || geno == 2)
                        engine.haplotypes[ i * 2 + 0 ].haplotype[ j ] = engine.haplotypes[ i * 2 + 1 ].haplotype[ j ] = geno / 2;
                else
                        hete_idx.push_back(j);
        }
}

void _Output_Wizzard::decide_hete(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx, int cur_id)
{
	if (hete_idx.size() <= 1) {
                for (vector<int>::iterator it = hete_idx.begin(); it != hete_idx.end(); ++it) {
                        engine.haplotypes[ i * 2 + 0 ].haplotype[*it] = 0;
                        engine.haplotypes[ i * 2 + 1 ].haplotype[*it] = 1;
                }
                return;
        }
        vector<int>::iterator it_a, it_b;
        it_b = hete_idx.begin(); it_a = it_b++;
        bool Turn = false; // true if previous hete site is 1/0, false if 0/1
        for (; it_b != hete_idx.end(); ++it_a, ++it_b) {
                int phase[2] = {0, 0}; // 0 for cis-phase, 1 for trans-phase
                for (int k = 0; k != cycle; ++k) {
			int base11 = GS.get( *it_a, k, 0, cur_id );
			int base12 = GS.get( *it_a, k, 1, cur_id );
			int base21 = GS.get( *it_b, k, 0, cur_id );
			int base22 = GS.get( *it_b, k, 1, cur_id );
                        if ( base11 == base21 ) {
                                phase[0]++;
                        } else {
                                phase[1]++;
                        }
                        if ( base12 == base22 ) {
                                phase[0]++;
                        } else {
                                phase[1]++;
                        }
                }
                if (it_a == hete_idx.begin()) {
                        engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_a ] = 0;
                        engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_a ] = 1;
                }
                if (phase[0] >= phase[1]) { // cis-phase 00/11
                        if (!Turn) {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 0;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 1;
                        } else {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 1;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 0;
                        }
                } else {
                        if (Turn) {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 0;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 1;
                                Turn = false;
                        } else {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 1;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 0;
                                Turn = true;
                        }
                }
        }
}

void _Output_Wizzard::decide_hete(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx)
{
	if (hete_idx.size() <= 1) {
                for (vector<int>::iterator it = hete_idx.begin(); it != hete_idx.end(); ++it) {
                        engine.haplotypes[ i * 2 + 0 ].haplotype[*it] = 0;
                        engine.haplotypes[ i * 2 + 1 ].haplotype[*it] = 1;
                }
                return;
        }
        vector<int>::iterator it_a, it_b;
        it_b = hete_idx.begin(); it_a = it_b++;
        bool Turn = false; // true if previous hete site is 1/0, false if 0/1
        for (; it_b != hete_idx.end(); ++it_a, ++it_b) {
                int phase[2] = {0, 0}; // 0 for cis-phase, 1 for trans-phase
                for (int k = 0; k != cycle; ++k) {
			int base11 = GS.get( *it_a, k, 0, i );
			int base12 = GS.get( *it_a, k, 1, i );
			int base21 = GS.get( *it_b, k, 0, i );
			int base22 = GS.get( *it_b, k, 1, i );
                        if ( base11 == base21 ) {
                                phase[0]++;
                        } else {
                                phase[1]++;
                        }
                        if ( base12 == base22 ) {
                                phase[0]++;
                        } else {
                                phase[1]++;
                        }
                }
                if (it_a == hete_idx.begin()) {
                        engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_a ] = 0;
                        engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_a ] = 1;
                }
                if (phase[0] >= phase[1]) { // cis-phase 00/11
                        if (!Turn) {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 0;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 1;
                        } else {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 1;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 0;
                        }
                } else {
                        if (Turn) {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 0;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 1;
                                Turn = false;
                        } else {
                                engine.haplotypes[ i * 2 + 0 ].haplotype[ *it_b ] = 1;
                                engine.haplotypes[ i * 2 + 1 ].haplotype[ *it_b ] = 0;
                                Turn = true;
                        }
                }
        }
}
