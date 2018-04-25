#include "Hap_Init.h"
#include <pthread.h>

bool node_t::conjunct( const node_t & _n )
{
	int pesudo_start = _n.start > start ? start : _n.start;
	int pesudo_end   = _n.end   < end   ? end   : _n.end;
#ifdef DEBUG
	cerr << "_n.start = " << _n.start << "\t_n.end = " << _n.end << endl;
	cerr << "start = " << start << "\tend = " << end << endl;
	cerr << "pesudo_start = " << pesudo_start << "\tpesudo_end = " << pesudo_end << endl;
#endif
	int * p = new int [ pesudo_end - pesudo_start + 1 ];
	for ( int i = _n.start; i <= _n.end; i++ ) {
		p[ i - pesudo_start ] = _n.piece[ i - _n.start ];
	}
	for ( int i = start; i <= end; i++ ) {
		p[ i - pesudo_start ] = piece[ i - start ];
	}
	delete [] piece;
	piece = p;
	start = pesudo_start;
	end   = pesudo_end;
	length = pesudo_end - pesudo_start + 1;
}

void _Hap_Init::haplotype_init( vector< vector< int > > & genotypes, hap_t * haplotypes )
{
	snp_number = genotypes[ 0 ].size();
	individuals = genotypes.size();
	int gap = 100;
	int distance = 10000;
	vector< pair< int, int > > pieces;
	int marker_start = 0;
	int marker_end   = 0;
	int snp_index = 0;
	while ( snp_index < snp_number ) {
		if ( snp_index == 0 ) {
			marker_start = snp_index;
			snp_index += distance;
			if ( snp_index >= snp_number ) {
				marker_end = snp_number - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
				break;
			} else {
				marker_end = snp_index - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
			}
		} else {
			snp_index -= gap;
			marker_start = snp_index;
			snp_index += distance;
			if ( snp_index >= snp_number ) {
				marker_end = snp_number - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
				break;
			} else {
				marker_end = snp_index - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
			}
		}
	}
	
	bool start_from_zero = true;
	
	for ( int i = 0; i < pieces.size(); i++ ) {
		marker_start = pieces[ i ].first;
		marker_end   = pieces[ i ].second;
		haplotype_init_piece( genotypes, haplotypes, marker_start, marker_end, gap, start_from_zero );
		start_from_zero = false;
	}
}

void _Hap_Init::haplotype_init_piece( vector< vector< int > > & genotypes, hap_t * haplotypes, int marker_start, int marker_end, \
						int gap, bool start_from_zero )
{
	snp_number = marker_end - marker_start + 1;
	individuals = genotypes.size();
	// alloc a matrix for store hom-genotype segments.
	int16_t *** matrix = new int16_t ** [ snp_number ];
	for ( int i = 0; i < snp_number; i++ ) {
		matrix[ i ] = new int16_t * [ snp_number ];
		for ( int j = 0; j < snp_number; j++ ) {
			matrix[ i ][ j ] = NULL;
		}
	}

	uint8_t ** mark = new uint8_t * [ individuals ];
	for ( int i = 0; i < individuals; i++ ) {
		mark[ i ] = new uint8_t [ snp_number ];
		memset( mark[ i ], 0, snp_number );
	}

	{
		vector< vector < int > > groups;
		for ( int i = 0; i < individuals; i++ ) {
			vector< int > group_t;
			for ( int j = 0; j < snp_number; j++ ) {
				if ( genotypes[ i ][ j + marker_start ] != 1 ) {
					group_t.push_back( j );
				}
			}
			groups.push_back( group_t );
		}		

		for ( int i = 0; i < individuals; i++ ) {
			int n = groups[ i ].size();
			for ( int j = 0; j < n; j++ ) {
				mark[ i ][ groups[ i ][ j ] ] = 1;
			}
		}
	}
#ifdef DEBUG
	for ( int i = 0; i < individuals; i++ ) {
		cerr << "B" << i << " show:\t|";
		for ( int j = 0; j < snp_number; j++ ) {
			if ( mark[ i ][ j ] == 1 ) {
				cerr << ".";
			} else {
				cerr << " ";
			}
		}
		cerr << "|\n";
	}
#endif
	int hold = 5;

	for ( int i = 0; i < individuals; i++ ) {
		vector< int > vec_t;
		bool first = true;
		for ( int k = 0; k < snp_number; k++ ) {
			if ( first && mark[ i ][ k ] == 1 ) {
				vec_t.push_back( k );
				first = false;
				if ( k == snp_number - 1 ) {
					if ( vec_t.size() >= hold ) {
						if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
							int16_t * p_sample = new int16_t [ threshold ];
							for ( int init = 0; init < threshold; init++ ) {
								p_sample[ init ] = -1;
							}
							p_sample[ 0 ] = i;
							matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
						} else {
							int it = -1;
							for ( int init = 0; init < threshold; init++ ) {
								if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
									it = init;
									break;
								}
							}
							if ( it > 0 && it < threshold ) {
								matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
							} else {
								it = rand() % threshold;
								matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
							}
						}
					} else {
						for ( int wr = 0; wr < vec_t.size(); wr++ ) {
							mark[ i ][ vec_t[ wr ] ] = 0;
						}
					}
					vec_t.clear();
				}
				continue;
			}
			if ( mark[ i ][ k ] == 1 ) {
				vec_t.push_back( k );
			} else if ( first && mark[ i ][ k ] == 0 ) {
				continue;
			} else {
				first = true;
				if ( vec_t.size() >= hold ) {
					if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
						int16_t * p_sample = new int16_t [ threshold ];
						for ( int init = 0; init < threshold; init++ ) {
							p_sample[ init ] = -1;
						}
						p_sample[ 0 ] = i;
						matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
					} else {
						int it = -1;
						for ( int init = 0; init < threshold; init++ ) {
							if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
								it = init;
								break;
							}
						}
						if ( it > 0 && it < threshold ) {
							matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
						} else {
							it = rand() % threshold;
							matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
						}
					}
				} else {
					for ( int wr = 0; wr < vec_t.size(); wr++ ) {
						mark[ i ][ vec_t[ wr ] ] = 0;
					}
				}
				vec_t.clear();
			}
			
			if ( k == snp_number - 1 ) {	
				if ( vec_t.size() >= hold ) {
					if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
						int16_t * p_sample = new int16_t [ threshold ];
						for ( int init = 0; init < threshold; init++ ) {
							p_sample[ init ] = -1;
						}
						p_sample[ 0 ] = i;
						matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
					} else {
						int it = -1;
						for ( int init = 0; init < threshold; init++ ) {
							if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
								it = init;
								break;
							}
						}
						if ( it > 0 && it < threshold ) {
							matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
						} else {
							it = rand() % threshold;
							matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
						}
					}
				} else {
					for ( int wr = 0; wr < vec_t.size(); wr++ ) {
						mark[ i ][ vec_t[ wr ] ] = 0;
					}
				}
				vec_t.clear();
			}
		}
	}
#ifdef DEBUG
	for ( int i = 0; i < individuals; i++ ) {
		cerr << "B" << i << " show:\t|";
		for ( int j = 0; j < snp_number; j++ ) {
			if ( mark[ i ][ j ] == 1 ) {
				cerr << ".";
			} else {
				cerr << " ";
			}
		}
		cerr << "|\n";
	}
#endif
	// generate haplotypes
	for ( int i = 0; i < individuals; i++ ) {
		int pesudo_end = -1;
		int start = -1;
		int end   = -1;
		if ( ! start_from_zero ) {
			pesudo_end = marker_start + gap - 1;
		}
		for ( int k = 0; k < snp_number; k++ ) {
			for ( int g = snp_number - 1; g > k; g-- ) {
				if ( matrix[ k ][ g ] != NULL ) {
					int indiv = 0;
					bool just_current_individual = false;
					for ( int init = 0; init < threshold; init++ ) {
						if ( matrix[ k ][ g ][ init ] != -1 ) {
							indiv++;
							if ( matrix[ k ][ g ][ init ] == i ) {
								just_current_individual = true;
							}
						} else {
							break;
						}
					}
					if ( indiv == 0 ) {
						continue;
					}
					if ( indiv == 1 && just_current_individual ) {
						continue;
					}
					vector< node_t > nodes;
					for ( int init = 0; init < indiv; init++ ) {
						int id = matrix[ k ][ g ][ init ];
						if ( id == i ) {
							continue;
						}
						node_t nt( g - k + 1, k, g );
						for ( int L = k; L <= g; L++ ) {
							nt.piece[ L - k ] = genotypes[ id ][ L + marker_start ];
						}
						if ( nodes.empty() ) {
							nt.count++;
							nt.indivs.push_back( id );
							nodes.push_back( nt );
						} else {
							bool exists = false;
							int log_L = -1;
							for ( int L = 0; L < nodes.size(); L++ ) {
								if ( nodes[ L ].compare( nt ) ) {
									exists = true;
									log_L = L;
									break;
								}
							}
							if ( exists ) {
								nodes[ log_L ].count++;
								nodes[ log_L ].indivs.push_back( id );
							} else {
								nt.count++;
								nt.indivs.push_back( id );
								nodes.push_back( nt );
							}
						}
					}

					int n_nodes = nodes.size();
					int sum = 0;
					for ( int init = 0; init < n_nodes; init++ ) {
						sum += nodes[ init ].count;
					}
					vector< double > probs;
					for ( int init = 0; init < n_nodes; init++ ) {
						double _t = double(nodes[ init ].count) / double(sum);
						probs.push_back( _t );
					}
					double r = drand48();
					double log_r = 0;
					int init_id = -1;
					for ( int init = 0; init < n_nodes; init++ ) {
						log_r += probs[ init ];
						if ( log_r >= r ) {
							init_id = init;
							break;
						}
					}
					int indiv_id = nodes[ init_id ].indivs[ 0 ];
					// init the haplotype.
					bool is_match = true;
					for ( int compare_i = k; compare_i <= g; compare_i++ ) {
						if ( abs( genotypes[ i ][ compare_i + marker_start ] - genotypes[ indiv_id ][ compare_i + marker_start ] ) > 1 ) {
							is_match = false;
							break;
						}
					}
					if ( is_match ) {
						if ( k > pesudo_end - marker_start ) {
							start = pesudo_end + 1;
							end   = k + marker_start;
							for ( int random_i = start; random_i < end; random_i++ ) {
								if ( genotypes[ i ][ random_i ] == 1 ) {
									haplotypes[ i ].first[ random_i ] = rand() % 2;
									haplotypes[ i ].second[ random_i ] = \
												1 - haplotypes[ i ].first[ random_i ] ;
								} else {
									haplotypes[ i ].first[ random_i ] = genotypes[ i ][ random_i ] / 2;
									haplotypes[ i ].second[ random_i ] = haplotypes[ i ].first[ random_i ];
								}
							//	cerr << genotypes[ i ][ random_i ] << haplotypes[ i ].first[ random_i ] << haplotypes[ i ].second[ random_i ] << "|" << endl;
							}
							pesudo_end = k - 1 + marker_start;
						}
						if ( g <= pesudo_end - marker_start ) {
							continue;
						}
						if ( k <= pesudo_end - marker_start && g > pesudo_end - marker_start ) {
							bool is_first = true;
							bool is_second = true;
							start = k + marker_start;
							end = pesudo_end;
							for ( int random_i = start; random_i <= end; random_i++ ) {
								if ( haplotypes[ i ].first[ random_i ] != \
									genotypes[ indiv_id ][ random_i ] / 2 ) {
									is_first = false;
								}
								if ( haplotypes[ i ].second[ random_i ] != \
									genotypes[ indiv_id ][ random_i ] / 2 ) {
									is_second = false;
								}
							}
							if ( (!is_first) && (!is_second ) ) {
								continue;
							}
							if ( is_first ) {
								start = k + marker_start;
								end   = g + marker_start;
								for ( int random_i = start; random_i <= end; random_i++ ) {
									haplotypes[ i ].first[ random_i ] = \
										genotypes[ indiv_id ][ random_i ] / 2;
									haplotypes[ i ].second[ random_i ] = \
										genotypes[ i ][ random_i ] - haplotypes[ i ].first[ random_i ];
								}
								pesudo_end = g + marker_start;
							}
							if ( is_second ) {
								start = k + marker_start;
								end   = g + marker_start;
								for ( int random_i = start; random_i <= end; random_i++ ) {
									haplotypes[ i ].second[ random_i ] = \
												genotypes[ indiv_id ][ random_i ] / 2;
									haplotypes[ i ].first[ random_i ] = \
										genotypes[ i ][ random_i ] - haplotypes[ i ].second[ random_i ];
								}
								pesudo_end = g + marker_start;
							}
						} else {
							start = k + marker_start;
							end   = g + marker_start;
							for ( int random_i = start; random_i <= end; random_i++ ) {
								haplotypes[ i ].first[ random_i ] = \
											genotypes[ indiv_id ][ random_i ] / 2;
								haplotypes[ i ].second[ random_i ] = \
										genotypes[ i ][ random_i ] - haplotypes[ i ].first[ random_i ];
							}
							pesudo_end = g + marker_start;
						}
					}
				}
			}
		}
	
		if ( pesudo_end - marker_start < snp_number - 1 ) {
			start = pesudo_end + 1;
			end   = snp_number - 1 + marker_start;
			for ( int init = start; init <= end; init++ ) {
				if ( genotypes[ i ][ init ] == 1 ) {
					haplotypes[ i ].first[ init ] = rand() % 2;
					haplotypes[ i ].second[ init ] = 1 - haplotypes[ i ].first[ init ];
				} else {
					haplotypes[ i ].first[ init ] = genotypes[ i ][ init ] / 2;
					haplotypes[ i ].second[ init ] = haplotypes[ i ].first[ init ];
				}
			}
		}
	}

	for ( int i = 0; i < individuals; i++ ) {
		delete [] mark[ i ];
	}
	delete [] mark;

	for ( int i = 0; i < snp_number; i++ ) {
                for ( int j = 0; j < snp_number; j++ ) {
                        delete [] matrix[ i ][ j ];
                }
                delete [] matrix[ i ];
        }
        delete [] matrix;
	matrix = NULL;
}

void _Hap_Init::split_works( int _marker_start, int _gap, int _start_from_zero, int16_t *** _matrix, vector< vector< int > > * _p_genotypes, \
			     hap_t * _haplotypes, _Real_Uniforms & RU )
{
	vector< pair< int, int > > pieces;
	int n_threads = RU.n_threads;
	int n_piece = individuals / n_threads;
	int remain  = individuals % n_threads;
	int locate = 0;
	for ( int i = 0; i < n_threads; i++ ) {
		if ( remain > 0 ) {
			pieces.push_back( pair< int, int >( locate, locate + n_piece + 1 ) );
			locate += n_piece + 1;
			remain--;
		} else {
			pieces.push_back( pair< int, int >( locate, locate + n_piece ) );
			locate += n_piece;
		}
	}
	if ( locate != individuals ) {
		cerr << "[_Hap_Init::split_works] fail to split works.\n";
		abort();
	}
	if ( n_threads != pieces.size() ) {
		cerr << "[_Hap_Init::split_works] fail to split works.\n";
		abort();
	}
	
	args = new _Init_Arg[ n_threads ];

	for ( int i = 0; i < n_threads; i++ ) {
		int _begin = pieces[ i ].first;
		int _end = pieces[ i ].second;
		
		args[ i ].init_args( _begin, _end, snp_number, individuals, _marker_start, _gap, _start_from_zero );
		args[ i ].init_pointers( _matrix, _p_genotypes, _haplotypes, &(RU.uniforms[ i ]) );
	}
}

void _Hap_Init::free_args()
{
	delete [] args;
	args = NULL;
}

_Hap_Init::~_Hap_Init()
{
	;
}

void _Hap_Init::hap_int2uint8_t( Haplotype * uint8_haps, hap_t * int_haps )
{
	for ( int i = 0; i < individuals; i++ ) {
		for ( int j = 0; j < snp_number; j++ ) {
			uint8_haps[ i * 2 + 0 ].haplotype[ j ] = int_haps[ i ].first[ j ];
			uint8_haps[ i * 2 + 1 ].haplotype[ j ] = int_haps[ i ].second[ j ];
		}
	}
}

void _Hap_Init::haplotype_init_rand(vector<vector<int> > & genotypes, hap_t * haplotypes, variate_generator< kreutzer1986 &, uniform_real<> > & uniform)
{
	snp_number = genotypes[ 0 ].size();
	individuals = genotypes.size();
	for (int i = 0; i != individuals; ++i) {
		for (int j = 0; j != snp_number; ++j) {
			if (genotypes[i][j] != 1) {
				haplotypes[i].first[j] = genotypes[i][j] / 2;
				haplotypes[i].second[j] = haplotypes[i].first[j];
			} else {
				double rand_num = uniform();
				if (rand_num < 0.5) {
					haplotypes[i].first[j] = 0;
					haplotypes[i].second[j] = 1;
				} else {
					haplotypes[i].first[j] = 1;
					haplotypes[i].second[j] = 0;
				}
			}
		}
	}
}

void * hap_init_function( void * arg )
{
	_Init_Arg * ia = ( _Init_Arg * )arg;
	int indiv_begin = ia->begin;
	int indiv_end   = ia->end;
	int snp_number  = ia->snp_number;
	int individuals = ia->individuals;
	int marker_start = ia->marker_start;
	int gap = ia->gap;
	bool start_from_zero = ia->start_from_zero;
	int16_t *** matrix = ia->matrix;
	variate_generator< kreutzer1986, uniform_real<> > * p_uniform = ia->p_uniform;
	hap_t * haplotypes = ia->haplotypes;
	vector< vector< int > > * p_genotypes = ia->p_genotypes;

	for ( int i = indiv_begin; i < indiv_end; i++ ) {
		int pesudo_end = -1;
		int start = -1;
		int end   = -1;
		if ( ! start_from_zero ) {
			pesudo_end = marker_start + gap - 1;
		}
		for ( int k = 0; k < snp_number; k++ ) {
			for ( int g = snp_number - 1; g > k; g-- ) {
				if ( matrix[ k ][ g ] != NULL ) {
					int indiv = 0;
					bool just_current_individual = false;
					for ( int init = 0; init < individuals; init++ ) {
						if ( matrix[ k ][ g ][ init ] != -1 ) {
							indiv++;
							if ( matrix[ k ][ g ][ init ] == i ) {
								just_current_individual = true;
							}
						} else {
							break;
						}
					}
					if ( indiv == 0 ) {
						continue;
					}
					if ( indiv == 1 && just_current_individual ) {
						continue;
					}
					vector< node_t > nodes;
					for ( int init = 0; init < indiv; init++ ) {
						int id = matrix[ k ][ g ][ init ];
						if ( id == i ) {
							continue;
						}
						node_t nt( g - k + 1, k, g );
						for ( int L = k; L <= g; L++ ) {
							nt.piece[ L - k ] = (*p_genotypes)[ id ][ L + marker_start ];
						}
						if ( nodes.empty() ) {
							nt.count++;
							nt.indivs.push_back( id );
							nodes.push_back( nt );
						} else {
							bool exists = false;
							int log_L = -1;
							for ( int L = 0; L < nodes.size(); L++ ) {
								if ( nodes[ L ].compare( nt ) ) {
									exists = true;
									log_L = L;
									break;
								}
							}
							if ( exists ) {
								nodes[ log_L ].count++;
								nodes[ log_L ].indivs.push_back( id );
							} else {
								nt.count++;
								nt.indivs.push_back( id );
								nodes.push_back( nt );
							}
						}
					}

					int n_nodes = nodes.size();
					int sum = 0;
					for ( int init = 0; init < n_nodes; init++ ) {
						sum += nodes[ init ].count;
					}
					vector< double > probs;
					for ( int init = 0; init < n_nodes; init++ ) {
						double _t = double(nodes[ init ].count) / double(sum);
						probs.push_back( _t );
					}
					double r = (*p_uniform)();
					double log_r = 0;
					int init_id = -1;
					for ( int init = 0; init < n_nodes; init++ ) {
						log_r += probs[ init ];
						if ( log_r >= r ) {
							init_id = init;
							break;
						}
					}
					int indiv_id = nodes[ init_id ].indivs[ 0 ];
					// init the haplotype.
					bool is_match = true;
					for ( int compare_i = k; compare_i <= g; compare_i++ ) {
						if ( abs( (*p_genotypes)[ i ][ compare_i + marker_start ] - (*p_genotypes)[ indiv_id ][ compare_i + marker_start ] ) > 1 ) {
							is_match = false;
							break;
						}
					}
					if ( is_match ) {
						if ( k > pesudo_end - marker_start ) {
					//		if ( ! start_from_zero ) {
					//			cerr << "k = " << k << "\tpesudo_end = " << pesudo_end << "\tmarker_start = " << \
					//				marker_start << "\tpesudo_end - marker_start = " << pesudo_end - marker_start << endl;
					//		}
							start = pesudo_end + 1;
							end   = k + marker_start;
							for ( int random_i = start; random_i < end; random_i++ ) {
								if ( (*p_genotypes)[ i ][ random_i ] == 1 ) {
									haplotypes[ i ].first[ random_i ] = ((*p_uniform)() > 0.5) ? 0 : 1;
									haplotypes[ i ].second[ random_i ] = \
												!haplotypes[ i ].first[ random_i ] ;
								} else {
									haplotypes[ i ].first[ random_i ] = (*p_genotypes)[ indiv_id ][ random_i ] / 2;
									haplotypes[ i ].second[ random_i ] = haplotypes[ i ].first[ random_i ];
								}
							}
							pesudo_end = k - 1 + marker_start;
						}
						if ( g <= pesudo_end - marker_start ) {
							continue;
						}
						if ( k <= pesudo_end - marker_start && g > pesudo_end - marker_start ) {
							bool is_first = true;
							bool is_second = true;
							start = k + marker_start;
							end = pesudo_end;
							for ( int random_i = start; random_i <= end; random_i++ ) {
								if ( haplotypes[ i ].first[ random_i ] != \
									(*p_genotypes)[ indiv_id ][ random_i ] / 2 ) {
									is_first = false;
								}
								if ( haplotypes[ i ].second[ random_i ] != \
									(*p_genotypes)[ indiv_id ][ random_i ] / 2 ) {
									is_second = false;
								}
							}
							if ( (!is_first) && (!is_second ) ) {
								continue;
							}
							if ( is_first ) {
								start = k + marker_start;
								end   = g + marker_start;
								for ( int random_i = start; random_i <= end; random_i++ ) {
									haplotypes[ i ].first[ random_i ] = \
										(*p_genotypes)[ indiv_id ][ random_i ] / 2;
									haplotypes[ i ].second[ random_i ] = \
										(*p_genotypes)[ i ][ random_i ] - haplotypes[ i ].first[ random_i ];
								}
								pesudo_end = g + marker_start;
							}
							if ( is_second ) {
								start = k + marker_start;
								end   = g + marker_start;
								for ( int random_i = start; random_i <= end; random_i++ ) {
									haplotypes[ i ].second[ random_i ] = \
												(*p_genotypes)[ indiv_id ][ random_i ] / 2;
									haplotypes[ i ].first[ random_i ] = \
										(*p_genotypes)[ i ][ random_i ] - haplotypes[ i ].second[ random_i ];
								}
								pesudo_end = g + marker_start;
							}
						} else {
							start = k + marker_start;
							end   = g + marker_start;
							for ( int random_i = start; random_i <= end; random_i++ ) {
								haplotypes[ i ].first[ random_i ] = \
											(*p_genotypes)[ indiv_id ][ random_i ] / 2;
								haplotypes[ i ].second[ random_i ] = \
										(*p_genotypes)[ i ][ random_i ] - haplotypes[ i ].first[ random_i ];
							}
							pesudo_end = g + marker_start;
						}
					}
				}
			}
		}
	
		if ( pesudo_end - marker_start < snp_number - 1 ) {
			start = pesudo_end + 1;
			end   = snp_number - 1 + marker_start;
			for ( int init = start; init <= end; init++ ) {
				if ( (*p_genotypes)[ i ][ init ] == 1 ) {
					haplotypes[ i ].first[ init ] = ((*p_uniform)() > 0.5) ? 0 : 1;
					haplotypes[ i ].second[ init ] = !haplotypes[ i ].first[ init ];
				} else {
					haplotypes[ i ].first[ init ] = (*p_genotypes)[ i ][ init ] / 2;
					haplotypes[ i ].second[ init ] = haplotypes[ i ].first[ init ];
				}
			}
		}
	}
}

void _Hap_Init::haplotype_init( vector< vector< int > > & genotypes, hap_t * haplotypes, _Real_Uniforms & RU )
{
	snp_number = genotypes[ 0 ].size();
	individuals = genotypes.size();
	int gap = 100;
	int distance = 10000;
	vector< pair< int, int > > pieces;
	int marker_start = 0;
	int marker_end   = 0;
	int snp_index = 0;
	while ( snp_index < snp_number ) {
		if ( snp_index == 0 ) {
			marker_start = snp_index;
			snp_index += distance;
			if ( snp_index >= snp_number ) {
				marker_end = snp_number - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
				break;
			} else {
				marker_end = snp_index - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
			}
		} else {
			snp_index -= gap;
			marker_start = snp_index;
			snp_index += distance;
			if ( snp_index >= snp_number ) {
				marker_end = snp_number - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
				break;
			} else {
				marker_end = snp_index - 1;
				pieces.push_back( pair< int, int >( marker_start, marker_end ) );
			}
		}
	}
	
	bool start_from_zero = true;
	
	for ( int i = 0; i < pieces.size(); i++ ) {
		marker_start = pieces[ i ].first;
		marker_end   = pieces[ i ].second;
		haplotype_init_piece( genotypes, haplotypes, marker_start, marker_end, gap, start_from_zero, RU );
		start_from_zero = false;
	}
}

void _Hap_Init::haplotype_init_piece( vector< vector< int > > & genotypes, hap_t * haplotypes, int marker_start, int marker_end, \
						int gap, bool start_from_zero, _Real_Uniforms & RU )
{
	snp_number = marker_end - marker_start + 1;
	individuals = genotypes.size();
	// alloc a matrix for store hom-genotype segments.
	int16_t *** matrix = new int16_t ** [ snp_number ];
	for ( int i = 0; i < snp_number; i++ ) {
		matrix[ i ] = new int16_t * [ snp_number ];
		for ( int j = 0; j < snp_number; j++ ) {
			matrix[ i ][ j ] = NULL;
		}
	}

	uint8_t ** mark = new uint8_t * [ individuals ];
	for ( int i = 0; i < individuals; i++ ) {
		mark[ i ] = new uint8_t [ snp_number ];
		memset( mark[ i ], 0, snp_number );
	}

	{
		vector< vector < int > > groups;
		for ( int i = 0; i < individuals; i++ ) {
			vector< int > group_t;
			for ( int j = 0; j < snp_number; j++ ) {
				if ( genotypes[ i ][ j + marker_start ] != 1 ) {
					group_t.push_back( j );
				}
			}
			groups.push_back( group_t );
		}		

		for ( int i = 0; i < individuals; i++ ) {
			int n = groups[ i ].size();
			for ( int j = 0; j < n; j++ ) {
				mark[ i ][ groups[ i ][ j ] ] = 1;
			}
		}
	}
#ifdef DEBUG
	for ( int i = 0; i < individuals; i++ ) {
		cerr << "B" << i << " show:\t|";
		for ( int j = 0; j < snp_number; j++ ) {
			if ( mark[ i ][ j ] == 1 ) {
				cerr << ".";
			} else {
				cerr << " ";
			}
		}
		cerr << "|\n";
	}
#endif
	int hold = 5;

	for ( int i = 0; i < individuals; i++ ) {
		vector< int > vec_t;
		bool first = true;
		for ( int k = 0; k < snp_number; k++ ) {
			if ( first && mark[ i ][ k ] == 1 ) {
				vec_t.push_back( k );
				first = false;
				if ( k == snp_number - 1 ) {
					if ( vec_t.size() >= hold ) {
						if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
							int16_t * p_sample = new int16_t [ individuals ];
							for ( int init = 0; init < individuals; init++ ) {
								p_sample[ init ] = -1;
							}
							p_sample[ 0 ] = i;
							matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
						} else {
							int it = -1;
							for ( int init = 0; init < individuals; init++ ) {
								if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
									it = init;
									break;
								}
							}
							matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
						}
					} else {
						for ( int wr = 0; wr < vec_t.size(); wr++ ) {
							mark[ i ][ vec_t[ wr ] ] = 0;
						}
					}
					vec_t.clear();
				}
				continue;
			}
			if ( mark[ i ][ k ] == 1 ) {
				vec_t.push_back( k );
			} else if ( first && mark[ i ][ k ] == 0 ) {
				continue;
			} else {
				first = true;
				if ( vec_t.size() >= hold ) {
					if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
						int16_t * p_sample = new int16_t [ individuals ];
						for ( int init = 0; init < individuals; init++ ) {
							p_sample[ init ] = -1;
						}
						p_sample[ 0 ] = i;
						matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
					} else {
						int it = -1;
						for ( int init = 0; init < individuals; init++ ) {
							if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
								it = init;
								break;
							}
						}
						matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
					}
				} else {
					for ( int wr = 0; wr < vec_t.size(); wr++ ) {
						mark[ i ][ vec_t[ wr ] ] = 0;
					}
				}
				vec_t.clear();
			}
			
			if ( k == snp_number - 1 ) {	
				if ( vec_t.size() >= hold ) {
					if ( (matrix[ vec_t[ 0 ] ][ vec_t.back() ]) == NULL ) {
						int16_t * p_sample = new int16_t [ individuals ];
						for ( int init = 0; init < individuals; init++ ) {
							p_sample[ init ] = -1;
						}
						p_sample[ 0 ] = i;
						matrix[ vec_t[ 0 ] ][ vec_t.back() ] = p_sample;
					} else {
						int it = -1;
						for ( int init = 0; init < individuals; init++ ) {
							if ( matrix[ vec_t[ 0 ] ][ vec_t.back() ][ init ] == -1 ) {
								it = init;
								break;
							}
						}
						matrix[ vec_t[ 0 ] ][ vec_t.back() ][ it ] = i;
					}
				} else {
					for ( int wr = 0; wr < vec_t.size(); wr++ ) {
						mark[ i ][ vec_t[ wr ] ] = 0;
					}
				}
				vec_t.clear();
			}
		}
	}
#ifdef DEBUG
	for ( int i = 0; i < individuals; i++ ) {
		cerr << "B" << i << " show:\t|";
		for ( int j = 0; j < snp_number; j++ ) {
			if ( mark[ i ][ j ] == 1 ) {
				cerr << ".";
			} else {
				cerr << " ";
			}
		}
		cerr << "|\n";
	}
#endif

	// generate haplotypes
	split_works( marker_start, gap, start_from_zero, matrix, &genotypes, haplotypes, RU );
	int n_threads = RU.n_threads;
	pthread_t * thread_ids = ( pthread_t * )malloc( sizeof( pthread_t ) * n_threads );
	for ( int i = 0; i < n_threads; i++ ) {
		if ( 0 != pthread_create( &(thread_ids[ i ]), NULL, hap_init_function, ( void * )(&(args[ i ])) ) ){
			cerr << "[_Hap_Init::haplotype_init_piece] fail to create " << i << " thread\n";
			abort();
		}
	}
	for ( int i = 0; i < n_threads; i++ ) {
		if ( 0 != pthread_join( thread_ids[ i ], NULL ) ) {
                        cerr << "[_Hap_Init::haplotype_init_piece] fail to wait thread " << i << " exit.\n";
                        abort();
                }
	}
	free( thread_ids );
	free_args();

	for ( int i = 0; i < individuals; i++ ) {
		delete [] mark[ i ];
	}
	delete [] mark;

	for ( int i = 0; i < snp_number; i++ ) {
                for ( int j = 0; j < snp_number; j++ ) {
                        delete [] matrix[ i ][ j ];
                }
                delete [] matrix[ i ];
        }
        delete [] matrix;
	matrix = NULL;
}
