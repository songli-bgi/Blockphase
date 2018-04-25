#include "Haplotype_graph.h"

using namespace std;

typedef map< ele , uint32_t , Comp > Map;
typedef map< ele , uint32_t , Comp >::iterator Map_Iterator;

void Haplo_graph::split_blocks( Haplotype * haplotypes , int state_k, Haplotype * reference, int _n_reference, \
				variate_generator< kreutzer1986, uniform_real<> > & uniform, bool reference_only )
{
	_nodes.clear();
	int core_snp = int(uniform() * snp_number);
	while( (core_snp <= 20) || (core_snp >= (snp_number - 20)) ){
		core_snp = int(uniform() * snp_number);
	}

	int max_hap_size = 500;

	int count_right = core_snp;
	bool right_flag = false;
	while( count_right < snp_number ){
		int old_count = count_right;
		int local_hap_size = 0;
		Map _t_map_;
		while( _t_map_.size() < static_cast< size_t >( state_k ) && local_hap_size <= max_hap_size ){
			_t_map_.clear();
			if ( count_right >= snp_number ){
				right_flag = true;
				_nodes.push_back( count_right - old_count );
				break;
			}
			count_right++;
			local_hap_size++;
			if ( ! reference_only ) {
				int length = count_right - old_count;
				ele hap_seg( length );
				for ( int sample_index = 0; sample_index < individuals; sample_index++ ){
					memcpy( hap_seg.p , haplotypes[ sample_index * 2 + 0 ].haplotype + old_count , length );
					_t_map_[ hap_seg ]++;
					memcpy( hap_seg.p , haplotypes[ sample_index * 2 + 1 ].haplotype + old_count , length );
					_t_map_[ hap_seg ]++;
				}
			}
			if ( NULL != reference ) {
				int length = count_right - old_count;
				ele hap_seg( length );
				for ( int hap_index = 0; hap_index < _n_reference; hap_index++ ) {
					memcpy( hap_seg.p, reference[ hap_index ].haplotype + old_count, length );
					_t_map_[ hap_seg ]++;
				}
			}
		}
		if ( count_right == snp_number && right_flag ){
			continue;
		}
		_nodes.push_back( count_right - old_count - 1 );
		count_right--;
	}

	// extern to left.
	int count_left = core_snp;
	bool left_flag = false;
	while( count_left > 0 ){
		int old_count_left = count_left;
		int local_hap_size = 0;
		Map _t_map_;
		while( _t_map_.size() < static_cast< size_t >( state_k ) && local_hap_size <= max_hap_size ){
			_t_map_.clear();
			if ( count_left <= 0 ){
				left_flag = true;
				_nodes.push_front( old_count_left - count_left );
				break;
			}
			count_left--;
			local_hap_size++;
			if ( ! reference_only ) {
				int length = old_count_left - count_left;
				ele hap_seg( length );
				for ( int sample_index = 0; sample_index < individuals; sample_index++ ){
					memcpy( hap_seg.p , haplotypes[ sample_index * 2 + 0 ].haplotype + count_left , length );
					_t_map_[ hap_seg ]++;
					memcpy( hap_seg.p , haplotypes[ sample_index * 2 + 1 ].haplotype + count_left , length );
					_t_map_[ hap_seg ]++;
				}
			}
			
			if ( NULL != reference ) {
				int length = old_count_left - count_left;
				ele hap_seg( length );
				for ( int hap_index = 0; hap_index < _n_reference; hap_index++ ) {
					memcpy( hap_seg.p, reference[ hap_index ].haplotype + count_left, length );
					_t_map_[ hap_seg ]++;
				}
			}
		}
		if ( count_left == 0 && left_flag ){
			continue;
		}
		_nodes.push_front( old_count_left - count_left - 1 );
		count_left++;
	}
	n_kjds = _nodes.size();
	_kjds = new _KJD[ n_kjds ];
}

void Haplo_graph::split_blocks(vector< int > & blocks, int size, variate_generator< kreutzer1986, uniform_real<> > & uniform )
{
	blocks.clear();
	int locate = 0;
	while (locate < snp_number) {
		if (blocks.empty()) {
			double p = uniform();
			int s = static_cast<int>(p * size + size / 2. + 0.499);
			blocks.push_back(s);
			locate += s;
		} else {
			blocks.push_back(size);
			locate += size;
		}
	}
	blocks.back() -= (locate - snp_number);
}

void Haplo_graph::list2vector( vector< int > & blocks )
{
	int locate = 0;
	list< int >::iterator iter1 = _nodes.begin();
	list< int >::iterator iter2 = _nodes.end();
	for ( ; iter1 != iter2; iter1++ ) {
		blocks.push_back( *iter1 );
		locate += *iter1;
	}
	
	if ( locate != snp_number ) {
	#ifdef DEBUG
		cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
	#endif
		cerr << "[Haplo_graph::list2vector] fail to push the information of list into a vector object.\n";
		abort();
	}
}

void Haplo_graph::construct_hg( Haplotype * haplotypes, int individuals_ID, vector< int > & blocks, Haplotype * reference, int _n_reference, \
				bool reference_only )
{
	int n_blocks = blocks.size();
	if ( n_blocks != n_kjds ) {
		delete [] _kjds;
		_kjds = new _KJD[ n_blocks ];
		n_kjds = n_blocks;
	}
	if ( ! reference_only ) {
		for ( int sample_index = 0; sample_index < individuals; sample_index++ ) {
			if ( sample_index == individuals_ID ) {
				continue;
			}
			for ( int phase = 0; phase < 2; phase++ ) {
				bool flag = true;
				int locate( 0 ), prev_id( 0 ), next_id( 0 );
				for ( int index = 0; index < n_blocks; index++ ) {
					int hap_len = blocks[ index ];
					ele _hap( hap_len );
					_kjds[ index ].length = hap_len;
					prev_id = next_id;
					memcpy( _hap.p, haplotypes[ sample_index * 2 + phase ].haplotype + locate, hap_len );
					if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
						_kjds[ index ].haps[ _hap ].second++;
					} else { // not exists, alloc a new id and add one.
						_kjds[ index ].add( _hap );
					}
					next_id = _kjds[ index ].getid( _hap );
					if ( flag ) {
						flag = false;
					} else {
						if ( _kjds[ index ].prev_edges.count( next_id ) ) {
							if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
								_kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
							} else {
								_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
							}
						} else {
							_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
						}
					}
					locate += hap_len;
				}
				if ( locate != snp_number ) {
				#ifdef DEBUG
					cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
				#endif
					cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
					abort();
				}
			}
		}
	}
	
	if ( reference != NULL ) {
                n_reference = _n_reference;
                for ( int sample_index = 0; sample_index < n_reference; sample_index++ ) {
                    	bool flag = true;
                        int locate( 0 ), prev_id( 0 ), next_id( 0 );
                       	for ( int index = 0; index < n_blocks; index++ ) {
				int hap_len = blocks[ index ];
                                _kjds[ index ].length = hap_len;
				ele _hap( hap_len );
                                prev_id = next_id;
                                memcpy( _hap.p, reference[ sample_index ].haplotype + locate, hap_len );
                                if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
                                     	_kjds[ index ].haps[ _hap ].second++;
                                } else { // not exists, alloc a new id and add one.
                                        _kjds[ index ].add( _hap );
                                }
                                next_id = _kjds[ index ].getid( _hap );
                                if ( flag ) {
                                        flag = false; // skip the first node.
                                } else {
                                        if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                                if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
                                                        _kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
                                                } else {
                                                        _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                                }
                                        } else {
                                                _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                        }
                               	}
                                locate += hap_len;
                     	}
                        if ( locate != snp_number ) {
                        #ifdef DEBUG
                              	cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                        #endif
                                cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
                               	abort();
                       	}
                }
        }
}

void Haplo_graph::construct_hg( Haplotype * haplotypes, map< int, int > & exclude_samples, vector< int > & blocks, Haplotype * reference, int _n_reference, \
				bool reference_only )
{
	int n_blocks = blocks.size();
	if ( n_blocks != n_kjds ) {
		delete [] _kjds;
		_kjds = new _KJD[ n_blocks ];
		n_kjds = n_blocks;
	}
	if ( ! reference_only ) {
		for ( int sample_index = 0; sample_index < individuals; sample_index++ ) {
			if ( exclude_samples.count( sample_index ) ) {
				continue;
			}
			for ( int phase = 0; phase < 2; phase++ ) {
				bool flag = true;
				int locate( 0 ), prev_id( 0 ), next_id( 0 );
				for ( int index = 0; index < n_blocks; index++ ) {
					int hap_len = blocks[ index ];
					ele _hap( hap_len );
					_kjds[ index ].length = hap_len;
					prev_id = next_id;
					memcpy( _hap.p, haplotypes[ sample_index * 2 + phase ].haplotype + locate, hap_len );
					if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
						_kjds[ index ].haps[ _hap ].second++;
					} else { // not exists, alloc a new id and add one.
						_kjds[ index ].add( _hap );
					}
					next_id = _kjds[ index ].getid( _hap );
					if ( flag ) {
						flag = false;
					} else {
						if ( _kjds[ index ].prev_edges.count( next_id ) ) {
							if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
								_kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
							} else {
								_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
							}
						} else {
							_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
						}
					}
					locate += hap_len;
				}
				if ( locate != snp_number ) {
				#ifdef DEBUG
					cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
				#endif
					cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
					abort();
				}
			}
		}
	}
	
	if ( reference != NULL ) {
                n_reference = _n_reference;
                for ( int sample_index = 0; sample_index < n_reference; sample_index++ ) {
                    	bool flag = true;
                        int locate( 0 ), prev_id( 0 ), next_id( 0 );
                       	for ( int index = 0; index < n_blocks; index++ ) {
				int hap_len = blocks[ index ];
                                _kjds[ index ].length = hap_len;
				ele _hap( hap_len );
                                prev_id = next_id;
                                memcpy( _hap.p, reference[ sample_index ].haplotype + locate, hap_len );
                                if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
                                     	_kjds[ index ].haps[ _hap ].second++;
                                } else { // not exists, alloc a new id and add one.
                                        _kjds[ index ].add( _hap );
                                }
                                next_id = _kjds[ index ].getid( _hap );
                                if ( flag ) {
                                        flag = false; // skip the first node.
                                } else {
                                        if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                                if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
                                                        _kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
                                                } else {
                                                        _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                                }
                                        } else {
                                                _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                        }
                               	}
                                locate += hap_len;
                     	}
                        if ( locate != snp_number ) {
                        #ifdef DEBUG
                              	cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                        #endif
                                cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
                               	abort();
                       	}
                }
        }
}

void Haplo_graph::modify_hg( Haplotype * haplotypes, int prev_individuals_ID, int next_individuals_ID, vector< int > & blocks )
{
        int n_blocks = blocks.size();
        if ( n_blocks != n_kjds ) {
                delete [] _kjds;
                _kjds = new _KJD[ n_blocks ];
                n_kjds = n_blocks;
        }

        for ( int phase = 0; phase < 2; phase++ ) {
                bool flag = true;
                int locate( 0 ), prev_id( 0 ), next_id( 0 );
                for ( int index = 0; index < n_blocks; index++ ) {
                        int hap_len = blocks[ index ];
                        _kjds[ index ].length = hap_len;
                        ele _hap( hap_len );
                        prev_id = next_id;
                        memcpy( _hap.p, haplotypes[ prev_individuals_ID * 2 + phase ].haplotype + locate, hap_len );
                        if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
                                _kjds[ index ].haps[ _hap ].second++;
                        } else { // not exists, alloc a new id and add one.
                                _kjds[ index ].add( _hap );
                        }
                        next_id = _kjds[ index ].getid( _hap );
                        if ( flag ) {
                                flag = false;
                        } else {
                                if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                        if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
                                                _kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
                                        } else {
                                                _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                        }
                                } else {
                                        _kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                }
                        }
                        locate += hap_len;
                }
                if ( locate != snp_number ) {
                #ifdef DEBUG
                        cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                #endif
                        cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
			abort();
                 }
        }

        for ( int phase = 0; phase < 2; phase++ ) {
                bool flag = true;
                int locate( 0 ), prev_id( 0 ), next_id( 0 );
                for ( int index = 0; index < n_blocks; index++ ) {
                        int hap_len = blocks[ index ];
                        _kjds[ index ].length = hap_len;
                        ele _hap( hap_len );
                        prev_id = next_id;
                        memcpy( _hap.p, haplotypes[ next_individuals_ID * 2 + phase ].haplotype + locate, hap_len );
                        next_id = _kjds[ index ].getid( _hap );
                        if ( _kjds[ index ].haps.count( _hap ) ) { // exists, sub one or erase it.
                                if ( _kjds[ index ].haps[ _hap ].second <= 1 ) {
                                        _KJD_TYPE_ITERATOR iter = _kjds[ index ].haps.find( _hap );
                                        if ( iter != _kjds[ index ].haps.end() )
                                                _kjds[ index ].haps.erase( iter );
                                        else
                                                abort();
                                } else {
                                        _kjds[ index ].haps[ _hap ].second--;
                                }
                        } else { // not exists, abort()
                                cerr << "[Haplo_graph::modify_hg] haplotype segment must be found in haplotype graph. \
					 When the word appear, it indicates un-normal.\n";
                                abort();
                        }
                        if ( flag ) {
                                flag = false;
                        } else {
                                if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                        if ( _kjds[ index ].prev_edges[ next_id ][ prev_id ] <= 1 ) {
                                                map< int, int >::iterator iter = _kjds[ index ].prev_edges[ next_id ].find( prev_id );
                                                if ( iter != _kjds[ index ].prev_edges[ next_id ].end() )
                                                        _kjds[ index ].prev_edges[ next_id ].erase( iter );
                                                else
                                                        abort();
                                                if ( _kjds[ index ].prev_edges[ next_id ].empty() ) {
                                                        map< int, map< int, int > >::iterator iter_t = _kjds[ index ].prev_edges.find( next_id );
                                                        _kjds[ index ].prev_edges.erase( iter_t );
                                                }
                                        } else {
						_kjds[ index ].prev_edges[ next_id ][ prev_id ]--;
                                        }
                                } else {
                                        cerr << "[Haplo_graph::modify_hg] enter a wrong code region.\n";
                                        abort();
                                }
                        }
                        locate += hap_len;
                }
                if ( locate != snp_number ) {
                #ifdef DEBUG
                        cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                #endif
                        cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
                        abort();
                }
        }

        modify_KJDS();
}

void Haplo_graph::modify_hg( Haplotype * haplotypes, vector< int > & prev_samples, vector< int > & next_samples, vector< int > & blocks )
{
        int n_blocks = blocks.size();
        if ( n_blocks != n_kjds ) {
                delete [] _kjds;
                _kjds = new _KJD[ n_blocks ];
                n_kjds = n_blocks;
        }

	for ( size_t spl_index = 0; spl_index < prev_samples.size(); spl_index++ ) {
		int current_individual_ID = prev_samples[ spl_index ];
        	for ( int phase = 0; phase < 2; phase++ ) {
                	bool flag = true;
                	int locate( 0 ), prev_id( 0 ), next_id( 0 );
                	for ( int index = 0; index < n_blocks; index++ ) {
                        	int hap_len = blocks[ index ];
                        	_kjds[ index ].length = hap_len;
                        	ele _hap( hap_len );
                        	prev_id = next_id;
                        	memcpy( _hap.p, haplotypes[ current_individual_ID * 2 + phase ].haplotype + locate, hap_len );
                        	if ( _kjds[ index ].haps.count( _hap ) ) { // exists, just add one.
                                	_kjds[ index ].haps[ _hap ].second++;
                        	} else { // not exists, alloc a new id and add one.
                                	_kjds[ index ].add( _hap );
                        	}
                        	next_id = _kjds[ index ].getid( _hap );
                        	if ( flag ) {
                                	flag = false;
                        	} else {
                                	if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                        	if ( _kjds[ index ].prev_edges[ next_id ].count( prev_id ) ) {
                                                	_kjds[ index ].prev_edges[ next_id ][ prev_id ]++;
                                        	} else {
                                                	_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                        	}
                                	} else {
                                        	_kjds[ index ].prev_edges[ next_id ][ prev_id ] = 1;
                                	}
                        	}
                        	locate += hap_len;
                	}
                	if ( locate != snp_number ) {
                	#ifdef DEBUG
                        	cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                	#endif
                        	cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
				abort();
                 	}
        	}
	}

	for ( size_t spl_index = 0; spl_index < next_samples.size(); spl_index++ ) {
		int current_individual_ID = next_samples[ spl_index ];
        	for ( int phase = 0; phase < 2; phase++ ) {
                	bool flag = true;
                	int locate( 0 ), prev_id( 0 ), next_id( 0 );
                	for ( int index = 0; index < n_blocks; index++ ) {
                        	int hap_len = blocks[ index ];
                        	_kjds[ index ].length = hap_len;
                        	ele _hap( hap_len );
                        	prev_id = next_id;
                        	memcpy( _hap.p, haplotypes[ current_individual_ID * 2 + phase ].haplotype + locate, hap_len );
                        	next_id = _kjds[ index ].getid( _hap );
                        	if ( _kjds[ index ].haps.count( _hap ) ) { // exists, sub one or erase it.
                                	if ( _kjds[ index ].haps[ _hap ].second <= 1 ) {
                                        	_KJD_TYPE_ITERATOR iter = _kjds[ index ].haps.find( _hap );
                                        	if ( iter != _kjds[ index ].haps.end() )
                                                	_kjds[ index ].haps.erase( iter );
                                        	else
                                                	abort();
                                	} else {
                                        	_kjds[ index ].haps[ _hap ].second--;
                                	}
                        	} else { // not exists, abort()
                                	cerr << "[Haplo_graph::modify_hg] haplotype segment must be found in haplotype graph. \
						 When the word appear, it indicates un-normal.\n";
                                	abort();
                        	}
                        	if ( flag ) {
                                	flag = false;
                        	} else {
                                	if ( _kjds[ index ].prev_edges.count( next_id ) ) {
                                        	if ( _kjds[ index ].prev_edges[ next_id ][ prev_id ] <= 1 ) {
                                                	map< int, int >::iterator iter = _kjds[ index ].prev_edges[ next_id ].find( prev_id );
                                                	if ( iter != _kjds[ index ].prev_edges[ next_id ].end() )
                                                        	_kjds[ index ].prev_edges[ next_id ].erase( iter );
                                                	else
                                                        	abort();
                                                	if ( _kjds[ index ].prev_edges[ next_id ].empty() ) {
                                                        	map< int, map< int, int > >::iterator iter_t = _kjds[ index ].prev_edges.find( next_id );
                                                        	_kjds[ index ].prev_edges.erase( iter_t );
                                                	}
                                        	} else {
							_kjds[ index ].prev_edges[ next_id ][ prev_id ]--;
                                        	}
                                	} else {
                                        	cerr << "[Haplo_graph::modify_hg] enter a wrong code region.\n";
                                        	abort();
                                	}
                        	}
                        	locate += hap_len;
                	}
                	if ( locate != snp_number ) {
                	#ifdef DEBUG
                        	cerr << "locate = " << locate << "\tsnp_number = " << snp_number << endl;
                	#endif
                        	cerr << "[Haplo_graph::construct_inter_edge] fail to construct haplotype graph.\n";
                        	abort();
                	}
        	}
	}

        modify_KJDS();
}

void Haplo_graph::modify_KJDS()
{
        map< int, int > prev_orders;
        map< int, int > next_orders;
        for ( int i = 0; i < n_kjds; i++ ) {
                vector< int > ids;
                if ( i != 0 ) {
                        prev_orders.swap( next_orders );
                        next_orders.clear();
                }
                _KJD_TYPE_ITERATOR iter1 = _kjds[ i ].haps.begin();
                _KJD_TYPE_ITERATOR iter2 = _kjds[ i ].haps.end();
                for ( ; iter1 != iter2; iter1++ ) {
                        ids.push_back( iter1->second.first );
                }

                sort( ids.begin(), ids.end() );
                for ( size_t j = 0; j < ids.size(); j++ ) {
                        next_orders[ ids[ j ] ] = j;
                }
                // have been build up a map to store ids. Then I will modify the ids.
                iter1 = _kjds[ i ].haps.begin();
                for ( ; iter1 != iter2; iter1++ ) {
                        iter1->second.first = next_orders[ iter1->second.first ];
                }
                _kjds[ i ].setid( ids.size() );
                // finish inter edges.
                if ( i == 0 ) { // skip the first node.
                        continue;
                }
                map< int, map< int, int > > prev_edges_t;
                map< int, map< int, int > >::iterator it_p1 = _kjds[ i ].prev_edges.begin();
                map< int, map< int, int > >::iterator it_p2 = _kjds[ i ].prev_edges.end();

                map< int, int >::iterator it_int1;
                map< int, int >::iterator it_int2;
		for ( ; it_p1 != it_p2; it_p1++ ) {
                        map< int, int > _t;
                        it_int1 = it_p1->second.begin();
                        it_int2 = it_p1->second.end();
                        for ( ; it_int1 != it_int2; it_int1++ ) {
                                _t[ prev_orders[ it_int1->first ] ] = it_int1->second;
                        }
                        prev_edges_t[ next_orders[ it_p1->first ] ] = _t;
                }
                 _kjds[ i ].prev_edges.swap( prev_edges_t );
        }
}

void Haplo_graph::free_kjds()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		_kjds[ i ].clear();
	}
	delete [] _kjds;
	n_kjds = 0;
}

void Haplo_graph::show_kjds()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		int len = _kjds[ i ].length;
		_KJD_TYPE_ITERATOR iter1 = _kjds[ i ].haps.begin();
		_KJD_TYPE_ITERATOR iter2 = _kjds[ i ].haps.end();
		cerr << "block " << i << ":\t";
		for ( ; iter1 != iter2; iter1++ ) {
			for ( int j = 0; j < len; j++ ) {
				const uint8_t * p = iter1->first.p;
				cerr << int( p[ j ] );
			}
			cerr << "(" << iter1->second.first << "|" << iter1->second.second << ")\t";
		}
		cerr << "\n";
	}
	
	for ( int i = 0; i < n_kjds; i++ ) {
		cerr << "block " << i << " has " << _kjds[ i ].prev_edges.size() << " haps.\n";
		int whole = 0;
		for( size_t id = 0; id < _kjds[ i ].prev_edges.size(); id++ ) {
			cerr << "\t\t hap " << id << " has " << _kjds[ i ].prev_edges[ id ].size() << " previous edges:";
			int total = 0;
			map< int, int >::iterator iter1 = _kjds[ i ].prev_edges[ id ].begin();
			map< int, int >::iterator iter2 = _kjds[ i ].prev_edges[ id ].end();
			for ( ; iter1 != iter2; iter1++ ) {
				cerr << "\t(" << iter1->first << "|" << iter1->second << ")";
				total += iter1->second;
			}
			cerr << "\ttotal haps = " << total << "\n";
			whole += total;
		}
		cerr << "**************** THE NUMBER OF ALL HAPLOTYPES IS " << whole << " ***************\n";
	}
}

void Haplo_graph::show_kjd_array()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		int s = _kjd_array[ i ]->size;
		int total = 0;
		cerr << "block " << i << " has " << s << " elements:\t";
		for ( int j = 0; j < s; j++ ) {
			cerr << "<" << j << ">";
			int l = _kjd_array[ i ]->length;
			uint8_t * q = _kjd_array[ i ]->piece[ j ].p;
			uint32_t * w = _kjd_array[ i ]->weigth;
			for ( int k = 0; k < l; k++ ) {
				cerr << int(q[ k ]);
			}
			cerr << "<" << w[ j ] << ">\t";
			total += w[ j ];
		}
		cerr << "\ttotal haplotype = " << total << "\n";
	}
}

void Haplo_graph::load_kjd_array()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		int l = _kjds[ i ].length;
		int s = _kjds[ i ].haps.size();
		Segment * p_seg = new Segment( s );
		p_seg->length = l;
		_KJD_TYPE_ITERATOR iter1 = _kjds[ i ].haps.begin();
		_KJD_TYPE_ITERATOR iter2 = _kjds[ i ].haps.end();
		for ( ; iter1 != iter2; iter1++ ) {
			int id = iter1->second.first;
			p_seg->piece[ id ].Alloc( l );
			memcpy( p_seg->piece[ id ].p, iter1->first.p, l );
			p_seg->weigth[ id ] = uint32_t(iter1->second.second);
		}
		_kjd_array.push_back( p_seg );
	}
}

void Haplo_graph::free_kjd_array()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		delete _kjd_array[ i ];
	}
	_kjd_array.clear();
}

void Haplo_graph::clear()
{
	for ( int i = 0; i < n_kjds; i++ ) {
		_kjds[ i ].clear();
	}
	free_kjd_array();
}
