#include "HMM_Algo.h"

void _HMM_Algo::transmit_state()
{
	int node_size = _hg.n_kjds;
	_start_array = new int[ node_size ];
	_end_array = new int[ node_size ];
	int locate = 0;
	for ( int i = 0; i < node_size; i++ ){
		_start_array[ i ] = locate;
		for ( int j = 0; j < _hg._kjd_array[ i ]->length; j++ ) {
			int marker = locate + j;
			if ( j == _hg._kjd_array[ i ]->length - 1 ) {
				_block_flags[ marker ] = true;
			} else {
				_block_flags[ marker ] = false;
			}
			_block_offsets[ marker ] = j;
			_block_IDs[ marker ] = i;
		}
		locate += _hg._kjd_array[ i ]->length;
		_end_array[ i ] = locate - 1;
	}

	_block_flags[ _snp_number - 1 ] = false;

	if ( locate != _snp_number ){
		cerr << "[_HMM_Algo::transmit_state] check locate and snp_number fail.\n";
        	abort();
	}
	_marginal_i = new float [ 2 * _individuals + _hg.n_reference ];
	_marginal_j = new float [ 2 * _individuals + _hg.n_reference ];
	_marginal_len += _hg.n_reference;
}

void _HMM_Algo::set_hapNumber( int n_hidden_states )
{
	_haplotype_number = float( n_hidden_states );
}

void _HMM_Algo::calculate_errors( float * mutation_rate )
{
	for ( int i = 0; i < _snp_number; i++ )
	{
		float rate = mutation_rate[ i ];
                _errors[ i ][ 0 ][ 0 ] = _errors[ i ][ 2 ][ 2 ] = (1. - rate ) * (1. - rate);
                _errors[ i ][ 0 ][ 1 ] = _errors[ i ][ 2 ][ 1 ] = 2 * (1. - rate) * rate;
                _errors[ i ][ 0 ][ 2 ] = _errors[ i ][ 2 ][ 0 ] = rate * rate;

                _errors[ i ][ 1 ][ 0 ] = _errors[ i ][ 1 ][ 2 ] =  (1. - rate) * rate;
                _errors[ i ][ 1 ][ 1 ] = (1. - rate) * (1. - rate) + rate * rate;
	}
}

void _HMM_Algo::show_leftmatrix()
{
	int node_size = _hg.n_kjds;
	int locate = 0;
	cerr << "***************************************************************LEFTMATRIX****************************************************************\n";
	for ( int i = 0; i < node_size; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		for ( int j = 0; j < n_site; j++ ) {
			int marker = locate + j;
			cerr << "\t\t\t\tmarker = " << marker << "\t\t\t\t\n";
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = 0; l < n_haps; l++ ) {
					cerr << "_LeftMatrix[ " << marker << " ][ " << k << " ][ " << l << " ] = " << \
						 _LeftMatrix[ marker ][ k * n_haps + l ] << endl;
				}
			}
		}
		locate += n_site;
	}
	cerr << "***************************************************************LEFTMATRIX****************************************************************\n";
}

void _HMM_Algo::alloc_leftmatrix()
{
	_LeftMatrix = new float * [ _snp_number ];
	int node_size = _hg.n_kjds;
	int locate = 0;
	for ( int i = 0; i < node_size; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		for ( int j = 0; j < n_site; j++ ) {
			int marker = locate + j;
			_LeftMatrix[ marker ] = new float [ n_haps * n_haps ];
		}
		locate += n_site;
	}
	if ( _snp_number != locate ) {
		cerr << "[_HMM_Algo::alloc_leftmatrix] alloc memory for HMM error.\n";
		abort();
	}
}

void _HMM_Algo::free_leftmatrix()
{
	for ( int i = 0; i < _snp_number; i++ ) {
		delete [] _LeftMatrix[ i ];
	}
	delete [] _LeftMatrix;
}

void _HMM_Algo::update_haplotype( Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
			Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, float * mutation_rate  )
{
	float * alpha = _LeftMatrix[ _snp_number - 1 ];
	double sum = 0;
	int cluster_id = _block_IDs[ _snp_number - 1 ];
	int n_haps = _hg._kjd_array[ cluster_id ]->size;
	for ( int i = 0; i < n_haps; i++ ) {
		sum += alpha[ i * n_haps + i ];
		for ( int j = i + 1; j < n_haps; j++ ) {
			sum += alpha[ i * n_haps + j ] * 2;
		}
	}
	double choice = uniform() * sum;
	int first = -1;
	int second = -1;
	double tmp = 0;
	for ( int i = 0; i < n_haps; i++ ) {
		 tmp -= alpha[ i * n_haps + i ];
		for ( int j = i; j < n_haps; j++ ) {
			tmp += alpha[ i * n_haps + j ] * 2;
			if ( tmp > choice ) {
				 if (uniform() > 0.5) {
					first = i;
					second = j;
				} else {
					first = j;
					second = i;
				}
				goto IBM;
			}
		}
	} // back samples and randomly select the haplotype pair

 IBM:
	for ( int marker = _snp_number - 2; marker >= 0; marker-- ) {
		update_marker( first, second, individuals_ID, marker + 1, error_model, genohap, uniform, haplotypes, mutation_rate );
		double theta = theta_array[ marker ];
		if ( _block_flags[ marker ] ) { // inter edge
			alpha = _LeftMatrix[ marker ];
			int prev_marker = marker;
			int next_marker = marker + 1;
			int prev_cluster = _block_IDs[ prev_marker ];
			int next_cluster = _block_IDs[ next_marker ];
			uint32_t * p_w_prev = _hg._kjd_array[ prev_cluster ]->weigth;
			uint32_t * p_w_next = _hg._kjd_array[ next_cluster ]->weigth;

			n_haps = _hg._kjd_array[ prev_cluster ]->size;
			int p_haps = _hg._kjd_array[ next_cluster ]->size; // poster haplotypes
			double sum11 = 0;
			double p00 = 0, p01 = 0, p10 = 0, p11 = 0;
			vector<double> vec11(n_haps * n_haps, 0.);
			vector<double> part1, part2;

			vector<double> v00, v01, v10, v11;
			vector< pair<int,int> > r00, r01, r10, r11;
			
			for ( int k = 0; k < n_haps; k++ ) {
				part1.push_back((1. - theta) * 1. / double(p_w_prev[ k ]));
			}
			for ( int k = 0; k < p_haps; k++ ) {
				part2.push_back(theta * double(p_w_next[ k ]) / _haplotype_number);
			}
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					double C_Km_Km1_A = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( k ) ) {
						C_Km_Km1_A = double( _hg._kjds[ next_cluster ].prev_edges[ first ][ k ] );
					 }
					double C_Km_Km1_B = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( l ) ) {
						C_Km_Km1_B = double( _hg._kjds[ next_cluster ].prev_edges[ second ][ l ] );
					 }
					double C_Km_Km1_C = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( l ) ) {
						C_Km_Km1_C = double( _hg._kjds[ next_cluster ].prev_edges[ first ][ l ] );
					 }
					double C_Km_Km1_D = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( k ) ) {
						C_Km_Km1_D = double( _hg._kjds[ next_cluster ].prev_edges[ second ][ k ] );
					 }
					/*
					double tmp1 = (1. - theta) * C_Km_Km1_A / double(p_w_prev[ k ]) + theta * double(p_w_next[ first ]) / _haplotype_number;
					tmp1 *= (1. - theta) * C_Km_Km1_B / double(p_w_prev[ l ]) + theta * double(p_w_next[ second ]) / _haplotype_number;
					double tmp2 = (1. - theta) * C_Km_Km1_C / double(p_w_prev[ l ]) + theta * double(p_w_next[ first ]) / _haplotype_number;
					tmp2 *= (1. - theta) * C_Km_Km1_D / double(p_w_prev[ k ]) + theta * double(p_w_next[ second ]) / _haplotype_number;
					*/
				//	double tmp1 = ( part1[k] * C_Km_Km1_A + part2[first] ) * ( part1[l] * C_Km_Km1_B + part2[second] );
				//	double tmp2 = ( part1[l] * C_Km_Km1_C + part2[first] ) * ( part1[k] * C_Km_Km1_D + part2[second] );
					
				//	tmp1 *= alpha[ k * n_haps + l ];
				//	tmp2 *= alpha[ l * n_haps + k ];
				//	sum11 += (tmp1 + tmp2);
				//	vec11[k * n_haps + l] = tmp1;
				//	vec11[l * n_haps + k] = tmp2;
					
					// add some code to update recommbination rate.
				//	double checking = 0;
					p11 += part2[first] * part2[second] * alpha[ k * n_haps + l ];
					v11.push_back( part2[first] * part2[second] * alpha[ k * n_haps + l ] );
					r11.push_back( pair<int,int>(k, l) );
				//	cerr << "p11 = " << part2[first] * part2[second] * alpha[ k * n_haps + l ] << endl;
				//	checking += part2[first] * part2[second] * alpha[ k * n_haps + l ];
					
					p01 += part1[k] * C_Km_Km1_A * part2[second] * alpha[ k * n_haps + l ];
					v01.push_back( part1[k] * C_Km_Km1_A * part2[second] * alpha[ k * n_haps + l ] );
					r01.push_back( pair<int,int>(k, l) );
				//	checking += part1[k] * C_Km_Km1_A * part2[second] * alpha[ k * n_haps + l ];

					
					p10 += part2[first] * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ];
					v10.push_back( part2[first] * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ] );
					r10.push_back( pair<int,int>(k, l) );
				//	checking += part2[first] * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ];

					
					p00 += part1[k] * C_Km_Km1_A * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ];
					v00.push_back( part1[k] * C_Km_Km1_A * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ] );
					r00.push_back( pair<int,int>(k, l) );
				//	checking += part1[k] * C_Km_Km1_A * part1[l] * C_Km_Km1_B * alpha[ k * n_haps + l ];
					
					if ( k == l ) {
						continue;
					}

					p11 += part2[first] * part2[second] * alpha[ l * n_haps + k ];
					v11.push_back( part2[first] * part2[second] * alpha[ l * n_haps + k ] );
					r11.push_back( pair<int,int>(l, k) );
				//	checking += part2[first] * part2[second] * alpha[ l * n_haps + k ];


					p01 += part1[l] * C_Km_Km1_C * part2[second] * alpha[ l * n_haps + k ];
					v01.push_back( part1[l] * C_Km_Km1_C * part2[second] * alpha[ l * n_haps + k ] );
					r01.push_back( pair<int,int>(l, k) );
				//	checking += part1[l] * C_Km_Km1_C * part2[second] * alpha[ l * n_haps + k ];


					p10 += part1[k] * C_Km_Km1_D *  part2[first] * alpha[ l * n_haps + k ];
					v10.push_back( part1[k] * C_Km_Km1_D *  part2[first] * alpha[ l * n_haps + k ] );
					r10.push_back( pair<int,int>(l, k) );
				//	checking += part1[k] * C_Km_Km1_D *  part2[first] * alpha[ l * n_haps + k ];


					p00 += part1[l] * C_Km_Km1_C * part1[k] * C_Km_Km1_D * alpha[ l * n_haps + k ];
					v00.push_back( part1[l] * C_Km_Km1_C * part1[k] * C_Km_Km1_D * alpha[ l * n_haps + k ] );
					r00.push_back( pair<int,int>(l, k) );
				//	checking += part1[l] * C_Km_Km1_C * part1[k] * C_Km_Km1_D * alpha[ l * n_haps + k ];
				//	cerr << "sum = " << tmp1 + tmp2 << "\tchecking = " << checking << endl;
				}
				sum11 -= vec11[k * n_haps + k];
			}
		//	choice = uniform() * sum11;
		//	tmp = 0;
		//	int first0 = -1;
		//	int second0 = -1;
		//	for ( size_t g = 0; g < vec11.size(); g++ ) {
		//		tmp += vec11[ g ];
		//		if ( tmp > choice ) {
		//			first0 = g / n_haps;
		//			second0 = g % n_haps;
		//			break;
		//		}
		//	}

			/// select weather recombination events occur
			/*
			double C_Km_Km1_A = 0;
			 if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( first0 ) ) {
				C_Km_Km1_A = double ( _hg._kjds[ next_cluster ].prev_edges[ first ][ first0 ] );
			 }
			double C_Km_Km1_B = 0;
			 if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( second0 ) ) {
				C_Km_Km1_B = double ( _hg._kjds[ next_cluster ].prev_edges[ second ][ second0 ] );
			 }
			double prob_A_one = (1. - theta) * C_Km_Km1_A / double(p_w_prev[ first0 ]);
			double prob_A_two = theta * double(p_w_next[ first ]) / _haplotype_number;
			double prob_B_one = (1. - theta) * C_Km_Km1_B / double(p_w_prev[ second0 ]);
			double prob_B_two = theta * double(p_w_next[ second ]) / _haplotype_number;
			
			vector< double > probs;
			probs.push_back( prob_A_one * prob_B_one );
			probs.push_back( prob_A_one * prob_B_two + prob_A_two * prob_B_one );
			probs.push_back( prob_A_two * prob_B_two );

			sum = probs[ 0 ] + probs[ 1 ] + probs[ 2 ];
			crossover[ marker ].cross += (probs[ 1 ] + probs[ 2 ] * 2) / sum;
			*/
		/*
			choice = uniform() * sum;
			tmp = 0;
			int change = 0;
			for ( int k = 0; k < 3; k++ ) {
				tmp += probs[ k ];
				if ( tmp > choice ) {
					change = k;
					break;
				}
			}
		*/
			double p_all = p00 + p01 + p10 + p11;
			crossover[ marker ].cross += (p01 + p10 + 2 * p11) / p_all;
			crossover[ marker ].total += 2;
		//	crossover[ marker ].cross += change;

		//	assert( sum11 == p_all );

			double rc = uniform() * p_all;

		//	cerr << "n_haps = " << n_haps << endl;

			double ac = 0;
			if ( rc < p00 ) {
				for ( size_t k = 0; k < v00.size(); k++ ) {
					ac += v00[ k ];
					if ( ac >= rc ) {
						first = r00[ k ].first;
						second = r00[ k ].second;
						break;
					}
				}
		//		cerr << "ac = " << ac << "\tfirst = " << first << "\tsecond = " << second << endl;
				continue;
			}

			if ( rc < p00 + p10 ) {
				ac = p00;
				for ( size_t k = 0; k < v10.size(); k++ ) {
					ac += v10[ k ];
					if ( ac >= rc ) {
						first = r10[ k ].first;
						second = r10[ k ].second;
						break;
					}
				}
		//		cerr << "ac = " << ac << "\tfirst = " << first << "\tsecond = " << second << endl;
				continue;
			}
				
			if ( rc < p00 + p10 + p01 ) {
				ac = p00 + p10;
				for ( size_t k = 0; k < v01.size(); k++ ) {
					ac += v01[ k ];
					if ( ac >= rc ) {
						first = r01[ k ].first;
						second = r01[ k ].second;
						break;
					}
				}
		//		cerr << "ac = " << ac << "\tfirst = " << first << "\tsecond = " << second << endl;
				continue;
			}
			
			ac = p00 + p01 + p10;
			for ( size_t k = 0; k < v11.size(); k++ ) {
				ac += v11[ k ];
				if ( ac >= rc ) {
					first = r11[ k ].first;
					second = r11[ k ].second;
					break;
				}
			}
		//	cerr << "ac = " << ac << "\tfirst = " << first << "\tsecond = " << second << endl;
			
		//	first = first0;
		//	second = second0;
		} else {
			alpha = _LeftMatrix[ marker ];
			cluster_id = _block_IDs[ marker ];
			n_haps = _hg._kjd_array[ cluster_id ]->size;
			uint32_t * p_w = _hg._kjd_array[ cluster_id ]->weigth;

			double sum00( 0 ), sum10( 0 ), sum01( 0 ), sum11( 0 );
			vector<double> vec10, vec01, vec11;
			vector< int >  ids10, ids01, ids11;
			double both_change = theta * theta / (_haplotype_number * _haplotype_number);
			double one_change = (1. - theta) * theta / _haplotype_number;
			double no_change = (1. - theta) * (1. - theta);
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					tmp = alpha[ k * n_haps + l ] * double(p_w[ first ] * p_w[ second ]) * both_change;
					tmp *= (1 + (k != l));
					vec11.push_back( tmp );
					ids11.push_back( k * n_haps + l );
					sum11 += tmp;
					if ( l == second ) {
						tmp = alpha[ k * n_haps + l ] * double(p_w[ first ]) * one_change;
						vec10.push_back( tmp );
						sum10 += tmp;
						ids10.push_back( k * n_haps + l );
					}
					if ( k == second && l != k ) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ first ]) * one_change;
						vec10.push_back( tmp );
						sum10 += tmp;
						ids10.push_back( l * n_haps + k );
					}
					if ( k == first) {
						tmp = alpha[ k * n_haps + l ] * double(p_w[ second ]) * one_change;
						vec01.push_back( tmp );
						sum01 += tmp;
						ids01.push_back( k * n_haps + l );
					}
					if ( l == first && l != k ) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ second ]) * one_change;
						vec01.push_back( tmp );
						sum01 += tmp;
						ids01.push_back( l * n_haps + k );
					}
					if ( k == first && l == second) {
						sum00 += (alpha[ k * n_haps + l ] * no_change);
					}
					if ( l == first && k == second && l != k) {
						sum00 += (alpha[ l * n_haps + k ] * no_change);
					}
				}
			}

			sum = sum00 + sum01 + sum10 + sum11;

			int pesudo_first = -1;
			int pesudo_second = -1;
			crossover[ marker ].total += 2;
			crossover[ marker ].cross += (sum01 + sum10 + sum11 * 2) / sum;

			choice = uniform() * sum;
			if ( choice < sum00 ) {
				continue;
			}
			choice -= sum00;
			if ( choice < sum01 ) {
				tmp = 0;
				for ( size_t k = 0; k < vec01.size(); k++ ) {
					tmp += vec01[ k ];
					if ( tmp > choice ) {
						pesudo_first = ids01[ k ] / n_haps;
						pesudo_second = ids01[ k ] % n_haps;
						break;
					}
				}
			//	crossover[ marker ].cross++;
				second = pesudo_second;
				continue;
			}
			choice -= sum01;
			if ( choice < sum10 ) {
				tmp = 0;
				for ( size_t k = 0; k < vec10.size(); k++ ) {
					tmp += vec10[ k ];
					if ( tmp > choice ) {
						pesudo_first = ids10[ k ] / n_haps;
						pesudo_second = ids10[ k ] % n_haps;
						break;
					}
				}
			//	crossover[ marker ].cross++;
				first = pesudo_first;
				continue;
			}
			choice -= sum10;
			if( choice < sum11 ){
				tmp = 0;
				for ( size_t k = 0; k < vec11.size(); k++ ) {
					tmp += vec11[ k ];
					if ( tmp > choice ) {
						if (uniform() > 0.5) {
							pesudo_first = ids11[ k ] / n_haps;
							pesudo_second = ids11[ k ] % n_haps;
					} else {
							pesudo_first = ids11[ k ] % n_haps;
							pesudo_second = ids11[ k ] / n_haps;
						}
						break;
					}
				}
			//	crossover[ marker ].cross += 2;
				first = pesudo_first;
				second = pesudo_second;
			}
		}
	}
	update_marker( first, second, individuals_ID, 0, error_model, genohap, uniform, haplotypes, mutation_rate );
}

void _HMM_Algo::update_marker( int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
				variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, float * mutation_rate )
{
	int offset = _block_offsets[ marker ];
	int cluster_id = _block_IDs[ marker ];
	ele * p_p = _hg._kjd_array[ cluster_id ]->piece;
	uint8_t * p_first = p_p[ first ].p;
	uint8_t * p_second = p_p[ second ].p;
	int base1 = p_first[ offset ];
	int base2 = p_second[ offset ];

	float * glikes = &genohap.glf[individuals_ID * _indiv_ost + (marker + _start) * 3 ];
	double posterior_11 = _errors[ marker ][ base1 + base2 ][ 0 ] * glikes[ 0 ];
	double posterior_12 = _errors[ marker ][ base1 + base2 ][ 1 ] * glikes[ 1 ];
	double posterior_22 = _errors[ marker ][ base1 + base2 ][ 2 ] * glikes[ 2 ];
	double sum = posterior_11 + posterior_12 + posterior_22;
	posterior_11 /= sum;
	posterior_22 /= sum;
	posterior_12 /= sum;
	
	double r = uniform();
	
	if ( r < posterior_11 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 0;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 0;
	} else if ( r < posterior_11 + posterior_22 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 1;
	} else if ( base1 != base2 ) {
		float rate = mutation_rate[ marker ];
		if ( uniform() < rate * rate / ( rate * rate + (1 - rate) * (1 - rate) ) ) {
			base1 = !base1;
			base2 = !base2;
		}
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = base1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = base2;
	} else {
		bool bit = (uniform() > 0.5) ? 1 : 0;
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = bit;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = bit ^ 1;
	}

	int imputed1 = haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ];
	int imputed2 = haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ];
	int differences = abs(base1 - imputed1) + abs(base2 - imputed2);
	error_model[ marker ].mismatch += differences;
	error_model[ marker ].match += 2 - differences;
}

void _HMM_Algo::update_haplotype_max( Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
			Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, \
			float * mutation_rate )
{
	float * alpha = _LeftMatrix[ _snp_number - 1 ];
	float max_posterior = 0.;
	int first = -1;
	int second = -1;
	int cluster_id = _block_IDs[ _snp_number - 1 ];
	int n_haps = _hg._kjd_array[ cluster_id ]->size;
	for ( int i = 0; i < n_haps; i++ ) {
		for ( int j = i; j < n_haps; j++ ) {
			double cur_prob = alpha[ i * n_haps + j ] * (1 + (i != j));
			if ( max_posterior < cur_prob ) {
				first = i;
				second = j;
				max_posterior = cur_prob;
			} else if (max_posterior >= cur_prob && max_posterior <= cur_prob) {
				if (0.5 < uniform()) {
					first = i;
					second = j;
					max_posterior = cur_prob;
				}
			}
		}
	} // back samples and select the haplotype pair with maximum probability

	double max = 0.;
	for ( int marker = _snp_number - 2; marker >= 0; marker-- ) {
		update_marker_max( first, second, individuals_ID, marker + 1, error_model, genohap, uniform, haplotypes, mutation_rate );
		double theta = theta_array[ marker ];
		max = 0.;
		if ( _block_flags[ marker ] ) { // inter edge
			alpha = _LeftMatrix[ marker ];
			int prev_marker = marker;
			int next_marker = marker + 1;
			int prev_cluster = _block_IDs[ prev_marker ];
			int next_cluster = _block_IDs[ next_marker ];
			uint32_t * p_w_prev = _hg._kjd_array[ prev_cluster ]->weigth;
			uint32_t * p_w_next = _hg._kjd_array[ next_cluster ]->weigth;

			int first0 = -1, second0 = -1;

			n_haps = _hg._kjd_array[ prev_cluster ]->size;
			int p_haps = _hg._kjd_array[ next_cluster ]->size; // poster haplotypes
			vector<double> part1, part2;
			for ( int k = 0; k < n_haps; k++ ) {
				part1.push_back((1. - theta) * 1. / double(p_w_prev[ k ]));
			}
			for ( int k = 0; k < p_haps; k++ ) {
				part2.push_back(theta * double(p_w_next[ k ]) / _haplotype_number);
			}
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					double C_Km_Km1_A = 0;
					if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( k ) ) {
						C_Km_Km1_A = double ( _hg._kjds[ next_cluster ].prev_edges[ first ][ k ] );
					}
					double C_Km_Km1_B = 0;
					if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( l ) ) {
						C_Km_Km1_B = double ( _hg._kjds[ next_cluster ].prev_edges[ second ][ l ] );
					}
					double C_Km_Km1_C = 0;
					if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( l ) ) {
						C_Km_Km1_C = double( _hg._kjds[ next_cluster ].prev_edges[ first ][ l ] );
					}
					double C_Km_Km1_D = 0;
					if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( k ) ) {
						C_Km_Km1_D = double( _hg._kjds[ next_cluster ].prev_edges[ second ][ k ] );
					}
					double tmp1 = ( part1[k] * C_Km_Km1_A + part2[first] ) * ( part1[l] * C_Km_Km1_B + part2[second] );
					double tmp2 = ( part1[l] * C_Km_Km1_C + part2[first] ) * ( part1[k] * C_Km_Km1_D + part2[second] );
					tmp1 *= alpha[ k * n_haps + l ];
					tmp2 *= alpha[ l * n_haps + k ];
					if (tmp1 > max) {
						first0 = k;
						second0 = l;
						max = tmp1;
					}
					if (tmp2 > max) {
						first0 = l;
						second0 = k;
						max = tmp2;
					}
				}
			}
			/// select weather recombination events occur
			double C_Km_Km1_A = 0;
			if ( _hg._kjds[ next_cluster ].prev_edges[ first ].count( first0 ) ) {
				C_Km_Km1_A = double ( _hg._kjds[ next_cluster ].prev_edges[ first ][ first0 ] );
			}
			double C_Km_Km1_B = 0;
			if ( _hg._kjds[ next_cluster ].prev_edges[ second ].count( second0 ) ) {
				C_Km_Km1_B = double ( _hg._kjds[ next_cluster ].prev_edges[ second ][ second0 ] );
			}
			crossover[ marker ].total += 2;
			if (C_Km_Km1_A >= 0 && C_Km_Km1_A <= 0) crossover[ marker ].cross++;
			if (C_Km_Km1_B >= 0 && C_Km_Km1_A <= 0) crossover[ marker ].cross++;

			first = first0;
			second = second0;
		} else {
			alpha = _LeftMatrix[ marker ];
			cluster_id = _block_IDs[ marker ];
			n_haps = _hg._kjd_array[ cluster_id ]->size;
			uint32_t * p_w = _hg._kjd_array[ cluster_id ]->weigth;
			double tmp = 0.;
			double sum[4] = {0.};
			vector< double > vec10, vec01, vec11;
			vector< int >    ids10, ids01, ids11;
			double both_change = theta * theta / (_haplotype_number * _haplotype_number);
			double one_change = (1. - theta) * theta / _haplotype_number;
			double no_change = (1. - theta) * (1. - theta);
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					tmp = alpha[ k * n_haps + l ] * double(p_w[ first ] * p_w[ second ]) * both_change;
					tmp *= (1 + (k != l));
					vec11.push_back( tmp );
					ids11.push_back( k * n_haps + l );
					sum[3] += tmp;
					if ( l == second ) { // recombinate on first haplotype
						tmp = alpha[ k * n_haps + l ] * double(p_w[ first ]) * one_change;
						vec10.push_back( tmp );
						ids10.push_back( k * n_haps + l );
						sum[1] += tmp;
					}
					if ( k == second && l != k ) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ first ]) * one_change;
						vec10.push_back( tmp );
						ids10.push_back( l * n_haps + k );
						sum[1] += tmp;
					}
					if ( k == first ) { // recombinate on second haplotype
						tmp = alpha[ k * n_haps + l ] * double(p_w[ second ]) * one_change;
						vec01.push_back( tmp );
						sum[2] += tmp;
						ids01.push_back( k * n_haps + l );
					}
					if (l == first && l != k) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ second ]) * one_change;
						vec01.push_back( tmp );
						sum[2] += tmp;
						ids01.push_back( l * n_haps + k );
					}
					if ( k == first && l == second ) {
						sum[0] += (alpha[ k * n_haps + l ] * no_change);
					}
					if (l == first && k == second && l != k) {
						sum[0] += (alpha[ l * n_haps + k ] * no_change);
					}
				}
			}
			int choice = -1;
			max = 0.;
			for (int i = 0; i != 4; ++i) {
				if (sum[i] > max) {
					choice = i;
					max = sum[i];
				} else if (sum[i] >= max && sum[i] <= max) {
					if (0.5 > uniform()) {
						choice = i;
						max = sum[i];
					}
				}
			}
			int pesudo_first = -1;
			int pesudo_second = -1;
			crossover[ marker ].total += 2;
			if ( choice == 0 ) {
				continue;
			} else if ( choice == 1 ) {
				max = 0.;
				for ( size_t k = 0; k < vec10.size(); k++ ) {
					if ( vec10[k] > max ) {
						pesudo_first = ids10[ k ] / n_haps;
						pesudo_second = ids10[ k ] % n_haps;
						max = vec10[k];
					} else if (vec10[k] >= max && vec10[k] <= max) {
						if (0.5 > uniform()) {
							pesudo_first = ids10[ k ] / n_haps;
							pesudo_second = ids10[ k ] % n_haps;
							max = vec10[k];
						}
					}
				}
				crossover[ marker ].cross++;
				first = pesudo_first;
				if ( pesudo_second != second ) {
					cerr << "[_HMM_Algo::update_haplotype] pesudo_first != first.\n";
					abort();
				}
				continue;
			} else if ( choice == 2 ) {
				max = 0.;
				for ( size_t k = 0; k < vec01.size(); k++ ) {
					if ( vec01[k] > max ) {
						pesudo_first = ids01[ k ] / n_haps;
						pesudo_second = ids01[ k ] % n_haps;
						max = vec01[k];
					} else if (vec01[k] > max ) {
						if (0.5 > uniform()) {
							pesudo_first = ids01[ k ] / n_haps;
							pesudo_second = ids01[ k ] % n_haps;
							max = vec01[k];
						}
					}
				}
				crossover[ marker ].cross++;
				second = pesudo_second;
				if ( pesudo_first != first ) {
					cerr << "[_HMM_Algo::update_haplotype] pesudo_second != second.\n";
					abort();
				}
				continue;
			} else if (choice == 3) {
				max = 0.;
				for ( size_t k = 0; k < vec11.size(); k++ ) {
					if ( vec11[k] > max ) {
						pesudo_first = ids11[ k ] / n_haps;
						pesudo_second = ids11[ k ] % n_haps;
						max = vec11[k];
					} else if (vec11[k] >= max && vec11[k] <= max) {
						if (0.5 > uniform()) {
							pesudo_first = ids11[ k ] / n_haps;
							pesudo_second = ids11[ k ] % n_haps;
							max = vec11[k];
						}
					}
				}
				crossover[ marker ].cross += 2;
				if (uniform() > 0.5) {
					first = pesudo_first;
					second = pesudo_second;
				} else {
					first = pesudo_second;
					second = pesudo_first;
				}
			}
		}
	}
	update_marker_max( first, second, individuals_ID, 0, error_model, genohap, uniform, haplotypes, mutation_rate );
}

void _HMM_Algo::update_marker_max( int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
				   variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, \
				   float * mutation_rate )
{
	int offset = _block_offsets[ marker ];
	int cluster_id = _block_IDs[ marker ];
	ele * p_p = _hg._kjd_array[ cluster_id ]->piece;
	uint8_t * p_first = p_p[ first ].p;
	uint8_t * p_second = p_p[ second ].p;
	int base1 = p_first[ offset ];
	int base2 = p_second[ offset ];

	float * glikes = &genohap.glf[individuals_ID * _indiv_ost + (marker + _start) * 3];
	
	vector< double > posteriors;
	posteriors.push_back( _errors[ marker ][ base1 + base2 ][ 0 ] * glikes[ 0 ] );
	posteriors.push_back( _errors[ marker ][ base1 + base2 ][ 1 ] * glikes[ 1 ] );
	posteriors.push_back( _errors[ marker ][ base1 + base2 ][ 2 ] * glikes[ 2 ] );

	int id = -1;
	double max_posteriors = 0;
	for ( size_t i = 0; i < posteriors.size(); i++ ) {
		if ( posteriors[ i ] > max_posteriors ) {
			max_posteriors = posteriors[ i ];
			id = i;
		}
	}
	
	if ( id == 0 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 0;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 0;
	} else if ( id == 2 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 1;
	} else if ( base1 != base2 ) {
		double rate = mutation_rate[ marker ];
		if ( uniform() < rate * rate / ( rate * rate + (1 - rate) * (1 - rate) ) ) {
			base1 = !base1;
			base2 = !base2;
		}
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = base1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = base2;
	} else {
		bool bit = (uniform() > 0.5) ? 1 : 0;
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = bit;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = bit ^ 1;
	}

	int imputed1 = haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ];
	int imputed2 = haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ];
	int differences = abs(base1 - imputed1) + abs(base2 - imputed2);
	error_model[ marker ].mismatch += differences;
	error_model[ marker ].match += 2 - differences;
}


void _HMM_Algo::check_haplotypes( Haplotype * haplotypes )
{
	for ( int i = 0; i < _individuals; i++ ) {
		cerr << "I" << i << " first  " << haplotypes[ i * 2 + 0 ];
		cerr << "I" << i << " second " << haplotypes[ i * 2 + 1 ];
		for ( int j = 0; j < _snp_number; j++ ) {
			if ( int(haplotypes[ i * 2 + 0 ].haplotype[ j ]) > 1 ) {
				abort();
			}
			if ( int(haplotypes[ i * 2 + 1 ].haplotype[ j ]) > 1 ) {
				abort();
			}
		}
	}
}


_HMM_Algo::~_HMM_Algo()
{
	delete [] _start_array;
	delete [] _end_array;
	delete [] _block_IDs;
	delete [] _block_offsets;
	delete [] _block_flags;
	for ( int i = 0; i < _snp_number; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			delete [] _errors[ i ][ j ];
		}
		delete [] _errors[ i ];
	}
	delete [] _errors;
	delete [] _marginal_i;
	delete [] _marginal_j;
}

void _HMM_Algo::prior()
{
	float * _tm = _LeftMatrix[ 0 ];
	int n_hidden_states = _hg._kjd_array[ 0 ]->size;
	uint32_t * p_w = _hg._kjd_array[ 0 ]->weigth;

	for ( int i = 0; i < n_hidden_states; i++ ) {
		float fi = static_cast< float >( double(p_w[ i ]) / _haplotype_number );
		for ( int j = i; j < n_hidden_states; j++ ) {
			_tm[i * n_hidden_states + j] = fi * static_cast< float >( double(p_w[ j ]) / _haplotype_number );
		}
	}
}

void _HMM_Algo::match_data( int marker, int individuals_ID, _GenoHap & genohap )
{
	int cluster_id = _block_IDs[ marker ];
	int n_haps = _hg._kjd_array[ cluster_id ]->size;
	int offset = _block_offsets[ marker ];
	float * _tm = _LeftMatrix[ marker ];
	ele * p_p = _hg._kjd_array[ cluster_id ]->piece;

	double conditional_probs[3];
	for (int i = 0; i < 3; i++)
		conditional_probs[i] = _errors[ marker ][ i ][ 0 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 0 ] + \
				       _errors[ marker ][ i ][ 1 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 1 ] + \
				       _errors[ marker ][ i ][ 2 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 2 ];

	double factors[ 2 ];
	for ( int i = 0; i < n_haps; i++ ) {
		uint8_t * p_i = p_p[ i ].p;
		factors[ 0 ] = conditional_probs[ p_i[ offset ] ];
		factors[ 1 ] = conditional_probs[ p_i[ offset ] + 1 ];
		for ( int j = i; j < n_haps; j++ ) {
			uint8_t * p_j = p_p[ j ].p;
			_tm[i * n_haps + j] *= static_cast< float >( factors[ p_j[ offset ] ] );
		}
	}
}

void _HMM_Algo::prepare_pointors( int prev_marker, int next_marker )
{
	int prev_cluster = _block_IDs[ prev_marker ];
	int next_cluster = _block_IDs[ next_marker ];
	int prev_size = _hg._kjd_array[ prev_cluster ]->size;
	int next_size = _hg._kjd_array[ next_cluster ]->size;
	
	// all prev pointors
	ele * prev_p = _hg._kjd_array[ prev_cluster ]->piece;
	prev_piece = new uint8_t * [ prev_size ];
	for ( int i = 0; i < prev_size; i++ ) {
		prev_piece[ i ] = prev_p[ i ].p;
	}
	prev_weigth = _hg._kjd_array[ prev_cluster ]->weigth;

	// all next pointors
	ele * next_p = _hg._kjd_array[ next_cluster ]->piece;
	next_piece = new uint8_t * [ next_size ];
	for ( int i = 0; i < next_size; i++ ) {
		next_piece[ i ] = next_p[ i ].p;
	}
	next_weigth = _hg._kjd_array[ next_cluster ]->weigth;
}

void _HMM_Algo::destory_pointors( int prev_marker, int next_marker )
{
	delete [] prev_piece;
	delete [] next_piece;
	prev_piece = NULL;
	next_piece = NULL;
	prev_weigth = NULL;
	next_weigth = NULL;
}

void _HMM_Algo::transpose( int prev_marker, int next_marker, float * theta_array, _GenoHap & genohap, int individuals_ID )
{
	clear_matrix();
	float * source = _LeftMatrix[ prev_marker ];
	float * dest   = _LeftMatrix[ next_marker ];
	float theta = theta_array[ prev_marker ];

	if ( _block_flags[ prev_marker ] ) {
		int prev_cluster = _block_IDs[ prev_marker ];
		int next_cluster = _block_IDs[ next_marker ];
		int prev_size = _hg._kjd_array[ prev_cluster ]->size;
		int next_size = _hg._kjd_array[ next_cluster ]->size;

		prepare_pointors( prev_marker, next_marker );

		double sum = 0;
		for ( int i = 0; i < prev_size; i++ ) {
			sum += source[ i * prev_size + i ];
			_marginal_i[ i ] += source[ i * prev_size + i ];
			for ( int j = i + 1; j < prev_size; j++ ) {
				sum += source[ i * prev_size + j ] * 2;
				_marginal_i[ i ] += source[ i * prev_size + j ];
				_marginal_i[ j ] += source[ j * prev_size + i ];
			}
		}
		for (int i = 0; i < prev_size; i++ ) {
			_marginal_j[ i ] = _marginal_i[ i ];
		}

		vector< vector< pair< int, double > > > freqs;
		vector< pair< int, double > > freq;
		// vector< double > next_freqs;
		// vector< pair< int, double > > & freq_k;
		// vector< pair< int, double > > & freq_l;

		double no_change = (1. - theta) * (1. - theta);
		double one_change = (1. - theta) * theta;
		double two_change = theta * theta * sum;
		for ( int k = 0; k < next_size; k++ ) {
			// next_freqs.push_back(double(next_weigth[ k ]) / _haplotype_number);
			freq.clear();
			map< int, int >::iterator iter1 = _hg._kjds[ next_cluster ].prev_edges[ k ].begin();
			map< int, int >::iterator iter2 = _hg._kjds[ next_cluster ].prev_edges[ k ].end();
			for ( ; iter1 != iter2; iter1++ ) {
				double f = double(iter1->second) / double(prev_weigth[ iter1->first ]);
				freq.push_back( pair< int, double >( iter1->first, f ) );
			}
			freqs.push_back(freq);
		}

		for ( int k = 0; k < next_size; k++ ) {
			for ( int l = k; l < next_size; l++ ) {
				// freq_k.clear();
				// freq_l.clear();
				double fl = /* next_freqs[k]; */ double(next_weigth[ l ]) / _haplotype_number;
				double fk = /* next_freqs[l]; */ double(next_weigth[ k ]) / _haplotype_number;
				// k
				vector< pair< int, double > > & freq_k = freqs[k];
				/*
				map< int, int >::iterator iter1 = _hg._kjds[ next_cluster ].prev_edges[ k ].begin();
				map< int, int >::iterator iter2 = _hg._kjds[ next_cluster ].prev_edges[ k ].end();
				for ( ; iter1 != iter2; iter1++ ) {
					double freq = double(iter1->second) / double(prev_weigth[ iter1->first ]);
					freq_k.push_back( pair< int, double >( iter1->first, freq ) );
				}
				*/
				// l
				/*
				iter1 = _hg._kjds[ next_cluster ].prev_edges[ l ].begin();
				iter2 = _hg._kjds[ next_cluster ].prev_edges[ l ].end();
				for ( ; iter1 != iter2; iter1++ ) {
					double freq = double(iter1->second) / double(prev_weigth[ iter1->first ]);
					freq_l.push_back( pair< int, double>( iter1->first, freq ) );
				}
				*/
				vector< pair< int, double > > & freq_l = freqs[l];

				double prob = 0;
				double prob_t = 0;
				// no recombination
				for ( size_t n_freq_k = 0; n_freq_k < freq_k.size(); n_freq_k++ ) {
					for ( size_t n_freq_l = 0; n_freq_l < freq_l.size(); n_freq_l++ ) {
						int id_i = freq_k[ n_freq_k ].first;
						int id_j = freq_l[ n_freq_l ].first;
						prob_t += source[ id_i * prev_size + id_j ] * freq_k[ n_freq_k ].second * \
							freq_l[ n_freq_l ].second;
					}
				}
				prob_t *= no_change;
				prob += prob_t;
				prob_t = 0;
				// j change
				for ( size_t n_freq_k = 0; n_freq_k < freq_k.size(); n_freq_k++ ) {
					int id_k = freq_k[ n_freq_k ].first;
					prob_t += _marginal_i[ id_k ] * freq_k[ n_freq_k ].second;
				}
				prob_t *= fl * one_change;
				prob += prob_t;
				prob_t = 0;
				// k change
				for ( size_t n_freq_l = 0; n_freq_l < freq_l.size(); n_freq_l++ ) {
					int id_l = freq_l[ n_freq_l ].first;
					prob_t += _marginal_j[ id_l ] * freq_l[ n_freq_l ].second;
				}
				prob_t *= fk * one_change;
				prob += prob_t;
				prob_t = 0;
				// k & j both changes
				prob_t = two_change * fk * fl;
				prob += prob_t;
				_LeftMatrix[ next_marker ][ k * next_size + l ] = static_cast< float >( prob );
			}
		}

		destory_pointors( prev_marker, next_marker );
	} else {
		int prev_cluster = _block_IDs[ prev_marker ];
		int n_haps = _hg._kjd_array[ prev_cluster ]->size;
		
		float * gl = &genohap.glf[ individuals_ID * _indiv_ost + (_start + next_marker) * 3 ];
	
		double sum = 0.;
		for ( int i = 0; i < n_haps; i++ ) {
			sum += source[ i * n_haps + i ];
			_marginal_i[ i ] += source[ i * n_haps + i ];
			for ( int j = i + 1; j < n_haps; j++ ) {
				sum += source[ i * n_haps + j ] * 2;
				_marginal_i[ i ] += source[ i * n_haps + j ];
				_marginal_i[ j ] += source[ j * n_haps + i ];
			}
		}
		for (int i = 0; i < n_haps; i++ ) {
			_marginal_j[ i ] = _marginal_i[ i ];
		}

		double no_change = (1. - theta) * (1. - theta);
		double one_change = (1. - theta) * theta / _haplotype_number;
		double two_change = sum * theta * theta  / (_haplotype_number * _haplotype_number);

		uint32_t * p_w_prev = _hg._kjd_array[ prev_cluster ]->weigth;

		for ( int i = 0; i < n_haps; i++ ) {
			for ( int j = i; j < n_haps; j++ ) {
				dest[i * n_haps + j] = source[i * n_haps + j] * no_change + \
					_marginal_i[ i ] * one_change * p_w_prev[ j ] + \
					_marginal_j[ j ] * one_change * p_w_prev[ i ] + \
					two_change * p_w_prev[ i ] * p_w_prev[ j ];
			}
		}
	}
}

void _HMM_Algo::scale( int marker )
{
	int cluster_id = _block_IDs[ marker ];
	int n_haps = _hg._kjd_array[ cluster_id ]->size;
	float max_prob = 0;
	for ( int i = 0; i < n_haps; i++ ) {
		for ( int j = i; j < n_haps; j++ ) {
			if ( _LeftMatrix[ marker ][ i * n_haps + j ] > max_prob ) {
				max_prob = _LeftMatrix[ marker ][ i * n_haps + j ];
			}
		}
	}

	double multiplier = 1. / max_prob;
	for ( int i = 0; i < n_haps; i++ ) {
		for ( int j = i; j < n_haps; j++ ) {
			_LeftMatrix[ marker ][ i * n_haps + j ] *= multiplier;
			_LeftMatrix[ marker ][ j * n_haps + i ] = _LeftMatrix[ marker ][ i * n_haps + j ];
		}
	}
}

void _HMM_Algo::process( float * theta_array, _GenoHap & genohap, int individuals_ID )
{
	prior();
	match_data( 0, individuals_ID, genohap );
	scale( 0 );

	for ( int marker = 1; marker < _snp_number; marker++ ) {
		int prev_marker = marker - 1;
		int next_marker = marker;
		transpose( prev_marker, next_marker, theta_array, genohap, individuals_ID );
		match_data( next_marker, individuals_ID, genohap );
		scale( next_marker );
	}
}

/**
 * Code for huge sample size algorithm
 * 12,05,2013
*/

void _HMM_Algo::shared_haplotypes( ele & hap, ele & state, vector< _XK > & _XKs, int id, bool output )
{
	assert( hap.length == state.length );
	int n_site = hap.length;
	int c = 1;
	bool new_spin, old_spin;
	if ( hap.p[ 0 ] == state.p[ 0 ] ) {
		new_spin = old_spin = true;
	} else {
		new_spin = old_spin = false;
	}
	if ( output )
		cerr << "hap\t" << id << "\t|";
	for ( int i = 1; i < n_site; i++ ) {
		new_spin = hap.p[ i ] == state.p[ i ] ? true : false;
		if ( state.p[ i ] ) {
			if ( output )
				cerr << "*";
		} else {
			if ( output )
				cerr << ".";
		}
		if ( new_spin == old_spin ) {
			c++;
		} else {
			if ( ! new_spin && c >= 2 ) {
				_XK xk( i - c, i, id );
				_XKs.push_back( xk );
			}
			old_spin = new_spin;
			c = 1;
		}
	}
	if (output )
		cerr << "\n";
	if ( new_spin && c >= 2 ) {
		_XK xk( n_site - c, n_site, id );
		_XKs.push_back( xk );
	}
}

void _HMM_Algo::shared_haplotypes( ele & hap, ele & state, vector< pair< int, int > > & blks, int & ibd_distance )
{
	blks.clear();
	ibd_distance = 0;
	int  n_site = hap.length;
	assert( hap.length == state.length );
	int  c = 1;
	bool new_spin, old_spin;
	if ( hap.p[ 0 ] == state.p[ 0 ] ) {
		new_spin = old_spin = true;
	} else {
		new_spin = old_spin = false;
	}
	for ( int i = 1; i < n_site; i++ ) {
		if ( hap.p[ i ] == state.p[ i ] ) {
			new_spin = true;
		} else {
			new_spin = false;
		}
		if ( new_spin == old_spin ) {
			c++;
		} else {
			if ( ! new_spin ) {
				if ( c >= 5 ) {
					blks.push_back( pair< int, int >( i - c, i ) );
					if ( c > ibd_distance ) {
						ibd_distance = c;
					}
				}
			}
			old_spin = new_spin;
			c = 1;
		}
	}
	if ( new_spin ) {
		if ( c >= 5 ) {
			blks.push_back( pair< int, int >( n_site - c, n_site ) );
			if ( c > ibd_distance ) {
				ibd_distance = c;
			}
		}
	}
}

void _HMM_Algo::prune_states_find( int individuals_ID, Haplotype * haplotypes, int maximum )
{
	int n_chunks = maximum / 10;

	for ( size_t i = 0; i < prune_states.size(); i++ ) {
		prune_states[ i ].clear();
	}
	prune_states.clear();
	
	int locate = 0;
	for ( int i = 0; i < _hg.n_kjds; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		
		ele * p = _hg._kjd_array[ i ]->piece;
		ele hap1( n_site );
		memcpy( hap1.p, haplotypes[ individuals_ID * 2 + 0 ].haplotype + locate, n_site );
		ele hap2( n_site );
		memcpy( hap2.p, haplotypes[ individuals_ID * 2 + 1 ].haplotype + locate, n_site );
		locate += n_site;
	
		if ( n_haps < maximum ) {
			vector< int > states_t;
			for ( int k = 0; k < n_haps; k++ ) {
				states_t.push_back( k );
			}
			prune_states.push_back( states_t );
			continue;
		}
	
		int n_gaps = n_site / n_chunks;
		int n_remain = n_site % n_chunks;
		vector< int > chunks;
		int P = 0;
		for ( int g = 0; g < n_chunks; g++ ) {
			chunks.push_back( P );
			if ( n_remain > 0 ) {
				P += n_gaps + 1;
				n_remain--;
			} else {
				P += n_gaps;
			}
		}
		chunks.push_back( n_site + 1 );
		vector< pair< int, int > > forward_chunks;
		vector< pair< int, int > > backward_chunks;
		if ( n_chunks % 2 == 1 ) {
			for ( int g = 0; g < n_chunks / 2 + 1; g++ ) {
				forward_chunks.push_back( pair< int, int >( chunks[ g ], chunks[ g + 1 ] ) );
			}
			int chunk_size = chunks.size();
			for ( int g = 0; g < n_chunks / 2; g++ ) {
				backward_chunks.push_back( pair< int, int>( chunks[ chunk_size - g - 2 ], chunks[ chunk_size - g - 1 ] ) );
			}
		} else {
			for ( int g = 0; g < n_chunks / 2; g++ ) {
				forward_chunks.push_back( pair< int, int >( chunks[ g ], chunks[ g + 1 ] ) );
			}
			int chunk_size = chunks.size();
			for ( int g = 0; g < n_chunks / 2; g++ ) {
				backward_chunks.push_back( pair< int, int>( chunks[ chunk_size - g - 2 ], chunks[ chunk_size - g  - 1 ] ) );
			}
		}
	
		vector< _XK > _xks1, _xks2;
		for ( int j = 0; j < n_haps; j++ ) {
			shared_haplotypes( hap1, p[ j ], _xks1, j, false );
			shared_haplotypes( hap2, p[ j ], _xks2, j, false );
		}

		// choose hidden states.
		sort( _xks1.begin(), _xks1.end(), _xk_compare_() );
		sort( _xks2.begin(), _xks2.end(), _xk_compare_() );
	
		map< int, int > hap1_map, hap2_map;
		size_t _xks1_index = 0, _xks2_index = 0;
		vector< _XK > _xks_t;
		for ( size_t k = 0; k < forward_chunks.size(); k++ ) {
			int pseudo_begin = forward_chunks[ k ].first;
			int pseudo_end   = forward_chunks[ k ].second;
			for ( ; _xks1_index < _xks1.size(); _xks1_index++ ) {
				if ( _xks1[ _xks1_index ].start >= pseudo_begin && _xks1[ _xks1_index ].start < pseudo_end ) {
					_xks_t.push_back( _xks1[ _xks1_index ] );
				} else if ( _xks1[ _xks1_index ].start >= pseudo_end ) {
					break;
				} else {
					;
				}
			}
			
			sort( _xks_t.begin(), _xks_t.end(), _xk_sort_() );
			if ( _xks_t.size() >= 5 ) {
				for ( size_t get = 0; get < 5; get++ ) {
					hap1_map[ _xks_t[ get ].id ]++;
				}
			} else {
				for ( size_t get = 0; get < _xks_t.size(); get++ ) {
					hap1_map[ _xks_t[ get ].id ]++;
				}
			}

			_xks_t.clear();
			for ( ; _xks2_index < _xks2.size(); _xks2_index++ ) {
				if ( _xks2[ _xks2_index ].start >= pseudo_begin && _xks2[ _xks2_index ].start < pseudo_end ) {
					_xks_t.push_back( _xks2[ _xks2_index ] );
				} else if ( _xks2[ _xks2_index ].start >= pseudo_end ) {
					break;
				} else {
					;
				}
			}
			
			sort( _xks_t.begin(), _xks_t.end(), _xk_sort_() );
			if ( _xks_t.size() >= 5 ) {
				for ( size_t get = 0; get < 5; get++ ) {
					hap2_map[ _xks_t[ get ].id ]++;
				}	
			} else {
				for ( size_t get = 0; get < _xks_t.size(); get++ ) {
					hap2_map[ _xks_t[ get ].id ]++;
				}
			}
			
			_xks_t.clear();
		}

		sort( _xks1.begin(), _xks1.end(), _xk_reverse_compare_() );
		sort( _xks2.begin(), _xks2.end(), _xk_reverse_compare_() );

		_xks1_index = 0; _xks2_index = 0;
		for ( size_t k = 0; k < backward_chunks.size(); k++ ) {
			int pseudo_begin = backward_chunks[ k ].first;
			int pseudo_end   = backward_chunks[ k ].second;
			for ( ; _xks1_index < _xks1.size(); _xks1_index++ ) {
				if ( _xks1[ _xks1_index ].end >= pseudo_begin && _xks1[ _xks1_index ].end < pseudo_end ) {
					_xks_t.push_back( _xks1[ _xks1_index ] );
				} else if ( _xks1[ _xks1_index ].end < pseudo_begin ) {
					break;
				} else {
					;
				}
			}
			
			sort( _xks_t.begin(), _xks_t.end(), _xk_sort_() );
			if ( _xks_t.size() >= 5 ) {
				for ( size_t get = 0; get < 5; get++ ) {
					hap1_map[ _xks_t[ get ].id ]++;
				}
			} else {
				for ( size_t get = 0; get < _xks_t.size(); get++ ) {
					hap1_map[ _xks_t[ get ].id ]++;
				}
			}
			_xks_t.clear();

			for ( ; _xks2_index < _xks2.size(); _xks2_index++ ) {
				if ( _xks2[ _xks2_index ].end >= pseudo_begin && _xks2[ _xks2_index ].end < pseudo_end ) {
					_xks_t.push_back( _xks2[ _xks2_index ] );
				} else if ( _xks2[ _xks2_index ].end < pseudo_begin ) {
					break;
				} else {
					;
				}
			}
			
			sort( _xks_t.begin(), _xks_t.end(), _xk_sort_() );
			if ( _xks_t.size() >= 5 ) {
				for ( size_t get = 0; get < 5; get++ ) {
					hap2_map[ _xks_t[ get ].id ]++;
				}
			} else {
				for ( size_t get = 0; get < _xks_t.size(); get++ ) {
					hap2_map[ _xks_t[ get ].id ]++;
				}
			}
			_xks_t.clear();
		}
		vector< int > states;
		map< int, int >::iterator it1 = hap1_map.begin();
		map< int, int >::iterator it2 = hap1_map.end();
		for ( ; it1 != it2; it1++ ) {
			states.push_back( it1->first );
		}
		
		it1 = hap2_map.begin();
		it2 = hap2_map.end();
		for ( ; it1 != it2; it1++ ) {
			if ( hap1_map.count( it1->first ) ) {
				continue;
			} else {
				states.push_back( it1->first );
			}
		}
		prune_states.push_back( states );
		if ( states.size() == 0 ) {
			abort();
		}
	}
}

void _HMM_Algo::prune_states_find( int individuals_ID, int maximum, Haplotype * haplotypes )
{
	for ( size_t i = 0; i < prune_states.size(); i++ ) {
		prune_states[ i ].clear();
	}
	prune_states.clear();
	
	int locate = 0;
	for ( int i = 0; i < _hg.n_kjds; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		
		ele * p_p = _hg._kjd_array[ i ]->piece;
		ele hap1( n_site );
		memcpy( hap1.p, haplotypes[ individuals_ID * 2 + 0 ].haplotype + locate, n_site );
		ele hap2( n_site );
		memcpy( hap2.p, haplotypes[ individuals_ID * 2 + 1 ].haplotype + locate, n_site );
		locate += n_site;
		
		if ( n_haps < maximum ) {
                        vector< int > states_t;
                        for ( int k = 0; k < n_haps; k++ ) {
                                states_t.push_back( k );
                        }
                        prune_states.push_back( states_t );
                        continue;
                }
	
		vector< vector< pair< int, int > > > hap1_store, hap2_store;
		vector< pair< int, int > > hap1_states, hap2_states;
		int ibd_distance;
		for ( int k = 0; k < n_haps; k++ ) {
			vector< pair< int, int > > _blks;
			shared_haplotypes( hap1, p_p[ k ], _blks, ibd_distance );
			hap1_store.push_back( _blks );
			hap1_states.push_back( pair< int, int >( k, ibd_distance ) );
		}
		for ( int k = 0; k < n_haps; k++ ) {
			vector< pair< int, int > > _blks;
			shared_haplotypes( hap2, p_p[ k ], _blks, ibd_distance );
			hap2_store.push_back( _blks );
			hap2_states.push_back( pair< int, int >( k, ibd_distance ) );
		}
		
		sort( hap1_states.begin(), hap1_states.end(), _reverse_pair_compare_() );
		sort( hap2_states.begin(), hap2_states.end(), _reverse_pair_compare_() );

		// conjunct haplotype 1.
		int hap_start = 0, hap_end = 0;
		vector< int > states;
		map< int, int > states_map;
		for ( size_t k = 0; k < hap1_states.size(); k++ ) {
			int id = hap1_states[ k ].first;
			if ( k == 0 ) {
				if ( hap1_store[ id ].size() == 1 ) {
					hap_start = hap1_store[ id ][ 0 ].first;
					hap_end   = hap1_store[ id ][ 0 ].second;
				} else {
					int l = 0;
					for ( size_t g = 0; g < hap1_store[ id ].size(); g++ ) {
						if ( hap1_states[ k ].second  == hap1_store[ id ][ g ].second - hap1_store[ id ][ g ].first ) {
							l = g;
							break;
						}
					}
					hap_start = hap1_store[ id ][ l ].first;
					hap_end   = hap1_store[ id ][ l ].second;
				}
				states.push_back( id );
				states_map[ id ] = 1;
				continue;
			}
			bool flg = false;
			for ( size_t g = 0; g < hap1_store[ id ].size(); g++ ) {
				int pseudo_start = hap1_store[ id ][ g ].first;
				int pseudo_end   = hap1_store[ id ][ g ].second;
				if ( pseudo_start < hap_start && pseudo_end > hap_start ) {
					hap_start = pseudo_start;
					flg = true;
				} else if ( pseudo_start < hap_end && pseudo_end > hap_end ) {
					hap_end   = pseudo_end;
					flg = true;
				} else {
					continue;
				}
			}
			if ( flg ) {
				states.push_back( id );
				states_map[ id ] = 1;
			}
			if ( hap_start <= 0 && hap_end >= n_site ) {
				break;
			}
		}
		// conjunct haplotype 2.
		hap_start = 0; hap_end = 0;
		vector< int > states2;
		for ( size_t k = 0; k < hap2_states.size(); k++ ) {
			int id = hap2_states[ k ].first;
			if ( k == 0 ) {
				if ( hap2_store[ id ].size() == 1 ) {
					hap_start = hap2_store[ id ][ 0 ].first;
					hap_end   = hap2_store[ id ][ 0 ].second;
				} else {
					int l = 0;
					for ( size_t g = 0; g < hap2_store[ id ].size(); g++ ) {
						if ( hap2_states[ k ].second  == hap2_store[ id ][ g ].second - hap2_store[ id ][ g ].first ) {
							l = g;
							break;
						}
					}
					hap_start = hap2_store[ id ][ l ].first;
					hap_end   = hap2_store[ id ][ l ].second;
				}
				states2.push_back( id );
				continue;
			}
			bool flg = false;
			for ( size_t g = 0; g < hap2_store[ id ].size(); g++ ) {
				int pseudo_start = hap2_store[ id ][ g ].first;
				int pseudo_end   = hap2_store[ id ][ g ].second;
				if ( pseudo_start < hap_start && pseudo_end > hap_start ) {
					hap_start = pseudo_start;
					flg = true;
				} else if ( pseudo_start < hap_end && pseudo_end > hap_end ) {
					hap_end   = pseudo_end;
					flg = true;
				} else {
					continue;
				}
			}
			if ( flg ) {
				states2.push_back( id );
			}
			if ( hap_start <= 0 && hap_end >= n_site ) {
				break;
			}
		}
		for ( size_t k = 0; k < states2.size(); k++ ) {
			if ( states_map.count( states2[ k ] ) ) {
				continue;
			}
			states.push_back( states2[ k ] );
		}
		prune_states.push_back( states );
	}
	
	if ( locate != _snp_number ) {
		cerr << "[_HMM_Algo::prune_states_find] fail to prune hidden states\n";
		abort();
	}
}

void _HMM_Algo::prune_states_find( Haplotype * haplotypes, int individuals_ID, int maximum )
{
	for ( size_t i = 0; i < prune_states.size(); i++ ) {
		prune_states[ i ].clear();
	}
	prune_states.clear();
	
	int locate = 0;
	int a, b, c, d, d1, d2, max_d1, max_d2;
	for ( int i = 0; i < _hg.n_kjds; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		ele hap1( n_site );
		memcpy( hap1.p, haplotypes[ individuals_ID * 2 + 0 ].haplotype + locate, n_site );
		ele hap2( n_site );
		memcpy( hap2.p, haplotypes[ individuals_ID * 2 + 1 ].haplotype + locate, n_site );
		locate += n_site;
		ele * p_p = _hg._kjd_array[ i ]->piece;
		vector< pair< int, int > > states;
		vector< pair< int, int > > IBDs;
		for ( int k = 0; k < n_haps; k++ ) {
			a = 0;
			b = 0;
			d1 = 0;
			d2 = 0;
			max_d1 = 0;
			max_d2 = 0;
			for ( int j = 0; j < n_site; j++ ) {
				d1++;
				d2++;
				if ( p_p[ k ].p[ j ] != hap1.p[ j ] ) {
					a++;
					d1 = 0;
				}
				if ( p_p[ k ].p[ j ] != hap2.p[ j ] ) {
					b++;
					d2 = 0;
				}
				if ( d1 > max_d1 ) {
					max_d1 = d1;
				}
				if ( d2 > max_d2 ) {
					max_d2 = d2;
				}
			}
			c = a > b ? b : a;
			d = max_d1 > max_d2 ? max_d1 : max_d2;
			states.push_back( pair< int, int >( k, c ) );
			IBDs.push_back( pair< int, int >( k, d ) );
		}
		vector< int > states_t;
		if ( states.size() <= maximum ) {
			for ( size_t g = 0; g < states.size(); g++ ) {
				states_t.push_back( states[ g ].first );
			}
			prune_states.push_back( states_t );
			continue;
		}
		
		map< int, int > state_counts;
		for ( size_t g = 0; g < states.size(); g++ ) {
			state_counts[ states[ g ].second ]++;
		}
		map< int, int >::iterator it1 = state_counts.begin();
                map< int, int >::iterator it2 = state_counts.end();
                int total_states = 0;
                int mismatch_id = 0;
                for ( ; it1 != it2; it1++ ) {
                        total_states += it1->second;
                        if ( total_states >= maximum * 2 / 3 ) {
				it1--;
                                mismatch_id = it1->first;
                                break;
                        }
                }

                sort( states.begin(), states.end(), _pair_compare_() );
		map< int, int > h_ids;
                for ( size_t k = 0; k < states.size(); k++ ) {
                        if ( states[ k ].second <= mismatch_id ) {
                                states_t.push_back( states[ k ].first );
				h_ids[ states[ k ].first ] = 1;
                        } else {
                                break;
                        }
                }
		sort( IBDs.begin(), IBDs.end(), _reverse_pair_compare_() );
		int n_index = 0;
		for ( size_t g = 0; g < IBDs.size(); g++ ) {
			if ( h_ids.count( IBDs[ g ].first ) ) {
				continue;
			} else {
				n_index++;
				if ( n_index >= maximum / 3 ) {
					break;
				}
				states_t.push_back( IBDs[ g ].first );
			}
		} 
                prune_states.push_back( states_t );
	}
	
	if ( locate != _snp_number ) {
		cerr << "[_HMM_Algo::prune_states_find] fail to prune hidden states\n";
		abort();
	}
}

void _HMM_Algo::prune_states_find( Haplotype * haplotypes, int individuals_ID, map< int, int > & freqs, int minimum, int maximum )
{
	for ( size_t i = 0; i < prune_states.size(); i++ ) {
		prune_states[ i ].clear();
	}
	prune_states.clear();

	vector< uint8_t > homs;
	vector< int > ids;
	uint8_t geno;
	int locate = 0;
	for ( int i = 0; i < _hg.n_kjds; i++ ) {
		int n_haps = _hg._kjd_array[ i ]->size;
		int n_site = _hg._kjd_array[ i ]->length;
		homs.clear();
		ids.clear();
		
		for ( int k = 0; k < n_site; k++ ) {
			if ( freqs.count( locate + k ) ) {
				continue;
			}
			geno = haplotypes[ individuals_ID * 2 + 0 ].haplotype[ locate + k ] + \
			       haplotypes[ individuals_ID * 2 + 1 ].haplotype[ locate + k ];
			if ( geno == 1 ) {
				continue;
			} else {
				homs.push_back( geno / 2 );
				ids.push_back( k );
			}
		}
		locate += n_site;
		
		vector< pair< int, int > > states;
		ele * p_p = _hg._kjd_array[ i ]->piece;
		for ( int k = 0; k < n_haps; k++ ) {
			uint8_t * p_hap = p_p[ k ].p;
			int c = 0;
			for ( size_t l = 0; l < ids.size(); l++ ) {
				if ( homs[ l ] != p_hap[ ids[ l ] ] ) {
					c++;
				}
			}
			states.push_back( pair< int, int >( k, c ) );
		}
		
		vector< int > states_t;
		if ( states.size() <= minimum ) {
			for ( size_t k = 0; k < states.size(); k++ ) {
				states_t.push_back( states[ k ].first );
			}
			prune_states.push_back( states_t );
			continue;
		}
		map< int, int > state_counts;
		for ( size_t k = 0; k < states.size(); k++ ) {
			state_counts[ states[ k ].second ]++;
		}
		map< int, int >::iterator it1 = state_counts.begin();
		map< int, int >::iterator it2 = state_counts.end();
		int total_states = 0;
		int mismatch_id = 0;
		for ( ; it1 != it2; it1++ ) {
			total_states += it1->second;
			if ( total_states >= maximum ) {
				it1--;
				mismatch_id = it1->first;
				break;
			}
		}

		states_t.clear();

		sort( states.begin(), states.end(), _pair_compare_() );
		for ( size_t k = 0; k < states.size(); k++ ) {
			if ( states[ k ].second <= mismatch_id ) {
				states_t.push_back( states[ k ].first );
			} else {
				break;
			}
		}
		prune_states.push_back( states_t );
		if ( states_t.size() > maximum ) {
			cerr << "states_t.size() = " << states_t.size() << endl;
			cerr << "total_states = " << total_states << endl;
			it1 = state_counts.begin();
			it2 = state_counts.end();
			cerr << "ids.size() = " << ids.size() << endl;
			total_states = 0;
			for ( ; it1 != it2; it1++ ) {
				total_states += it1->second;
				cerr << "total_states = " << total_states << "\tit1->second = " << it1->second << "\tit1->first = " << it1->first << endl;
			}
			abort();
		}
	}
	if ( locate != _snp_number ) {
		cerr << "[_HMM_Algo::prune_states_find] fail to prune hidden states.\n";
		abort();
	}
}

void _HMM_Algo::prune_alloc_leftmatrix()
{
	int locate = 0;
	_LeftMatrix = new float * [ _snp_number ];
	for ( size_t i = 0; i < prune_states.size(); i++ ) {
		int n_site = _hg._kjd_array[ i ]->length;
		int n_haps = prune_states[ i ].size();
		for ( int j = 0; j < n_site; j++ ) {
			_LeftMatrix[ locate + j ] = new float [ n_haps * n_haps ];
		}
		locate += n_site;
	}
	if ( locate != _snp_number ) {
		cerr << "[_HMM_Algo::prune_alloc_leftmatrix] fail to alloc memory for individual-specific HMM.\n";
		abort();
	}
}

void _HMM_Algo::prune_free_leftmatrix()
{
	for ( int i = 0; i < _snp_number; i++ ) {
		delete [] _LeftMatrix[ i ];
	}
	delete [] _LeftMatrix;
}

void _HMM_Algo::prune_prior()
{
	float * _tm = _LeftMatrix[ 0 ];
	int n_hidden_states = prune_states[ 0 ].size();
	uint32_t * p_w = _hg._kjd_array[ 0 ]->weigth;
	
	for ( int i = 0; i < n_hidden_states; i++ ) {
		float fi = static_cast< float >( double(p_w[ prune_states[ 0 ][ i ] ]) / _haplotype_number );
		for ( int j = i; j < n_hidden_states; j++ ) {
			_tm[i * n_hidden_states + j] = fi * static_cast< float >( double(p_w[ prune_states[ 0 ][ j ] ]) / _haplotype_number );
		}
	}
}

void _HMM_Algo::prune_match_data( int marker, int individuals_ID, _GenoHap & genohap )
{
	int cluster_id = _block_IDs[ marker ];
	int n_haps = prune_states[ cluster_id ].size();
	int offset = _block_offsets[ marker ];
	float * _tm = _LeftMatrix[ marker ];
	ele * p_p = _hg._kjd_array[ cluster_id ]->piece;

	double conditional_probs[3];
	for (int i = 0; i < 3; i++)
		conditional_probs[i] = _errors[ marker ][ i ][ 0 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 0 ] + \
				       _errors[ marker ][ i ][ 1 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 1 ] + \
				       _errors[ marker ][ i ][ 2 ] * genohap.glf[ individuals_ID * _indiv_ost + (marker + _start) * 3 + 2 ];
	double factors[ 2 ];
	for ( int i = 0; i < n_haps; i++ ) {
		uint8_t * p_i = p_p[ prune_states[ cluster_id ][ i ] ].p;
		factors[ 0 ] = conditional_probs[ p_i[ offset ] ];
		factors[ 1 ] = conditional_probs[ p_i[ offset ] + 1 ];
		for ( int j = i; j < n_haps; j++ ) {
			uint8_t * p_j = p_p[ prune_states[ cluster_id ][ j ] ].p;
			_tm[i * n_haps + j] *= static_cast< float >( factors[ p_j[ offset ] ] );
		}
	}
}

void _HMM_Algo::prune_prepare_pointors( int prev_marker, int next_marker )
{
	int prev_cluster = _block_IDs[ prev_marker ];
	int next_cluster = _block_IDs[ next_marker ];
	int prev_size    = prune_states[ prev_cluster ].size();
	int next_size    = prune_states[ next_cluster ].size();
	
	ele * prev_p = _hg._kjd_array[ prev_cluster ]->piece;
	prev_piece = new uint8_t * [ prev_size ];
	for ( int i = 0; i < prev_size; i++ ) {
		prev_piece[ i ] = prev_p[ prune_states[ prev_cluster ][ i ] ].p;
	}
	prev_weigth = _hg._kjd_array[ prev_cluster ]->weigth;

	ele * next_p = _hg._kjd_array[ next_cluster ]->piece;
	next_piece = new uint8_t * [ next_size ];
	for ( int i = 0; i < next_size; i++ ) {
		next_piece[ i ] = next_p[ prune_states[ next_cluster ][ i ] ].p;
	}
	next_weigth = _hg._kjd_array[ next_cluster ]->weigth;
}

void _HMM_Algo::prune_destory_pointors()
{
	delete [] prev_piece;
	delete [] next_piece;
	prev_piece = NULL;
	next_piece = NULL;
	prev_weigth = NULL;
	next_weigth = NULL;
}

void _HMM_Algo::prune_transpose( int prev_marker, int next_marker, float * theta_array, _GenoHap & genohap, int individuals_ID )
{
	clear_matrix();
	float * source = _LeftMatrix[ prev_marker ];
	float * dest   = _LeftMatrix[ next_marker ];
	float   theta  = theta_array[ prev_marker ];

	if ( _block_flags[ prev_marker ] ) {
		int prev_cluster = _block_IDs[ prev_marker ];
		int next_cluster = _block_IDs[ next_marker ];
		int prev_size    = prune_states[ prev_cluster ].size();
		int next_size    = prune_states[ next_cluster ].size();
		
		prune_prepare_pointors( prev_marker, next_marker );
		double sum = 0;
		for ( int i = 0; i < prev_size; i++ ) {
			sum += source[ i * prev_size + i ];
			_marginal_i[ i ] += source[ i * prev_size + i ];
			for ( int j = i + 1; j < prev_size; j++ ) {
				sum += source[ i * prev_size + j ] * 2;
				_marginal_i[ i ] += source[ i * prev_size + j ];
				_marginal_i[ j ] += source[ j * prev_size + i ];
			}
		}
		for (int i = 0; i < prev_size; i++ ) {
			_marginal_j[ i ] = _marginal_i[ i ];
		}
		
		vector< vector< pair< int, double > > > freqs;
		vector< pair< int, double > > freq;

		double no_change  = (1. - theta) * (1. - theta);
		double one_change = (1. - theta) * theta;
		double two_change = theta * theta * sum;

		map< int, int > prev_prune_haps;
		for ( size_t prune_haps_index = 0; prune_haps_index < prune_states[ prev_cluster ].size(); prune_haps_index++ ) {
			prev_prune_haps[ prune_states[ prev_cluster ][ prune_haps_index ] ] = prune_haps_index;
		}
		for ( int k = 0; k < next_size; k++ ) {
			freq.clear();
			map< int, int >::iterator iter1 = _hg._kjds[ next_cluster ].prev_edges[ prune_states[ next_cluster ][ k ] ].begin();
			map< int, int >::iterator iter2 = _hg._kjds[ next_cluster ].prev_edges[ prune_states[ next_cluster ][ k ] ].end();
			for ( ; iter1 != iter2; iter1++ ) {
				if ( ! prev_prune_haps.count( iter1->first ) ) {
					continue;
				}
				double f = double(iter1->second) / double(prev_weigth[ iter1->first ]);
				freq.push_back( pair< int, double >( prev_prune_haps[ iter1->first ], f ) );
			}
			freqs.push_back(freq);
		}

		for ( int k = 0; k < next_size; k++ ) {
			for ( int l = k; l < next_size; l++ ) {
				double fl = double(next_weigth[ prune_states[ next_cluster ][ l ] ]) / _haplotype_number;
				double fk = double(next_weigth[ prune_states[ next_cluster ][ k ] ]) / _haplotype_number;
				
				vector< pair< int, double > > & freq_k = freqs[k];
				vector< pair< int, double > > & freq_l = freqs[l];
				
				double prob = 0;
				double prob_t = 0;
				// no recombination
				for ( size_t n_freq_k = 0; n_freq_k < freq_k.size(); n_freq_k++ ) {
					for ( size_t n_freq_l = 0; n_freq_l < freq_l.size(); n_freq_l++ ) {
						int id_i = freq_k[ n_freq_k ].first;
						int id_j = freq_l[ n_freq_l ].first;
						prob_t += source[ id_i * prev_size + id_j ] * freq_k[ n_freq_k ].second * \
							freq_l[ n_freq_l ].second;
					}
				}
				prob_t *= no_change;
				prob += prob_t;
				prob_t = 0;
	
				// j change
				for ( size_t n_freq_k = 0; n_freq_k < freq_k.size(); n_freq_k++ ) {
					int id_k = freq_k[ n_freq_k ].first;
					prob_t += _marginal_i[ id_k ] * freq_k[ n_freq_k ].second;
				}
				prob_t *= fl * one_change;
				prob += prob_t;
				prob_t = 0;
				// k change
				for ( size_t n_freq_l = 0; n_freq_l < freq_l.size(); n_freq_l++ ) {
					int id_l = freq_l[ n_freq_l ].first;
					prob_t += _marginal_j[ id_l ] * freq_l[ n_freq_l ].second;
				}
				prob_t *= fk * one_change;
				prob += prob_t;
				prob_t = 0;
				// k & j both changes
				prob_t = two_change * fk * fl;
				prob += prob_t;
				_LeftMatrix[ next_marker ][ k * next_size + l ] = static_cast< float >( prob );
			}
		}
		prune_destory_pointors();
	} else {
		int prev_cluster = _block_IDs[ prev_marker ];
		int n_haps = prune_states[ prev_cluster ].size();
		
		float * gl = &genohap.glf[ individuals_ID * _indiv_ost + (_start + next_marker) * 3 ];
	
		double sum = 0.;
		for ( int i = 0; i < n_haps; i++ ) {
			sum += source[ i * n_haps + i ];
			_marginal_i[ i ] += source[ i * n_haps + i ];
			for ( int j = i + 1; j < n_haps; j++ ) {
				sum += source[ i * n_haps + j ] * 2;
				_marginal_i[ i ] += source[ i * n_haps + j ];
				_marginal_i[ j ] += source[ j * n_haps + i ];
			}
		}
		for (int i = 0; i < n_haps; i++ ) {
			_marginal_j[ i ] = _marginal_i[ i ];
		}
		
		double no_change  = (1. - theta) * (1. - theta);
		double one_change = (1. - theta) * theta / _haplotype_number;
		double two_change = sum * theta * theta  / (_haplotype_number * _haplotype_number);
		
		uint32_t * p_w_prev = _hg._kjd_array[ prev_cluster ]->weigth;

		for ( int i = 0; i < n_haps; i++ ) {
			for ( int j = i; j < n_haps; j++ ) {
				dest[i * n_haps + j] = source[i * n_haps + j] * no_change + \
					_marginal_i[ i ] * one_change * p_w_prev[ prune_states[ prev_cluster ][ j ] ] + \
					_marginal_j[ j ] * one_change * p_w_prev[ prune_states[ prev_cluster ][ i ] ] + \
					two_change * p_w_prev[ prune_states[ prev_cluster ][ i ] ] * p_w_prev[ prune_states[ prev_cluster ][ j ] ];
			}
		}
	}
}

void _HMM_Algo::prune_scale( int marker )
{
	int cluster_id = _block_IDs[ marker ];
	int n_haps = prune_states[ cluster_id ].size();
	
	float max_prob = 0;
	for ( int i = 0; i < n_haps; i++ ) {
		for ( int j = i; j < n_haps; j++ ) {
			if ( _LeftMatrix[ marker ][ i * n_haps + j ] > max_prob ) {
				max_prob = _LeftMatrix[ marker ][ i * n_haps + j ];
			}
		}
	}

	double multiplier = 1. / max_prob;
	for ( int i = 0; i < n_haps; i++ ) {
		for ( int j = i; j < n_haps; j++ ) {
			_LeftMatrix[ marker ][ i * n_haps + j ] *= multiplier;
			_LeftMatrix[ marker ][ j * n_haps + i ] = _LeftMatrix[ marker ][ i * n_haps + j ];
		}
	}
}

void _HMM_Algo::prune_process( float * theta_array, _GenoHap & genohap, int individuals_ID )
{
	prune_prior();
	prune_match_data( 0, individuals_ID, genohap );
	prune_scale( 0 );

	for ( int marker = 1; marker < _snp_number; marker++ ) {
		int prev_marker = marker - 1;
		int next_marker = marker;
		prune_transpose( prev_marker, next_marker, theta_array, genohap, individuals_ID );
		prune_match_data( next_marker, individuals_ID, genohap );
		prune_scale( next_marker );
	}
}

void _HMM_Algo::prune_update_haplotype( Haplotype * haplotypes, int individuals_ID, _GenoHap & genohap, float * theta_array, \
			Error * error_model, variate_generator< kreutzer1986, uniform_real<> > & uniform, Theta * crossover, float * mutation_rate  )
{
	float * alpha = _LeftMatrix[ _snp_number - 1 ];
	double sum = 0;
	int cluster_id = _block_IDs[ _snp_number - 1 ];
	int n_haps = prune_states[ cluster_id ].size();
	for ( int i = 0; i < n_haps; i++ ) {
		sum += alpha[ i * n_haps + i ];
		for ( int j = i + 1; j < n_haps; j++ ) {
			sum += alpha[ i * n_haps + j ] * 2;
		}
	}
	double choice = uniform() * sum;
	int first = -1;
	int second = -1;
	double tmp = 0;
	for ( int i = 0; i < n_haps; i++ ) {
		 tmp -= alpha[ i * n_haps + i ];
		for ( int j = i; j < n_haps; j++ ) {
			tmp += alpha[ i * n_haps + j ] * 2;
			if ( tmp > choice ) {
				 if (uniform() > 0.5) {
					first = i;
					second = j;
				} else {
					first = j;
					second = i;
				}
				goto HP;
			}
		}
	} // back samples and randomly select the haplotype pair

HP:
	for ( int marker = _snp_number - 2; marker >= 0; marker-- ) {
		prune_update_marker( first, second, individuals_ID, marker + 1, error_model, genohap, uniform, haplotypes, mutation_rate );
		double theta = theta_array[ marker ];
		if ( _block_flags[ marker ] ) { // inter edge
			alpha = _LeftMatrix[ marker ];
			int prev_marker = marker;
			int next_marker = marker + 1;
			int prev_cluster = _block_IDs[ prev_marker ];
			int next_cluster = _block_IDs[ next_marker ];
			uint32_t * p_w_prev = _hg._kjd_array[ prev_cluster ]->weigth;
			uint32_t * p_w_next = _hg._kjd_array[ next_cluster ]->weigth;

			n_haps = prune_states[ prev_cluster ].size();
			int p_haps = prune_states[ next_cluster ].size(); // poster haplotypes
			vector< int > & prev_states = prune_states[ prev_cluster ];
			vector< int > & next_states = prune_states[ next_cluster ];
			double sum11 = 0;
			vector<double> vec11(n_haps * n_haps, 0.);
			vector<double> part1, part2;
			
			for ( int k = 0; k < n_haps; k++ ) {
				part1.push_back((1. - theta) * 1. / double(p_w_prev[ prev_states[ k ] ]));
			}
			for ( int k = 0; k < p_haps; k++ ) {
				part2.push_back(theta * double(p_w_next[ next_states[ k ] ]) / _haplotype_number);
			}

			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					double C_Km_Km1_A = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ].count( prev_states[ k ] ) ) {
						C_Km_Km1_A = double( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ][ prev_states[ k ] ] );
					 }
					double C_Km_Km1_B = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ].count( prev_states[ l ] ) ) {
						C_Km_Km1_B = double( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ][ prev_states[ l ] ] );
					 }
					double C_Km_Km1_C = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ].count( prev_states[ l ] ) ) {
						C_Km_Km1_C = double( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ][ prev_states[ l ] ] );
					 }
					double C_Km_Km1_D = 0;
					 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ].count( prev_states[ k ] ) ) {
						C_Km_Km1_D = double( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ][ prev_states[ k ] ] );
					 }
					/*
					double tmp1 = (1. - theta) * C_Km_Km1_A / double(p_w_prev[ k ]) + theta * double(p_w_next[ first ]) / _haplotype_number;
					tmp1 *= (1. - theta) * C_Km_Km1_B / double(p_w_prev[ l ]) + theta * double(p_w_next[ second ]) / _haplotype_number;
					double tmp2 = (1. - theta) * C_Km_Km1_C / double(p_w_prev[ l ]) + theta * double(p_w_next[ first ]) / _haplotype_number;
					tmp2 *= (1. - theta) * C_Km_Km1_D / double(p_w_prev[ k ]) + theta * double(p_w_next[ second ]) / _haplotype_number;
					*/
					double tmp1 = ( part1[k] * C_Km_Km1_A + part2[first] ) * ( part1[l] * C_Km_Km1_B + part2[second] );
					double tmp2 = ( part1[l] * C_Km_Km1_C + part2[first] ) * ( part1[k] * C_Km_Km1_D + part2[second] );
					
					tmp1 *= alpha[ k * n_haps + l ];
					tmp2 *= alpha[ l * n_haps + k ];
					sum11 += (tmp1 + tmp2);
					vec11[k * n_haps + l] = tmp1;
					vec11[l * n_haps + k] = tmp2;
				}
				sum11 -= vec11[k * n_haps + k];
			}
			choice = uniform() * sum11;
			tmp = 0;
			int first0 = -1;
			int second0 = -1;
			for ( size_t g = 0; g < vec11.size(); g++ ) {
				tmp += vec11[ g ];
				if ( tmp > choice ) {
					first0 = g / n_haps;
					second0 = g % n_haps;
					break;
				}
			}

			/// select weather recombination events occur
			double C_Km_Km1_A = 0;
			 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ].count( prev_states[ first0 ] ) ) {
				C_Km_Km1_A = double ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ first ] ][ prev_states[ first0 ] ] );
			 }
			double C_Km_Km1_B = 0;
			 if ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ].count( prev_states[ second0 ] ) ) {
				C_Km_Km1_B = double ( _hg._kjds[ next_cluster ].prev_edges[ next_states[ second ] ][ prev_states[ second0 ] ] );
			 }
			double prob_A_one = (1. - theta) * C_Km_Km1_A / double(p_w_prev[ prev_states[ first0 ] ]);
			double prob_A_two = theta * double(p_w_next[ next_states[ first ] ]) / _haplotype_number;
			double prob_B_one = (1. - theta) * C_Km_Km1_B / double(p_w_prev[ prev_states[ second0 ] ]);
			double prob_B_two = theta * double(p_w_next[ next_states[ second ] ]) / _haplotype_number;
			
			vector< double > probs;
			probs.push_back( prob_A_one * prob_B_one );
			probs.push_back( prob_A_one * prob_B_two + prob_A_two * prob_B_one );
			probs.push_back( prob_A_two * prob_B_two );

			sum = probs[ 0 ] + probs[ 1 ] + probs[ 2 ];

			crossover[ marker ].cross += (probs[ 1 ] + probs[ 2 ] * 2) / sum;
		/*
			choice = uniform() * sum;
			tmp = 0;
			int change = 0;
			for ( int k = 0; k < 3; k++ ) {
				tmp += probs[ k ];
				if ( tmp > choice ) {
					change = k;
					break;
				}
			}
		*/
			crossover[ marker ].total += 2;
		//	crossover[ marker ].cross += change;

			first = first0;
			second = second0;
		} else {
			alpha = _LeftMatrix[ marker ];
			cluster_id = _block_IDs[ marker ];
			n_haps = prune_states[ cluster_id ].size();
			uint32_t * p_w = _hg._kjd_array[ cluster_id ]->weigth;
			vector< int > & states = prune_states[ cluster_id ];

			double sum00( 0 ), sum10( 0 ), sum01( 0 ), sum11( 0 );
			vector<double> vec10, vec01, vec11;
			vector< int >  ids10, ids01, ids11;
			double both_change = theta * theta / (_haplotype_number * _haplotype_number);
			double one_change = (1. - theta) * theta / _haplotype_number;
			double no_change = (1. - theta) * (1. - theta);
			for ( int k = 0; k < n_haps; k++ ) {
				for ( int l = k; l < n_haps; l++ ) {
					tmp = alpha[ k * n_haps + l ] * double(p_w[ states[ first ] ] * p_w[ states[ second ] ]) * both_change;
					tmp *= (1 + (k != l));
					vec11.push_back( tmp );
					ids11.push_back( k * n_haps + l );
					sum11 += tmp;
					if ( l == second ) {
						tmp = alpha[ k * n_haps + l ] * double(p_w[ states[ first ] ]) * one_change;
						vec10.push_back( tmp );
						sum10 += tmp;
						ids10.push_back( k * n_haps + l );
					}
					if ( k == second && l != k ) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ states[ first ] ]) * one_change;
						vec10.push_back( tmp );
						sum10 += tmp;
						ids10.push_back( l * n_haps + k );
					}
					if ( k == first ) {
						tmp = alpha[ k * n_haps + l ] * double(p_w[ states[ second ] ]) * one_change;
						vec01.push_back( tmp );
						sum01 += tmp;
						ids01.push_back( k * n_haps + l );
					}
					if ( l == first && l != k ) {
						tmp = alpha[ l * n_haps + k ] * double(p_w[ states[ second ] ]) * one_change;
						vec01.push_back( tmp );
						sum01 += tmp;
						ids01.push_back( l * n_haps + k );
					}
					if ( k == first && l == second) {
						sum00 += (alpha[ k * n_haps + l ] * no_change);
					}
					if ( l == first && k == second && l != k) {
						sum00 += (alpha[ l * n_haps + k ] * no_change);
					}
				}
			}

			sum = sum00 + sum01 + sum10 + sum11;

			int pesudo_first = -1;
			int pesudo_second = -1;
			crossover[ marker ].total += 2;
			crossover[ marker ].cross += (sum01 + sum10 + sum11 * 2) / sum;

			choice = uniform() * sum;
			if ( choice < sum00 ) {
				continue;
			}
			choice -= sum00;
			if ( choice < sum01 ) {
				tmp = 0;
				for ( size_t k = 0; k < vec01.size(); k++ ) {
					tmp += vec01[ k ];
					if ( tmp > choice ) {
						pesudo_first = ids01[ k ] / n_haps;
						pesudo_second = ids01[ k ] % n_haps;
						break;
					}
				}
			//	crossover[ marker ].cross++;
				second = pesudo_second;
				continue;
			}
			choice -= sum01;
			if ( choice < sum10 ) {
				tmp = 0;
				for ( size_t k = 0; k < vec10.size(); k++ ) {
					tmp += vec10[ k ];
					if ( tmp > choice ) {
						pesudo_first = ids10[ k ] / n_haps;
						pesudo_second = ids10[ k ] % n_haps;
						break;
					}
				}
			//	crossover[ marker ].cross++;
				first = pesudo_first;
				continue;
			}
			choice -= sum10;
			if( choice < sum11 ){
				tmp = 0;
				for ( size_t k = 0; k < vec11.size(); k++ ) {
					tmp += vec11[ k ];
					if ( tmp > choice ) {
						if (uniform() > 0.5) {
							pesudo_first = ids11[ k ] / n_haps;
							pesudo_second = ids11[ k ] % n_haps;
						} else {
							pesudo_first = ids11[ k ] % n_haps;
							pesudo_second = ids11[ k ] / n_haps;
						}
						break;
					}
				}
			//	crossover[ marker ].cross += 2;
				first = pesudo_first;
				second = pesudo_second;
			}
		}
	}
	prune_update_marker( first, second, individuals_ID, 0, error_model, genohap, uniform, haplotypes, mutation_rate );
}

void _HMM_Algo::prune_update_marker( int first, int second, int individuals_ID, int marker, Error * error_model, _GenoHap & genohap, \
				variate_generator< kreutzer1986, uniform_real<> > & uniform, Haplotype * haplotypes, float * mutation_rate )
{
	int offset = _block_offsets[ marker ];
	int cluster_id = _block_IDs[ marker ];
	ele * p_p = _hg._kjd_array[ cluster_id ]->piece;
	uint8_t * p_first  = p_p[ prune_states[ cluster_id ][ first ] ].p;
	uint8_t * p_second = p_p[ prune_states[ cluster_id ][ second ] ].p;
	int base1 = p_first[ offset ];
	int base2 = p_second[ offset ];

	float * glikes = &genohap.glf[individuals_ID * _indiv_ost + (marker + _start) * 3 ];
	double posterior_11 = _errors[ marker ][ base1 + base2 ][ 0 ] * glikes[ 0 ];
	double posterior_12 = _errors[ marker ][ base1 + base2 ][ 1 ] * glikes[ 1 ];
	double posterior_22 = _errors[ marker ][ base1 + base2 ][ 2 ] * glikes[ 2 ];
	double sum = posterior_11 + posterior_12 + posterior_22;
	posterior_11 /= sum;
	posterior_22 /= sum;
	posterior_12 /= sum;
	
	double r = uniform();
	
	if ( r < posterior_11 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 0;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 0;
	} else if ( r < posterior_11 + posterior_22 ) {
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = 1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = 1;
	} else if ( base1 != base2 ) {
		float rate = mutation_rate[ marker ];
		if ( uniform() < rate * rate / ( rate * rate + (1 - rate) * (1 - rate) ) ) {
			base1 = !base1;
			base2 = !base2;
		}
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = base1;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = base2;
	} else {
		bool bit = (uniform() > 0.5) ? 1 : 0;
		haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ] = bit;
		haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ] = bit ^ 1;
	}

	int imputed1 = haplotypes[ individuals_ID * 2 + 0 ].haplotype[ marker ];
	int imputed2 = haplotypes[ individuals_ID * 2 + 1 ].haplotype[ marker ];
	
	int differences = abs(base1 - imputed1) + abs(base2 - imputed2);
	
	error_model[ marker ].mismatch += differences;
	error_model[ marker ].match += 2 - differences;
}
