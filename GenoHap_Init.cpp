#include "GenoHap_Init.h"

#include <stdlib.h>
#include <cmath>

static int string_split(string & str, char * delim, vector<string> & vec_str)
{
	vec_str.clear();
	size_t len = str.length();
	if (len == 0) return 0;
	if (delim == NULL) {
		vec_str.push_back(str);
		return 1;
	}
	size_t start = 0, found;
	string line;
	while (start < len) {
		found = str.find_first_of(delim, start);
		if (found == string::npos) found = len;
		if (found > start) {
			line = str.substr(start, found - start);
			vec_str.push_back(line);
		}
		start = found + 1;
	}
	return vec_str.size();
}

void _GenoHap::selfinit_glf(string & infile, _Run_Wizzard & RW )
{
	freader ifs;
	int ret = ifs.open( infile.c_str() );
	if ( ret == false ) {
		cerr << "[_GenoHap::selfinit_glf] fail to open " << infile << "\n";
		abort();
	}

	const char * s = NULL;
	string line;

	if ( NULL != (s = ifs.getline( '\n' )) ) {
		line = s;
	} else {
		cerr << "[_GenoHap::selfinit_glf] fail to read " << infile << "\n";
		abort();
	}
	vector<string> sample_ID_vector;
	int indiv_n = (string_split(line, " \t", sample_ID_vector) - 3) / 3;
	for ( int i = 3; i < sample_ID_vector.size(); i = i + 3 ) {
		RW.samples.push_back( sample_ID_vector[ i ] );
	}
	int snp_n = 0;
	while (NULL != (s = ifs.getline( '\n' ))) ++snp_n;
	ifs.close();
	glf = new float[ indiv_n * snp_n * 3 ];
	int indiv_ost = snp_n * 3;
	RW.ref_Alt_bases = new char * [ snp_n ];
	for ( int i = 0; i < snp_n; i++ ) {
		RW.ref_Alt_bases[ i ] = new char [ 2 ];
	}
	ifs.open(infile.c_str());
	s = ifs.getline( '\n' );
	int snp_index = 0;
	vector<string> gl_line;
	while (NULL != (s = ifs.getline( '\n' ))) {
		line = s;
		string_split(line, " \t", gl_line);
		RW.snpids.push_back( atoi(gl_line[0].c_str()) );
		RW.ref_Alt_bases[ snp_index ][ 0 ] = gl_line[ 1 ].at( 0 );
		RW.ref_Alt_bases[ snp_index ][ 1 ] = gl_line[ 2 ].at( 0 );
		for (int i = 0; i != indiv_n; ++i) {
			glf[ i * indiv_ost + snp_index * 3 + 0 ] = atof( gl_line[ 3 + i * 3 + 0 ].c_str() );
			glf[ i * indiv_ost + snp_index * 3 + 1 ] = atof( gl_line[ 3 + i * 3 + 1 ].c_str() );
			glf[ i * indiv_ost + snp_index * 3 + 2 ] = atof( gl_line[ 3 + i * 3 + 2 ].c_str() );
		}
		snp_index++;
	}
	ifs.close();
	individuals = indiv_n;
	snp_number = snp_n;
}

void _GenoHap::random_createGeno( int start, int stop )
{
	int pseudo_snp_number = stop - start;
	for ( int i = 0; i < genotypes.size(); i++ ) {
		genotypes[ i ].clear();
	}
	genotypes.clear();
	
	double * freqs = new double[ pseudo_snp_number ];
	memset( freqs, 0 , sizeof( double ) * pseudo_snp_number );
	int indiv_ost = snp_number * 3;

	for ( int i = 0; i < pseudo_snp_number; i++ ) {
		double pu, diff;
		double p = 0.05;
		double accu = 1e-6;
                int n_iter = 100;
                double p11, p12, p22;
                double w0, w1, w2;
                for (int it = 0; it != n_iter; ++it) {
                        p11 = (1.0 - p) * (1.0 - p);
                        p12 = 2 * p * (1.0 - p);
                        p22 = p * p;
                        double sum = .0;
                        for (int j = 0; j != individuals; ++j) {
                                w0 = glf[ j * indiv_ost + (i + start) * 3 + 0 ] * p11;
                                w1 = glf[ j * indiv_ost + (i + start) * 3 + 1 ] * p12;
                                w2 = glf[ j * indiv_ost + (i + start) * 3 + 2 ] * p22;
                                sum += (w1 + 2 * w2) / (w0 + w1 + w2);
                        }
                        pu = sum / (2 * individuals);
                        diff = pu - p > 0 ? pu - p : p - pu;
                        p = pu;
                        if (diff < accu) break;
                }
                freqs[ i ] = p;
	}

	genotypes.resize(individuals);
	for ( int i = 0; i != individuals; i++ ) {
		genotypes[ i ].resize(pseudo_snp_number);
	}

	for ( int j = 0; j < pseudo_snp_number; j++ ) {
		double prior00 = (1. - freqs[ j ]) * (1. - freqs[ j ]);
		double prior01 = 2 * freqs[ j ] * (1. - freqs[ j ]);
		double prior11 = freqs[ j ] * freqs[ j ];
		for ( int i = 0; i < individuals; i++ ) {
			double poster00 = glf[ i * indiv_ost + (j + start) * 3 + 0 ] * prior00;
			double poster01 = glf[ i * indiv_ost + (j + start) * 3 + 1 ] * prior01;
			double poster11 = glf[ i * indiv_ost + (j + start) * 3 + 2 ] * prior11;
			double sum = poster00 + poster01 + poster11;
			poster00 /= sum;
			poster01 /= sum;
			double rand_num = (rand() + 0.) / (RAND_MAX + 0.);
			if ( rand_num < poster00 ) {
				genotypes[ i ][ j ] = 0;
			} else if ( rand_num < poster00 + poster01 ) {
				genotypes[ i ][ j ] = 1;
			} else {
				genotypes[ i ][ j ] = 2;
			}
		}
	}

	delete [] freqs;
}

void _GenoHap::em_createGeno( int start, int stop )
{
	int pseudo_snp_number = stop - start;
	for ( int i = 0; i < genotypes.size(); i++ ) {
		genotypes[ i ].clear();
	}
	genotypes.clear();
	genotypes.resize(individuals);
	for (int i = 0; i != individuals; ++i)
		genotypes[ i ].resize(pseudo_snp_number);
}

void _GenoHap::createHap( int _snp_number )
{
	hap_t _h( _snp_number );
	haplotypes = new hap_t [ individuals ];
	for ( int i = 0; i < individuals; i++ ) {
		haplotypes[ i ] = _h;
	}
}

void _GenoHap::randomHap( int _snp_number )
{
	for ( int i = 0; i < individuals; i++ ) {
		for ( int j = 0; j < _snp_number; j++ ) {
			if ( genotypes[ i ][ j ] == 1 ) {
				if ( (rand() + 0.) / (RAND_MAX + 0.) > 0.5 ) {
					haplotypes[ i ].first[ j ] = 0;
					haplotypes[ i ].second[ j ] = 1;
				} else {
					haplotypes[ i ].first[ j ] = 1;
					haplotypes[ i ].second[ j ] = 0;
				}
			} else if ( genotypes[ i ][ j ] == 0 ) {
				haplotypes[ i ].first[ j ] = 0;
				haplotypes[ i ].second[ j ] = 0;
			} else {
				haplotypes[ i ].first[ j ] = 1;
				haplotypes[ i ].second[ j ] = 1;
			}
		}
	}
}

void _GenoHap::free_hap()
{
	if (haplotypes != NULL) delete [] haplotypes;
	haplotypes = NULL;
}

void _GenoHap::destroyHap()
{
	if (haplotypes != NULL) delete [] haplotypes;
	haplotypes = NULL;
	if (glf != NULL) {
		delete [] glf;
		glf = NULL;
	}
	if ( whole_reference != NULL ) {
		delete [] whole_reference;
		whole_reference = NULL;
		whole_n_reference = 0;
	}
}

void _GenoHap::load_reference( const char * ref_file, vector< int > & snpids, char ** ref_Alt_bases )
{
	if ( ref_file == NULL ) {
		whole_reference = NULL;
		whole_n_reference = 0;
		return;
	}
	
	freader ifs;
	int ret = ifs.open( ref_file );
	if ( ret == false ) {
		cerr << "[_GenoHap::load_reference] fail to open " << ref_file << endl;
		abort();
	}
	const char * s = NULL;
	string line;
	if ( NULL != (s = ifs.getline( '\n' )) ) {
		line = s;
	} else {
		cerr << "[_GenoHap::load_reference] fail to read " << ref_file << endl;
		abort();
	}
	vector<string> ref_ID_vector;
	whole_n_reference = string_split(line, " \t", ref_ID_vector) - 1;
	int pseudo_snp_number = 0;
	while (NULL != (s = ifs.getline( '\n' ))) ++pseudo_snp_number;
	ifs.close();

	if ( pseudo_snp_number != snp_number ) {
		fprintf( stderr, "No. GLF markers should be equal to reference panel markers.\n" );
		abort();
	}

	Haplotype hap( pseudo_snp_number );
	whole_reference = new Haplotype [ whole_n_reference ];
	for ( int i = 0; i < whole_n_reference; i++ ) {
		whole_reference[ i ] = hap;
	}

	ifs.open( ref_file );
	s = ifs.getline( '\n' );

	int snp_index = 0;
	while (NULL != (s = ifs.getline( '\n' ))) {
		line = s;
		istringstream is_str( line );
		int pos;
		char base;
		uint8_t code;
		is_str >> pos;
		if ( snpids[ snp_index ] != pos ) {
			cerr << "[_GenoHap::load_reference] read wrong snp, while physical position don't match.\n";
			abort();
		}
		for ( int hap_id = 0; hap_id < whole_n_reference; hap_id++ ) {
			is_str >> base;
			if ( base == ref_Alt_bases[ snp_index ][ 0 ] ) {
				code = 0;
			} else if ( base == ref_Alt_bases[ snp_index ][ 1 ] ) {
				code = 1;
			} else {
				cerr << "[_GenoHap::load_reference] triallelic marker.\n";
				abort();
			}
			whole_reference[ hap_id ].haplotype[ snp_index ] = code;
		}
		snp_index++;
	}
	ifs.close();
}

void _GenoHap::show_variables( FILE * fp, _Program_Variables & pv )
{
	if ( whole_reference != NULL ) {
		fprintf( fp, "\nReference :\n" );
		fprintf( fp, "  * %d Reference SNPs included\n", snp_number );
		fprintf( fp, "  * %d Reference haplotypes included\n", whole_n_reference );
		fprintf( fp, "  * Reference File : %s\n", pv.reference.c_str() );
	}
}
