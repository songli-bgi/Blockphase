#include "Run_Wizzard.h"

// all regions are [*,*).
void _Run_Wizzard::split( _Program_Variables & pv )
{
	assert( pv.piece > pv.overlap );
	assert( pv.overlap > 0 );
	assert( pv.piece > 0 );

	int start_id = 0;
	if ( pv.begin > 0 ) {
		start_id = get_index( pv.begin );
	}
	int end_id = snp_number;
	if ( pv.end > 0 ) {
		end_id = get_index( pv.end );
	}
	assert( end_id > start_id );

	int piece = pv.piece;
	int overlap = pv.overlap;

	int locate = start_id;
	while ( locate < end_id ) {
		if ( locate + piece < end_id ) {
			Run_Pieces.push_back( pair<int, int>( locate, locate + piece ) );
		} else {
			Run_Pieces.push_back( pair<int, int>( locate, end_id ) );
			break;
		}
		locate += piece - overlap;
	}
}

void _Run_Wizzard::show_variables( _Program_Variables & pv, FILE * fp )
{
	char * p = strrchr( pv.in.c_str(), '/' );
	if ( p == NULL ) {
		fprintf( fp, "Reading SNPs in [%s] in Genotypes Likelihoods\n", pv.in.c_str() );
	} else {
		fprintf( fp, "Reading SNPs in [%s] in Genotypes Likelihoods\n", p + 1 );
	}
	fprintf( fp, "  * %d SNPs included\n", snp_number );
	fprintf( fp, "  * %d Individuals included\n\n", individuals );

	fprintf( fp, "Others :\n" );
	
	int begin_pos = 0;
	int end_pos   = 0;

	if ( pv.begin > 0 ) {
		begin_pos = pv.begin;
	} else {
		begin_pos = snpids[ 0 ];
	}

	if ( pv.end > 0 ) {
		end_pos = pv.end;
	} else {
		end_pos = snpids[ snp_number -1 ];
	}
	
	fprintf( fp, "  * Begin marker : %d\n", begin_pos );
	fprintf( fp, "  * End   marker : %d\n", end_pos );
	fprintf( fp, "  * Input file   : %s\n", pv.in.c_str() );
	if ( pv.out.empty() ) {
		fprintf( fp, "  * Output file  : phase_haplotypes\n" );
		fprintf( fp, "  * Log file     : phase_haplotypes.log\n" );
	} else {
		fprintf( fp, "  * Output file  : %s\n", pv.out.c_str() );
		fprintf( fp, "  * Log file     : %s.log\n", pv.out.c_str() );
	}

	if ( pv.reference.empty() ) {
		fprintf( fp, "  * Reference    : No Reference panel used\n" );
	}
}

int _Run_Wizzard::get_index( int pos )
{
	assert( pos > 0 );
	int ret = 0;
	for ( int i = 0; i < snp_number; i++ ) {
		if ( snpids[ i ] > pos ) {
			ret = i - 1;
			break;
		}
	}
	return ret;
}
