#include "Gibbsgeno_Store.h"

void _Gibbsgeno_Store::set( int marker, int flag, int cycle, int phase, int individuals_ID )
{
	int pos = marker / 8;
	int off = marker % 8;
	switch( off )
	{
		case 0:
			if ( flag )
				gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x01;
			else
				gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xFE;
			break;
		case 1:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x02;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xFD;
			break;
		case 2:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x04;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xFB;
			break;
		case 3:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x08;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xF7;
			break;
		case 4:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x10;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xEF;
			break;
		case 5:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x20;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xDF;
			break;
		case 6:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x40;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0xBF;
			break;
		default:
			if ( flag )
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] |= 0x80;
                        else
                                gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] &= 0x7F;
			break;
	}
}

int _Gibbsgeno_Store::get( int marker, int cycle, int phase, int individuals_ID )
{
	int pos = marker / 8;
	int off = marker % 8;
	uint8_t ret;
	switch ( off )
	{
		case 0:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x01;
			break;
		case 1:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x02;
			break;
		case 2:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x04;
			break;
		case 3:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x08;
			break;
		case 4:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x10;
			break;
		case 5:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x20;
			break;
		case 6:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x40;
			break;
		default:
			ret = gibbs_genotypes[ cycle ][ 2 * individuals_ID + phase ][ pos ] & 0x80;
	}
	return ret > 0 ? 1 : 0;
}

void _Gibbsgeno_Store::load( Init_Type & engine, int cycle )
{
	if ( snp_number != engine.snp_number || individuals != engine.individuals ) {
		cerr << "[_Gibbsgeno_Store::load] markers and samples don't match.\n";
		abort();
	}
	if ( cycle < 0 || cycle >= circles ) {
		cerr << "[_Gibbsgeno_Store::load] cycle overflow.\n";
		abort();
	}
	for ( int i = 0; i < individuals; i++ ) {
		for ( int j = 0; j < snp_number; j++ ) {
			for ( int p = 0; p < 2; p++ ) {
				int flag = engine.haplotypes[ 2 * i + p ].haplotype[ j ];
				set( j, flag, cycle, p, i );
			}
		}
	}
}

void _Gibbsgeno_Store::load( Haplotype * haplotypes, int individuals_ID, int cycle )
{
        for ( int j = 0; j < snp_number; j++ ) {
                for ( int p = 0; p < 2; p++ ) {
                        int flag = haplotypes[ 2 * individuals_ID + p ].haplotype[ j ];
                        set( j, flag, cycle, p, 0 );
                }
        }
}

void _Gibbsgeno_Store::binary2uint8_t( uint8_t *** _gibbs_genotypes )
{
	for ( int c = 0; c < circles; c++ ) {
		for ( int i = 0; i < individuals; i++ ) {
			for ( int j = 0; j < snp_number; j++ ) {
				int base1 = get( j, c, 0, i );
				int base2 = get( j, c, 1, i );
				_gibbs_genotypes[ i ][ j ][ c ] = base1 * 2 + base2;
			}
		}
	}
}

/*

#include <iostream>
using namespace std;

int main()
{
	// marker flag cycle phase individuals_ID
	_Gibbsgeno_Store GS( 100, 10000, 40 );
	for ( int c = 0; c < 40; c++ ) {
		for ( int i = 0; i < 100; i++ ) {
			for ( int j = 0; j < 10000; j++ ) {
				GS.set( j, 1, c, 1, i );
				GS.set( j, 1, c, 0, i );
			}
		}
	}
	
	for ( int c = 15; c < 16; c++ ) {
		for ( int i = 0; i < 100; i++ ) {
			cout << "Indiv " << i << " Haplo 0:";
			for ( int j = 0; j < 10000; j++ ) {
				cout << " " << GS.get( j, c, 0, i );
			}
			cout << "\n";
			cout << "Indiv " << i << " Haplo 1:";
			for ( int j = 0; j < 10000; j++ ) {
				cout << " " << GS.get( j, c, 1, i );
			}
			cout << "\n";
		}
	}
}

*/
