#ifndef __HAPLOTYPE_H__

#define __HAPLOTYPE_H__

#include <stdint.h>
#include <iostream>
#include <string.h>

using namespace std;

struct Haplotype
{
        uint8_t * haplotype;
        int length;
        Haplotype()
        {
                haplotype = NULL;
                length = 0;
        }

        Haplotype( int snp_number )
        {
                length = snp_number;
                haplotype = new uint8_t[ length ];
                memset( haplotype , 0 , length );
        }

        Haplotype( const Haplotype & h )
        {
                length = h.length;
                haplotype = new uint8_t[ length ];
                memcpy( haplotype , h.haplotype , length );
        }

        inline bool operator==( const Haplotype & hap )
        {
                for ( int i = 0; i < this->length; i++ ){
                        if ( this->haplotype[i] ^ hap.haplotype[i] != 0 ){
                                return false;
                        }
                }
                return true;
        }

        inline void operator=( const Haplotype & hap )
        {
                this->length = hap.length;
                this->haplotype = new uint8_t[ hap.length ];
                memcpy( this->haplotype , hap.haplotype , hap.length );
        }

	friend ostream & operator << ( ostream & os, Haplotype & h )
        {
                os << "Haplotype show:\t|";
                for ( int i = 0; i < h.length; i++ ) {
                        os << int(h.haplotype[ i ]);
                }
                os << "\n";
                return os;
        }

        ~Haplotype()
        {
                delete [] haplotype;
        }
};

#endif
