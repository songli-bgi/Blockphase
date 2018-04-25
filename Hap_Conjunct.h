#ifndef __HAP_CONJUNCT_H__

#define __HAP_CONJUNCT_H__

#include "Haplotype.h"
#include "Init_Type.h"

struct _Hap_Conjunct
{
	int overlap;
	int individuals;
	bool * orders;
	Haplotype * prev_haplotypes;
	Haplotype * next_haplotypes;

	_Hap_Conjunct( int _individuals, int _overlap )
	{
		overlap = _overlap;
		individuals = _individuals;
		prev_haplotypes = new Haplotype [ individuals * 2 ];
		next_haplotypes = new Haplotype [ individuals * 2 ];
		
		int snp_number = overlap / 2;
		Haplotype hap( snp_number );
		for ( int i = 0; i < 2 * individuals; i++ ) {
			prev_haplotypes[ i ] = hap;
			next_haplotypes[ i ] = hap;
		}

		/* true -> (0,1) | false -> (1,0) */
		orders = new bool [ individuals ];
		for ( int i = 0; i < individuals; i++ ) {
			orders[ i ] = true;
		}
	}

	void load_prev( Init_Type & engine, int begin );
	void load_next( Init_Type & engine, int begin );
	void update_orders();
	void show();

	~_Hap_Conjunct()
	{
		delete [] orders;
		delete [] prev_haplotypes;
		delete [] next_haplotypes;
	}
};

#endif
