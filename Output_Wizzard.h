#ifndef __OUTPUT_WIZZARD_H__

#define __OUTPUT_WIZZARD_H__

#include <string>
#include <zlib.h>
#include "Init_Type.h"
#include "Run_Wizzard.h"
#include "sq_Pearson_cc.h"
#include "Hap_Conjunct.h"
#include "Gibbsgeno_Store.h"

struct _Output_Wizzard
{
	string outfile;
	int individuals;
	_Hap_Conjunct * _HC;
	
	_Output_Wizzard( string & _outfile, int _individuals, _Hap_Conjunct * _hc )
	{
		outfile = _outfile;
		individuals = _individuals;
		_HC = _hc;
	}

	_Output_Wizzard()
	{
		_HC = NULL;
		individuals = 0;
	}

	void phased_out( _Gibbsgeno_Store & GS, Init_Type & engine, int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start, int cycle );
	void gprobs_out( _Gibbsgeno_Store & GS, int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start );
	void dosage_out( _Gibbsgeno_Store & GS, \
			 int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start );
	void r2_out( _Gibbsgeno_Store & GS, Init_Type & engine, \
		     int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start );

	void process( _Gibbsgeno_Store & GS, Init_Type & engine, \
		      int out_begin, int out_end, bool add, _Run_Wizzard & RW, int start, int cycle );
	
	void decide_all(_Gibbsgeno_Store & GS, Init_Type & engine, int cycle);
	void decide_homo(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx);
	void decide_hete(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx);

	void decide_one( _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, int individuals_ID, \
			 variate_generator< kreutzer1986, uniform_real<> > & uniform );
	void decide_homo(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx, int cur_id, \
			 variate_generator< kreutzer1986, uniform_real<> > & uniform);
	void decide_hete(int i, _Gibbsgeno_Store & GS, Init_Type & engine, int cycle, vector<int> & hete_idx, int cur_id);

	int max_id( uint8_t * g );
};

#endif
