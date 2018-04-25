#ifndef __RNADOM_WIZZARD_H__

#define __RNADOM_WIZZARD_H__

#include <vector>
#include <boost/random.hpp>

using namespace std;
using namespace boost;

struct _Real_Uniforms
{
	int n_threads;
	unsigned long seed;
	vector<  variate_generator< kreutzer1986, uniform_real<> > > uniforms;

	_Real_Uniforms( int _n_threads, unsigned long _seed );
};

#endif
