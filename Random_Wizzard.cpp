#include "Random_Wizzard.h"

_Real_Uniforms::_Real_Uniforms( int _n_threads, unsigned long _seed )
{
	n_threads = _n_threads;
	seed = _seed;
	for ( int i = 0; i < n_threads; i++ ) {
		kreutzer1986 p( seed * (i + 1) );
		uniform_real<> r(0, 1);
		variate_generator< kreutzer1986, uniform_real<> > u( p, r );
		uniforms.push_back( u );
	}
}
