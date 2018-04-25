#ifndef __SAMPLE_ORDERS_H__

#define __SAMPLE_ORDERS_H__

#include <vector>
#include <boost/random.hpp>

using namespace std;
using namespace boost;

struct _Sample_Orders
{
private:
	int _individuals;
	variate_generator< kreutzer1986, uniform_real<> > uniform;
	int _circle;
public:
	vector< int > orders;
	_Sample_Orders( int individuals );
	
	void randOrder();
	void showOrder();
};


#endif
