#include "Sample_Orders.h"

_Sample_Orders::_Sample_Orders( int individuals ):uniform( kreutzer1986(2000), uniform_real<>(0,1) )
{
	_individuals = individuals;
	for ( int i = 0; i < _individuals; i++ ) {
		orders.push_back( i );
	}
	_circle = 100000;
}

void _Sample_Orders::randOrder()
{
	int a, b, t;
	for ( int i = 0; i < _circle; i++ ) {
		a = static_cast< int >( uniform() * _individuals );
		b = static_cast< int >( uniform() * _individuals );
		t = orders[ a ];
		orders[ a ] = orders[ b ];
		orders[ b ] = t;
	}
}

void _Sample_Orders::showOrder()
{
	for ( int i = 0; i < _individuals; i++ ) {
		cerr << "\t" << orders[ i ];
	}
	cerr << "\n";
}

/*

int main()
{
	_Sample_Orders SO( 100 );
	SO.showOrder();
	SO.randOrder();
	SO.showOrder();
	SO.randOrder();
	SO.showOrder();
	return 0;
}

*/
