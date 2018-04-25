#ifndef __THREAD_SYNCHRONIZATION_H__

#define __THREAD_SYNCHRONIZATION_H__

struct _Thread_Syn
{
	bool * flags;
	int n_threads;
	
	_Thread_Syn( int _n_threads )
	{
		n_threads = _n_threads;
		flags = new bool [ n_threads ];
		for ( int i = 0; i < n_threads; i++ ) {
			flags[ i ] = false;
		}
	}

	void reset()
	{
		for ( int i = 0; i < n_threads; i++ ) {
			flags[ i ] = false;
		}
	}

	void set( int i )
	{
		flags[ i ] = true;
	}

	bool is_Can_go()
	{
		for ( int i = 0; i < n_threads; i++ ) {
			if ( ! flags[ i ] ) {
				return false;
			}
		}
		return true;
	}

	~_Thread_Syn()
	{
		delete [] flags;
	}
};

#endif
