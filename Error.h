#ifndef __ERROR_H__

#define __ERROR_H__

struct Error
{
	int match;
	int mismatch;
	float rate;
	Error()
	{
		match = 0;
		mismatch = 0;
		rate = 0;
	}
	Error & operator+=( Error & a )
	{
		this->match += a.match;
		this->mismatch += a.mismatch;
		return *this;
	}
	float update(){
		rate = (float)mismatch / (float)( mismatch + match );
		return rate;
	}

	void reset(){
		match = 0;
		mismatch = 0;
		rate = 0;
	}
};

struct Theta
{
	float total;
	float cross;
	Theta()
	{
		total = 0;
		cross = 0;
	}
	void reset()
	{
		total = 0;
		cross = 0;
	}

	Theta & operator +=( Theta & a )
	{
		this->total += a.total;
		this->cross += a.cross;
		return *this;
	}
	float update()
	{
		return cross / total;
	}
};

#endif
