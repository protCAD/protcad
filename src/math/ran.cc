#include "ran.h"

ran::ran()
{	itsRan1 = new ran1(0);
}

ran::ran(UInt _seed)
{//	cout << "ran constructor called:"
//	     << "ran::ran()" << endl;
	itsRan1 = new ran1(_seed);
}

ran::~ran()
{//	cout<< "ran destructor called " << endl;
	delete itsRan1;
}

void ran::setSeed(unsigned int _seed)
{	itsRan1->setSeed(_seed);
}

double ran::getNext()
{	return itsRan1->getNext();
}

int ran::getNext(const int _lowerBound, const int _upperBound)
{	return itsRan1->getNext(_lowerBound, _upperBound);
}

double ran::getNext(const double _lowerBound, const double _upperBound)
{	return itsRan1->getNext(_lowerBound, _upperBound);
}
