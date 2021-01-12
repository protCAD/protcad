#include "ran.h"
#include "math.h"

unsigned int ran1::howMany=0;

ran1::ran1(unsigned int _seed)
{	
#ifdef RAN1_DEBUG
	cout << "ran1 constructor called:"
	     << "ran1::ran1()" << endl;
#endif
	if (_seed > 54773 )
	{
		cout << "Error from ran1: random seed may not exceed 54773" << endl;
		exit(1);
	}
	initialize(_seed);
	howMany++;
}

ran1::~ran1()
{
#ifdef RAN1_DEBUG
	cout<< "ran1 destructor called " << endl;
#endif
        howMany--;
}

void ran1::setSeed(unsigned int _seed)
{
	if (_seed > 54773 )
	{
		cout << "Error from ran1: random seed may not exceed 54773" << endl;
		exit(1);
	}
    initialize(_seed);
}

void ran1::initialize(unsigned int _seed)
{
	idum = _seed;
	M1 = 259200;
	IA1 = 7141;
	IC1 = 54773;
	RM1 = 1.0 / M1;
	M2 = 134456;
	IA2 = 8121;
	IC2 = 28411;
	RM2 = 1.0 / M2;
	M3 = 243000;
	IA3 = 4561;
	IC3 = 52349;

	ix1 = (IC1 - idum) % M1;
	ix1 = (IA1 * ix1 + IC1) % M1;
	ix2 = ix1 % M2;
	ix1 = (IA1 * ix1 + IC1) % M1;
	ix3 = ix1 % M3;
	for (int j=1; j<=97; j++)
	{
		ix1 = (IA1 * ix1 + IC1) % M1;
		ix2 = (IA2 * ix2 + IC2) % M2;
		r[j] = (ix1+ix2*RM2) * RM1;
	}
	idum = 1;
}

double ran1::getNext()
{
	ix1 = (IA1 * ix1 + IC1) % M1;
	ix2 = (IA2 * ix2 + IC2) % M2;
	ix3 = (IA3 * ix3 + IC3) % M3;
	int j;
	j = 1 + ((97*ix3)/M3);
	if (j>97 || j<1)
	{
		cout  << "ran1: fatal error!" << endl;
	}
	double temp = r[j];
	r[j] = (ix1+ix2*RM2) * RM1;
	return temp;
}

int ran1::getNext(int _lowerBound, int _upperBound)
{
	double tempreal = getNext();
	int tempupper = _upperBound;
	int templower = _lowerBound;
	int temprange = 0;
	int tempint = 0;
	if (tempupper == templower)
	{
	//	cout << "Error returned by ran1::getNext(int,int): " << endl;
	//	cout << "Upper and lower bounds are equal!" << endl;
		return tempupper;
	}

	if (tempupper < templower)
	{
		// Swap 'em!
		int temptemp = templower;
		templower = tempupper;
		tempupper = temptemp;
	}

	temprange = tempupper - templower;
	tempreal = (tempreal * (temprange + 1.0)) + templower;
	tempint = (int)tempreal;
	return tempint;
}

double ran1::getNext(double _lowerBound, double _upperBound)
{
	double tempreal = getNext();
	double tempupper = _upperBound;
	double templower = _lowerBound;
	double temprange = 0;
	double tempdouble = 0;
	if (tempupper == templower)
	{
	//	cout << "Error returned by ran1::getNext(double,double): " << endl;
	//	cout << "Upper and lower bounds are equal!" << endl;
		return tempupper;
	}

	if (tempupper < templower)
	{
		// Swap 'em!
		double temptemp = templower;
		templower = tempupper;
		tempupper = temptemp;
	}

	temprange = tempupper - templower;
	tempdouble = (tempreal * temprange) + templower;
	return tempdouble;
}

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

