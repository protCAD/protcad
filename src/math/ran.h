// filename: ran.h

#include "typedef.h"
#include <iostream>
#include "ran1.h"

#ifndef RAN_H
#define RAN_H
// Important: when using ran1, seed may not exceed the value 54773 !!!!
// Generates uniform deviate between 0.0 and 1.0
class ran
{
public:
	// Constructor and Destructor declaration
	ran();
	ran(UInt _seed);
	~ran();

	// Accessors
	void setSeed(unsigned int _seed);
    //static UInt getHowMany() {return ran::getHowMany();}
	double getNext();
	double getNext(double _lowerBound, double _upperBound);
	int getNext(int _lowerBound, int _upperBound);

private:
	ran1* itsRan1;
};
#endif
