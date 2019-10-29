// filename: ran.h

#include "typedef.h"
#include <iostream>

#ifndef RAN_H
#define RAN_H
// Important: when using ran1, seed may not exceed the value 54773 !!!!
// Generates uniform deviate between 0.0 and 1.0
class ran1
{
public:
	
	// Constructor and Destructor declaration
	ran1(unsigned int _seed);
	~ran1();

	// Accessors
	void setSeed(unsigned int _seed);
	void initialize(unsigned int _seed);
	static unsigned int getHowMany() {return howMany; }
	double getNext();
	double getNext(double _lowerBound, double _upperBound);
	int getNext(int _lowerBound, int _upperBound);

private:
	//variable declarations
	static unsigned int howMany;
	unsigned int idum;
	unsigned int M1;
	unsigned int IA1;
	unsigned int IC1;
	double RM1;
	unsigned int M2;
	unsigned int IA2;
	unsigned int IC2;
	double RM2;
	unsigned int M3;
	unsigned int IA3;
	unsigned int IC3;
	double r[98];
	long ix1;
	long ix2;
	long ix3;
};
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
