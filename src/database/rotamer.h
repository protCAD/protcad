#include "assert.h"
#include "typedef.h"

#ifndef ROTAMER_H
#define ROTAMER_H

class rotamer
{
public:
	friend class rotamerLib;

	rotamer();
	rotamer(const UIntVec _angles, const double _energy);
	rotamer(const rotamer& _otherRotamer); //deep copy
	~rotamer();

	DouVec getAngles() const;
	double getEnergy() const;

	void print() const;

	static double itsScaleFactor;
	static void setScaleFactor(const double _scale)
		{ itsScaleFactor = _scale; }
	static double getScaleFactor()
		{ return itsScaleFactor; }

private:
	void   setEnergy(double _energy) {itsEnergy = _energy;}
	void   setAngles(UIntVec _angles) {itsAngles = _angles;}

private:
	UIntVec itsAngles;
	double itsEnergy;
};
#endif
