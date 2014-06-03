// filename: microEnvironment.h

#include "assert.h"
#include "pmf.h"
#include "enums.h"
#include <string.h>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include "molecule.h"
#include "protein.h"
#include "microEnvDB.h"
#include "atomIterator.h"

#ifndef MICRO_ENVIRONMENT_H
#define MICRO_ENVIRONMENT_H

class microEnvironment
{
public:
	microEnvironment();
	microEnvironment(const double _rad, const UInt _skip);
	~microEnvironment();

	void initialize();
	
	double calculateEnergy();
	double calculateEnergy(molecule* _pMol);

	void setDistanceDependenceOn() {itsDistDepFlag = true;}
	void setDistanceDependenceOff() {itsDistDepFlag = false;}

	void printFlagStatus() const;

	vector< double >  getEnergyPerAtom() const;
	vector< vector <UInt> > getEnvPerAtom() const { return environments;}
	vector< UInt > getEnvAssignedAtomTypes() const;

	microEnvDB* pItsDB;
	pmf* pItsLoResPMF;

	static double itsScaleFactor;
	static void setScaleFactor(const double _scale)
		{ itsScaleFactor = _scale; }
	static double getScaleFactor()
		{ return itsScaleFactor; }

private:
	molecule* pItsMol;
	vector<double> itsEnergyPerAtom;
	vector< vector <UInt> > environments;
	vector< vector <int> > theAtomTypes;
	double itsEnergySum;
	double itsCriticalRadius;
	UInt itsResidueSkippingNumber; 
	bool itsDistDepFlag;

};
#endif
