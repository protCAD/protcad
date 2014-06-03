// filename: annealer.h

#include "assert.h"
#include "typedef.h"

#ifndef ANNEALER_H
#define ANNEALER_H

#include "ran.h"

#ifndef ENSEMBLE_H
#include "ensemble.h"
#endif

class annealer
{
public:

	// Constructor and Destructor declaration
	annealer(ensemble* _pEnsemble);
	~annealer();

	// Accessors
	void run(double _tempH, double _tempL, UInt _iter, UInt _seed);
	void run(double _tempH, double _tempL, UInt _iter);
	static UInt getHowMany() {return howMany; }
	void writeEveryState();
	void dontWriteEveryState();

private:
	int   modifyEnsemble(ran& _ran);
	float firstrun;
	double energy();
	void   decrementTemp();
	void   acceptModification();
	void   rejectModification();
	void   setupSystem(ran& ran);
	void   saveState(string& _filename);

private:
	//variable declarations
	static UInt howMany;
	double itsTemp;
	double itsTempDecrement;
	UInt itsNumIterations;
	UInt itsCurrentIteration;
	UInt itsSeed;
	UInt itsNumNormalKeep;
	UInt itsNumMetropolisKeep;
	double itsNewEnergy;
	double itsOldEnergy;
	double itsProbAccept;
	double e;
	double itsLowestEnergy;
	ensemble* pItsEnsemble;
	ran itsRan;
	bool writeEveryStateFlag;
};
#endif
