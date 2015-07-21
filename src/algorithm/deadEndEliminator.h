// filename: deadEndEliminator.h

#include "assert.h"
#include "typedef.h"
#include "atomIterator.h"

#ifndef DEAD_END_ELIMINATOR_H
#define DEAD_END_ELIMINATOR_H

#ifndef PROTEIN_H
#include "protein.h"
#endif

class deadEndEliminator
{
public:
	
	// Constructor and Destructor declaration
	deadEndEliminator(protein* _pProtein);
	~deadEndEliminator();

	// Accessors
	void run(UInt _level, double _cutoff);

	double min( vector< double > );
	double max( vector< double > );
	void makeRotamerDeadEnding(UInt _chain, UInt _residueIndex,
		UInt _residueType, UInt _bpt, UInt _rotamer);

private:
	//variable declarations
	UInt itsNumIterations;
	UInt itsCurrentIteration;
	protein* pItsProtein;
	UInt itsRunLevel;
	double itsEnergyCutoff;
};
#endif
