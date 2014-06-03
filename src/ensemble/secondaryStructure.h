// filename secondaryStructure.h

#include "assert.h"
#include "typedef.h"
#include "ramachandranMap.h"

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef SECONDARYSTRUCTURE_H
#define SECONDARYSTRUCTURE_H

class secondaryStructure 
{
public:
	secondaryStructure();
	secondaryStructure(int _indexInChain);
	secondaryStructure(const secondaryStructure& _rhs);
	~secondaryStructure();
	
	int getSecondaryStructure() {return itsSecondaryStructureIndex;}
	UInt getStartingResidue() {return itsStartingResidueIndex;}
	UInt getEndingResidue() {return itsEndingResidueIndex;}
	void setStartingResidue(UInt _index);
	void setEndingResidue(UInt _index);

	void assign(chain* _pTheChain);
	
	static ramachandranMap theMap;

private:
	void assignViaPhiPsi(chain* _pTheChain);
	UInt itsStartingResidueIndex;
	UInt itsEndingResidueIndex;
	int itsSecondaryStructureIndex;
	int itsIndexInChain;
};
#endif
