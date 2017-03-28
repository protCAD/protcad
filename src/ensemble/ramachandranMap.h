// filename ramachandranMap.h

#include "assert.h"
#include <vector>
#include "typedef.h"

#ifndef RAMACHANDRANMAP_H
#define RAMACHANDRANMAP_H

class ramachandranMap 
{
public:
	ramachandranMap();
	~ramachandranMap();
	
	void initialize();
	int getSecondaryStructureIndex(double _phi, double _psi);

private:
	UIntVec itsAllowedRegions;
	vector< vector<double> > itsRegionsUpperRight;
	vector< vector<double> > itsRegionsLowerLeft;
};
#endif
