#include "ramachandranMap.h"

ramachandranMap::ramachandranMap()
{	initialize();
}

ramachandranMap::~ramachandranMap()
{
}

void ramachandranMap::initialize()
{
	// 0 = alpha right	-60,-40 +/- 20
	// 1 = gamma		-80,0	+/- 20
	// 2 = beta			-90,120	+/- 40
	// 3 = alpha left	90,0	+/- 30
	// 4 = other

	itsAllowedRegions.resize(0);
	itsRegionsUpperRight.resize(0);
	itsRegionsLowerLeft.resize(0);

	for (UInt i=0; i<5; i++)
	{	itsAllowedRegions.push_back(i);
	}

	// alpha right
	vector<double> tempvec1;
	tempvec1.resize(0);
	tempvec1.push_back(-40.0);
	tempvec1.push_back(-20.0);
	itsRegionsUpperRight.push_back(tempvec1);
	vector<double> tempvec2;
	tempvec2.resize(0);
	tempvec2.push_back(-80.0);
	tempvec2.push_back(-60.0);
	itsRegionsLowerLeft.push_back(tempvec2);
	// gamma
	vector<double> tempvec3;
	tempvec3.resize(0);
	tempvec3.push_back(-60.0);
	tempvec3.push_back(20.0);
	itsRegionsUpperRight.push_back(tempvec3);
	vector<double> tempvec4;
	tempvec4.resize(0);
	tempvec4.push_back(-100.0);
	tempvec4.push_back(-20.0);
	itsRegionsLowerLeft.push_back(tempvec4);
	// beta
	vector<double> tempvec5;
	tempvec5.resize(0);
	tempvec5.push_back(-50.0);
	tempvec5.push_back(160.0);
	itsRegionsUpperRight.push_back(tempvec5);
	vector<double> tempvec6;
	tempvec6.resize(0);
	tempvec6.push_back(-130.0);
	tempvec6.push_back(80.0);
	itsRegionsLowerLeft.push_back(tempvec6);
	// alpha left 
	vector<double> tempvec7;
	tempvec7.resize(0);
	tempvec7.push_back(120.0);
	tempvec7.push_back(30.0);
	itsRegionsUpperRight.push_back(tempvec7);
	vector<double> tempvec8;
	tempvec8.resize(0);
	tempvec8.push_back(60.0);
	tempvec8.push_back(-30.0);
	itsRegionsLowerLeft.push_back(tempvec8);
}

int ramachandranMap::getSecondaryStructureIndex(double _phi, double _psi)
{
	int index = 4;
	for (UInt i=0;i<4;i++)
	{
		if (_phi <= itsRegionsUpperRight[i][0] &&
		    _phi >= itsRegionsLowerLeft[i][0])
		{
			if (_psi <= itsRegionsUpperRight[i][1] &&
			    _psi >= itsRegionsLowerLeft[i][1])
			{
				index = i;
				return index;
			}
		}
	}				
	return index;
}
