// filename: microEnvDB.h

#include "assert.h"
#include "enums.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "typedef.h"
#include "generalio.h"

#include <stdlib.h>
#include <stdio.h>

#ifndef MICRO_ENV_DB_H
#define MICRO_ENV_DB_H

class microEnvDB
{
public:
	microEnvDB();
	microEnvDB(double _rad, UInt const _skip);
	~microEnvDB();
	
	double getEnergy(const UInt _type1, const UIntVec& cluster) const;

	static vector< vector< vector< UInt > > > itsMicroEnvMap;
	static	void printMicroEnvMap();
	void printEnergyData();

	static void setUseSingleBodyEnergyOn()
		{useSingleBodyEnergy = true; 
		//cout << "Using singleBodyEnergy for unspecified microenvironments" << endl;
		}
	static void setUseSingleBodyEnergyOff()
		{useSingleBodyEnergy = false;
		//cout << "NOT using singleBodyEnergy for unspecified microenvironments" << endl;
		}
	
	static bool useSingleBodyEnergy;
private:
	void readClusterData();
	void readEnergyData();
	void readSingleBodyEnergyData();

	double findEnergyInDB(const UInt _type, const UIntVec& _cluster) const;

private:
	vector< vector< vector< double> > > energyData;
	vector< vector< double > > itsSingleBodyEnergy;
	double itsCriticalRadius;
	UInt itsResidueSkippingNumber; 
};
#endif
