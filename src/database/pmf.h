// filename: pmf.h

#include "assert.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include "generalio.h"
#include "parse.h"

#ifndef PMF_H
#define PMF_H

class pmf
{
public:
	pmf(string _fname);
	pmf(const pmf& _other);
	~pmf();
	
	double getEnergy(const UInt _type1,const UInt _type2,const double _distance) const;

	static double itsScaleFactor;
	static void setScaleFactor(const double _scale)
		{ itsScaleFactor = _scale;}
	static double getScaleFactor()
		{ return itsScaleFactor; }
	int getIndexFromNameString(const string _name, const UInt _searchField) const;
	static vector< vector< string > > atomTypeNameStrings;
	static bool atomTypeNameStringsBuilt;
	
private:
	int findDistanceBin(const double _distance) const;
	void buildDataBase();
	// StrVec is intended non const reference
	void convertToDistanceBinUpperLimits(const StrVec& _parsedStrings);
	void convertToEnergyDataBase(const StrVec& _parsedStrings);
	UInt binarySearch(const UInt low, const UInt high, const double _distance) const;
	UInt linearSearch(const double _distance) const;
	void readAtomTypeData();

private:
	vector< vector<DouVec> > pairwiseEnergyData;
	DouVec distanceBinUpperLimits;
	double stepSize;
	double lowerDistanceLimit;
	double upperDistanceLimit;
	string itsFileName;
};

#endif
