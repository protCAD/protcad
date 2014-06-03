//filename:  referenceSolvation.h

#include "assert.h"
#include <string.h>
#include <fstream.h>
#include <vector.h>
#include "typedef.h"
#include <stdio.h>
#include <math.h>
#include "generalio.h"
#include "parse.h"

#ifndef REFERENCESOLVATION_H
#define REFERENCESOLVATION_H
class referenceSolvation
{
public:
	referenceSolvation();
	~referenceSolvation();

	void buildReferenceSolvationTable();
	void printTable();
	double getReferenceSolvationEnergy(UInt _resType, UInt _param);
	
private: 
	void convertToDataElements(StrVec& _strings);

private:
	vector < vector <double> > its Params;
	StrVec resTypeList;
	string dataFileName;

};

#endif
