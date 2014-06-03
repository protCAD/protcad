// filename:  solvation.h

#include "assert.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include <math.h>
#include "generalio.h"
#include "parse.h"

#ifndef SOLVATION_H
#define SOLVATION_H
class solvation 
{
	
public:
	solvation();
	~solvation();

	void buildSolvationDataBase();	
	int getIndexFromNameString(string _name);
	void printParamTable();
	
	double getSolvationEnergy(double _surfaceArea, UInt _atomType, UInt _paramSet);

	
	static double getItsScaleFactor() { return itsScaleFactor; }
	static void setItsScaleFactor(double _scaleFactor) { itsScaleFactor = _scaleFactor; }


private:  //operations
	
	void convertToDataElements(StrVec& _strings);
	

private:  //variables
	static double itsScaleFactor;
	vector < vector <double> > itsParams;
	StrVec atomTypeList;
	string dataFileName;
};

#endif
