// filename: aaBaseline.h

#include "assert.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include "generalio.h"
#include "parse.h"
//#include "CMath.h"

#ifndef AABASELINE_H
#define AABASELINE_H

class aaBaseline
{
public:
	aaBaseline();
	aaBaseline(int _dummy);
	aaBaseline(const aaBaseline& _other);
	~aaBaseline();
	
	double getEnergy(const string& _name) const;
    vector <string> list() const;
	static double itsScaleFactor;
	static void setScaleFactor(const double _scale)
		{ itsScaleFactor = _scale;}
	static double getScaleFactor()
		{ return itsScaleFactor; }
	
private:
	int getIndexFromNameString(string _name);
	void buildDataBase();
	void convertToDataElements(const StrVec& _parsedStrings);

private:
	vector< double > baselineData;
	vector< string > residueNameStrings;
	string itsFileName;
};

#endif
