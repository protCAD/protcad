// helixPropensity.h
//
//  Using experimentally measured delta delta G's of amino acids, the helix propensities
//  measured relative to glycine can be used to estimate the change in stability when
//  mutating from one position to another.
//
//  needs file $REPACKDIR/data/helixPropensity.frc

#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include "generalio.h"
#include "parse.h"

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef HELIX_PROPENSITY_H
#define HELIX_PROPENSITY_H

class helixPropensity
{
	public:
		
		helixPropensity();
		~helixPropensity();

		double getEnergy(const UInt _resType);
		void setHelixPropensityScaleFactor(const double _scaleFactor) {itsHelixPropensityScaleFactor = _scaleFactor; }
		double getHelixPropensityScaleFactor() {return itsHelixPropensityScaleFactor; }
		void buildDatabase();
		void printEnergyMap();

	private:

		void addToDatabase(StrVec _parsedStrings);	

	private:
		
		double itsHelixPropensityScaleFactor;
		string itsFileName;
		vector <double> itsEnergyMap;		
};

#endif		
