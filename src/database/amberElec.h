//*****************    -amberElec.h-    *********************
//* AMBER partial charges for all atom and united           *
//* atom molecules found in uni_mod.in and                  *
//* all_amino94_mod.in.  In uni_mod, charges for            *
//* explicit hydrogens and lone pairs summed with           *
//* and assigned to the heavy atoms that coordinate         *
//* them. buildWithHydrogens() and buildWithoutHydrogens()  *
//* determine which charge set is used.                     *
//* getScaleFactor() and setScaleFactor() are accessors     *
//* for the scalar multipliers of the energy.               *
//*                                                         *
//* distanceDependenceOff() and distanceDependanceOn()      *
//* toggle the dielectric constant between a constant and   *
//* 40R.  The constant can be set if distance dependence    *
//* is off using setDielectricConstant(), and checked using *
//* getDielectricConstant().                                *
//***********************************************************

#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include "assert.h"
#include <string.h>

#ifndef RESIDUE_H
#include "residue.h"
#endif

#include "generalio.h"
#include "parse.h"

#ifndef AMBER_ELEC_H
#define AMBER_ELEC_H

class amberElec
{
	public:
		amberElec();
		amberElec(int _dummy);
		amberElec(const amberElec& _other);
		~amberElec();

		double getEnergy(const UInt _restype1, const UInt _atomtype1, const UInt _restype2, const UInt _atomtype2, const double _distance) const;
		double getSoluteEnergy(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distance, const double _dielectric) const;
		double getSoluteEnergySQ(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distanceSquared, const double _dielectric) const;
		double getEnergySQ(const UInt _restype1, const UInt _atomtype1, const UInt _restype2, const UInt _atomtype2, const double _distanceSquared) const;
		double getItsCharge(const UInt _restype, const UInt _atomtype) const;
		string getItsAtomName(const UInt _restype, const UInt _atomtype) const;

		void buildElectrostatics();

		void dummy();
		static void distanceDependanceOn()
			{ distanceDependance = true; cout << " amberElec distance dependance turned on" << endl;}
		static void distanceDependanceOff()
			{ distanceDependance = false; cout << " amberElec distance dependance turned off" << endl;}
		static bool isDistanceDependanceOn()
			{ return distanceDependance; }
		static void setScaleFactor (const double _scale)
			{ itsScaleFactor = _scale; }
		static double getScaleFactor()
			{ return itsScaleFactor; }
		static void setDielectricConstant(const double _dielectricConstant)
			{ itsDielectricConstant = _dielectricConstant; }
		static double getDielectricConstant()
			{ return itsDielectricConstant; }
		static void setHighEnergyCutOff(const bool _doCutOff)
			{ highEnergyCutOff = _doCutOff; return;}

		static bool distanceDependance;
		static double itsScaleFactor;
		static double itsDielectricConstant;
		static bool highEnergyCutOff;

	private:

		void buildDataBase();
		void convertToDataElements(const StrVec& _parsedStrings);
 		void orderDataElements(); 

	private:
		vector< vector< double > > charges;
		vector< vector< string > > atomNames;
		vector< string > resNames;
		string itsFileName;
};

#endif
