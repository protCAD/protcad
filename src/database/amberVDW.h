// filename: amberVDW.h

#include "assert.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "typedef.h"
#include <stdio.h>
#include <math.h>
#include "generalio.h"
#include "parse.h"

#ifndef _AMBER_VDW_H
#define _AMBER_VDW_H

class amberVDW
{
public:
	amberVDW();
	amberVDW(int _dumy);
	amberVDW(const amberVDW& _otherAmberVDW);
	~amberVDW();
	
	double getEnergy(const UInt _type1,const UInt _type2,const double _distance) const;
    double getWaterEnergy(const UInt _type1) const;
	double getEnergySQ(const UInt _type1, const UInt _type2, const double _distanceSquared) const;
	int getIndexFromNameString(string _name);
	bool isClash(const UInt _type1, const UInt _type2, const double _distance);
    double getRadius(const UInt _type1);
    double getPolarizability(const UInt _type1);
    double getVolume(const UInt _type1);

	static double itsScaleFactor;
	static void setScaleFactor(const double _scale)
		{ itsScaleFactor = _scale;
                  // cout << " vdw scale factor set to: " << _scale << endl;
                 }
	static double getScaleFactor()
		{ return itsScaleFactor; }

	static double itsRadiusScaleFactor;
	static void setRadiusScaleFactor(const double _scale)
		{ itsRadiusScaleFactor = _scale;
                //  cout << " vdW radius scale factor set to: " << _scale << endl;

                }
	static double getRadiusScaleFactor()
		{ return itsRadiusScaleFactor; }

	static bool linearRepulsionDampening;	// Kuhlman & Baker PNAS v97 p10383
	static void setLinearRepulsionDampeningOn()
		{ linearRepulsionDampening = true; 
                //  cout << " vdW linear repulsion dampening activated." << endl; 

		}
	static void setLinearRepulsionDampeningOff()
		{ linearRepulsionDampening = false; 
		//cout << " vdw linear repulsion dampening deactivated." << endl; 
		}

	static double itsAttractionScaleFactor;
	static void setAttractionScaleFactor( const double _scale)
		{ itsAttractionScaleFactor = _scale; 
	//		cout << " vdW attraction scale factor set to: " << itsAttractionScaleFactor << endl; 
		}
	static double getAttractionScaleFactor()
		{ return itsAttractionScaleFactor; }
	
	static double itsRepulsionScaleFactor;
	static void setRepulsionScaleFactor( const double _scale)
		{ itsRepulsionScaleFactor = _scale; 
		//cout << " vdW repulsion scale factor set to: " << itsRepulsionScaleFactor << endl; 
		}
	static double getRepulsionScaleFactor()
		{ return itsRepulsionScaleFactor; }
	
private:
	void buildDataBase();
	// StrVec is intended non const reference
	void convertToDataElements(const StrVec& _parsedStrings);

private:
	vector< double > R_ref;
	vector< double > EPS;
    vector< double > Pol_ref;
    vector< double > Vol_ref;
	vector< string > amberAtomTypeNames;
	string itsFileName;
};

#endif
