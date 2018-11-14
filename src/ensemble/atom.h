// ******************************************************************
// ******************************************************************
// 	filename: atom.h
// 	contents: class atom is defined
// 	Last Modified: 01/05/2000
// ******************************************************************
// ******************************************************************

#include "assert.h"
#include <vector>
#include <string>
#include <iostream>
#include "typedef.h"
#include "pdbData.h"
#include "PDBAtomRecord.h"
#include "treeNode.h"
#include "generalio.h"
#include "point.h"
#include "unitSphere.h"

#ifndef RESIDUE_H
class residue;
#endif

#ifndef ATOMITERATOR_H
class atomIterator;
#endif

#ifndef ATOM_H
#define ATOM_H

class atom : public treeNode, public point
{
// ***********************************************************************
// ***********************************************************************
// 	Nested Types Under Atom ---- Declarations
// ***********************************************************************
// ***********************************************************************

protected:
	class typeInfo; // defined at the end

// ***********************************************************************
// ***********************************************************************
// 	Making friends with these fellows
// ***********************************************************************
// ***********************************************************************

	friend class residue;
	friend class atomIterator;
	friend class unitSphere;

// ***********************************************************************
// ***********************************************************************
//	Contructors and Destructors
// ***********************************************************************
// ***********************************************************************

public:
	// Constructor and Destructor declaration
	atom();
	atom(const string& _atomType);
	atom(const UInt _atomType);
	atom(const pdbAtom& _pdbAtomData);
	atom(const PDBAtomRecord& _theRecord, bool _hetflag);
	atom(const PDBAtomRecord& _theRecord);
	atom(const atom& _otherAtom);
	~atom();

protected:
	void atomDefaultValues();
	// auxilary function for the contructors, so set as private


// ***********************************************************************
// ***********************************************************************
//	Geometry Related Operations
// ***********************************************************************
// ***********************************************************************



public:
	// All Coordinate Related Operations are same as in point.h
	double distance(const atom& otherAtom) const;
	double distance(atom* pOtherAtom) const;
	double distanceSquared(const atom& otherAtom) const;
	double distanceSquared(atom* pOtherAtom) const;
	bool inCutoff (const atom* pOtherAtom, double _cutoffDistance);
	double inCubeWithDist (const atom* pOtherAtom, double _cutoff);
    double inCubeWithDistSQ (const atom* pOtherAtom, double _cutoff);
	bool inCube (const atom* pOtherAtom, double _cutoffDistance);
	bool inCutoffSQ (const atom* pOtherAtom, double _cutoff, double _cutoffSquared);
	bool inCutoff (const atom& OtherAtom, double _cutoffDistance);
	bool inCutoffSQ (const atom& OtherAtom, double _cutoff, double _cutoffSquared);
	// Other
 	bool isFullySpecified() const {return fullySpecified;}
 	void setIsFullySpecified() {fullySpecified = true;}
 	void setIsNotFullySpecified() {fullySpecified = false;}
	void sethetatmFlag(bool _flag) {hetatmFlag=_flag;}
// ***********************************************************************
// ***********************************************************************
//	Atom Properties`Accessors And Modifiers
// ***********************************************************************
// ***********************************************************************

public:
	// Radius Related Operations
	double getRadius() const {return itsRadius; }
	double getEpsilon() const {return itsEpsilon; }
    double getPolarizability() const {return itsPolarizability; }
    double getVolume() const {return itsVolume; }
	void setRadius(const double _radius);
	void setAtomicRadius(double _radius);
	double getSolvationEnergy() {return itsSolvationEnergy;}
	double getDielectric() {return itsDielectric;}
    double getNumberofWaters() {return itsWaters;}
	double getMaxDielectric() {return itsMaxDielectric;}
	double getMinDielectric() {return itsMinDielectric;}
	void setSolvationEnergy(double _solvationEnergy);
	void setDielectric(double _dielectric);
    void setNumberofWaters(double _waters);
	void setMaxDielectric(double _maxDielectric);
	void setMinDielectric(double _minDielectric);

	// Charge
	int   getSerialNumber() const {return itsSerialNumber;}
	void  setSerialNumber(const UInt _s);
	double getOccupancy() const {return itsOccupancy;}
	double getTempFactor() const  {return itsTempFactor;}
	void  setCharge(const double _c);
	double getCharge() const {return itsCharge; }
	char  getChainID() const {return itsChainID;}
	string getLigChainID() const {return itsLigChainID;}
	void  setChainID(char _ID) {itsChainID = _ID; }
        void  setLigChainID(string _ID){itsLigChainID= _ID;}
	string getResType() const;
	bool getSilentStatus() { return isSilent; }

protected:
	void  setOccupancy(const double _o);
	void  setTempFactor(const double _tf);

// ***********************************************************************
// ***********************************************************************
//	Atom Type Operations
// ***********************************************************************
// ***********************************************************************

protected:
	void defineType(const UInt _atomType );
	void defineType(const string& _atomType );
	void setType(const string& _atomType );
	void setType(const UInt _atomType );
	void makeAtomSilent();
public:
	// Note: this function returns a chemical "type",
	// i.e., C,O,N,S, etc.
	string getType() const;

protected:
	// Key module that builds the atom type Base
	// Not meant to be used by users and one copy suffices
	static void buildDataBase();


// ***********************************************************************
// ***********************************************************************
//	Atom Name Operations
// ***********************************************************************
// ***********************************************************************

public:
	//Note: this function returns an atom "name" which depends
	// on its position, etc. in the molecule
	// i.e. CA,CB,ND1, etc.
	string getName() const;
protected:
	void setName(const int _atomName) {itsName = _atomName;}
	// Name should only be set from residue which has access to this
	// function


// ***********************************************************************
// ***********************************************************************
// 	Static Variable Accessors
// ***********************************************************************
// ***********************************************************************

public:
	// It can be called even when no object is instantiated
	static UInt getHowMany() {return howMany; }


// ***********************************************************************
// ***********************************************************************
// 	Variables Declartions
// ***********************************************************************
// ***********************************************************************

protected:
	// non-static variable declaration
	double itsRadius;
	double itsSolvationEnergy;
	double itsDielectric;
    double itsWaters;
	double itsMaxDielectric;
	double itsMinDielectric;
	double itsEpsilon;
    double itsPolarizability;
    double itsVolume;
	bool isSilent;
	UInt itsType; // e.q. N, C, P ..
	int itsName; // e. q. NH1, CA, CB ...
	int itsAtomEnergyType;
	bool fullySpecified;
	
	//hetatm variables, not used for proteins
	bool hetatmFlag;
	string itsLigName;  //used in place of itsName
	string itsResName;   //used in place of itsResType
	string itsLigChainID; //used in place of itsChainID
        int itsAmberAllType;
        int itsAmberUnitedType;
        double itsAmberAllCharge;
        double itsAmberUnitedCharge;
	

	UInt itsSerialNumber;

	// double precision is good enough and the size of the atom
	// construct can be reduced as well
	double itsOccupancy;
	double itsTempFactor;
	double itsCharge;
	int itsResType;
	char  itsChainID;

	//static variable declaration
	static UInt howMany;
	static vector<typeInfo> dataBase;
	static bool dataBaseBuilt;

// **********************************************************************
// **********************************************************************
// spherePoint and surface area operations
// **********************************************************************
// **********************************************************************

public:
	void initializeSpherePoints();
	UInt getHowManySpherePoints() { return itsSpherePointFlags.size(); }
	void removeSpherePoints(atom* _pOtherAtom);
	double getItsProbeRadius() { return itsProbeRadius; }
	void setItsProbeRadius(double _probeRadius) { itsProbeRadius = _probeRadius; }
	double calculateTotalSASA();
	double calculateExposedSASA();	
	double calculateBuriedSASA();
private:
	vector <bool> itsSpherePointFlags;
	static double itsProbeRadius;

	void removeSpherePoint(UInt _index);
	UInt countTrueSpherePoints(); 
	void generateNewSpherePoints(dblVec _center, double _radius);
	void generateNewSpherePoints(); 
};

// ***********************************************************************
// ***********************************************************************
// 	
//  Nested Types Under Atom ------ Definitions
// ***********************************************************************
// ***********************************************************************

class atom::typeInfo
{
public:
	typeInfo()
    {	vdwRadius.resize(2);
	}

	string typeName;
	vector<double> vdwRadius;
	// 0: radius, 1:lower bound, 2:upper bound
};

#endif

// ***********************************************************************
// ***********************************************************************
// END OF FILE
// ***********************************************************************
// ***********************************************************************
