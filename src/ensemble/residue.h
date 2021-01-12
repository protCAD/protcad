// ***********************************************************************
// ***********************************************************************
// 	filename: residue.h
// 	contents: class residue is defined
// 	Last Modified: 1/5/2000
// ***********************************************************************
// ***********************************************************************

#include "assert.h"
#include <string.h>
#include <vector>
#include "generalio.h"
#include "atom.h"
#include "typedef.h"
#include "PDBAtomRecord.h"

//#ifndef RESIDUETEMPLATE_H
//#include "residueTemplate.h"
//#endif

#ifndef MOLECULE_H
#include "molecule.h"
#endif

class residueTemplate;

#ifndef RESIDUE_H
#define RESIDUE_H
#ifndef CHAINPOSITION_H
#include "chainPosition.h"
#endif

#ifndef ROTAMERLIB_H
#include "rotamerLib.h"
#endif

#ifndef CHAIN_H
class chain;
#endif

#ifndef ATOMITERATOR_H
class atomIterator;
#endif

#ifndef PROTEIN_H
class protein;
#endif

class residue
{
// ***********************************************************************
// ***********************************************************************
// 	Making friends with these folks
// ***********************************************************************
// ***********************************************************************

	friend class atom;
	friend class chain;
	friend class atomIterator;
	friend class chainPosition;
	friend class allowedResidue;
	friend class rotamerLib;
	friend class protein;
	friend class ruler;
	//friend molecule* pdbReader(const string& _filename);

// ***********************************************************************
// ***********************************************************************
// 	Constructors and Destructors
// ***********************************************************************
// ***********************************************************************

public:
	residue();
	residue(const string& _aaType);
	residue(const string& _aaType, const bool _hFlag);
	residue(const UInt _itsType);
	residue(const UInt _itsType, const bool _hFlag);
	residue(const residue& _rhs);
	residue(const string& _aaType, const bool _hFlag, const bool _hPFlag);
	residue(const UInt _itsType, const bool _hFlag, const bool _hPFlag);
	~residue();

private:
	void initializeAtomsAndConnectivity();
	void initializeSidechainDihedralAngles();
	void initializePolarHDihedralAngle();
	// generateAtoms and buildConnectivity are auxilary functions
	// for constructors, so PRIVATE
	void generateAtoms();
	void buildConnectivity(); // defunct
	void buildConnectivityNew();

// ***********************************************************************
// ***********************************************************************
//	Mutation
// ***********************************************************************
// ***********************************************************************

private:
	// mutate should only be called from the 'chain' level object
	// So it is defined to be PRIVATE
	residue* mutate(const UInt _aaType);
	residue* mutateNew(const UInt _newTypeIndex);
	residue* fixBrokenResidue();

public:
	// Atoms Related
	void addAtom(PDBAtomRecord& _theRecord);
	void addAtom(pdbAtom& _pdbAtom);
	void makeAtomSilent(const UInt _atomIndex);
	void makeResidueSilent();

	// Residue Types
	bool isL(unsigned int resType){if(resType >= 0 && resType < 26){return true;}else{return false;}}        // L amino acids
	bool isG(unsigned int resType){if(resType == 26){return true;}else{return false;}}                       // Glycine
	bool isD(unsigned int resType){if(resType > 26 && resType < 53){return true;}else{return false;}}        // D amino acids
	bool isLNterm(unsigned int resType){if(resType >= 53 && resType < 79){return true;}else{return false;}}  // N terminal L amino acids
	bool isGNterm(unsigned int resType){if(resType == 79){return true;}else{return false;}}                  // N terminal Glycine
	bool isDNterm(unsigned int resType){if(resType > 79 && resType < 106){return true;}else{return false;}}  // N terminal D amino acids
	bool isLCterm(unsigned int resType){if(resType >= 106 && resType < 132){return true;}else{return false;}}// C terminal L amino acids
	bool isGCterm(unsigned int resType){if(resType == 132){return true;}else{return false;}}                 // C terminal Glycine
	bool isDCterm(unsigned int resType){if(resType > 132 && resType < 159){return true;}else{return false;}} // C terminal D amino acids
	bool isCofactor(unsigned int resType){if(resType >= 159){return true;}else{return false;}}               // Cofactor
	bool isCofactor(){if(itsType >= 159){return true;}else{return false;}}   

// ***********************************************************************
// ***********************************************************************
//	Dihedral Angle Calculation
// ***********************************************************************
// ***********************************************************************

	static UInt getNumBpt(const UInt _index);
	static UInt getNumDihedralAngles(const UInt _index, const UInt _bpt);
	double calculateDihedral(const vector<UInt>& _quad) const;
	double calculateDihedral(vector<atom*>& _quad) const;
	void calculateSidechainDihedralAngles();
	void calculatePolarHDihedralAngle();
	vector< vector< double > > getSidechainDihedralAngles();
	vector< vector< double > > randContinuousSidechainConformation();
	double getPhi();
	double getPsi();
	vector <double> getBackboneAngles();
	double getAngle(UInt angleType);
	double getOmega();
	double getAmide();
	double getAtomCharge(UInt _atomNum) {return residueTemplate::itsAmberElec.getItsCharge(itsType, itsAtoms[_atomNum]->itsType); }
	string getAtomName(UInt _atomNum) {return itsAtoms[_atomNum]->getName(); }
	void setOmega(double _omega);
	int setPhi(double _phi);
	int setPsi(double _psi);
	int setDihedral(double _dihedral, UInt _angleType, UInt _direction);
	void setRotamer(const UInt _lib, const UInt _bpt, const UInt _rotamer);
	void setRotamer(const UInt _bpt, const DouVec _chis);
	void setRotamerWithCheck(const UInt _lib, const UInt _bpt, const UInt _rotamer);
	void setPolarHRotamer(const UInt _rotamerIndex);
	void setPolarHRotamerWithCheck(const UInt _rotamerIndex);
	DouVec setRotamerWithCheckTest(const UInt _lib, const UInt _bpt, const UInt _rotamer);
	void setChiByDelta(const UInt _bpt, const UInt _index, const double _angleDelta);
	void setChi(const UInt _bpt, const UInt _index, const double _angle);
	void setChi(const UInt _index, const double _angle);
	void setBetaChi(const double _angle);
	void setPolarHChi(const UInt _rotamerIndex);
	void setPolarHChiByDelta(const UInt _atom1, const UInt _atom2, const double _angle);
	double getChi(const UInt _bpt, const UInt _index) const;
	double getChi(const UInt _index) const;
    double getBetaChi();
	double getBetaChiR();
    double getPolarHChi() const;
    double netCharge();

public:
	atom* getMainChain(UInt _index);
	bool inCube(const residue* _other, double _cutoff);
	bool isBonded(atom* _pAtom1, atom* _pAtom2);
	bool isBonded(UInt _index1, UInt _index2);
	bool isSeparatedByFewBonds(UInt _index1, UInt _index2);
    UInt getBondSeparation(UInt _index1, UInt _index2);
	bool isSeparatedByFewBonds(residue* _pRes1,UInt _index1,
		       residue* _pRes2, UInt _index2);
    UInt getBondSeparation(residue* _pRes1,UInt _index1,
               residue* _pRes2, UInt _index2);
	atom* getAtom(UInt _index)
		{	if (_index < itsAtoms.size())
			{	return itsAtoms[_index];
			}
			else
			{	return 0;
			}
		}
	void deleteAtom(const UInt atomIndex);
private:
	
	bool isSeparatedByOneOrTwoBonds(UInt _index1, UInt _index2);
	bool isSeparatedByOneOrTwoBonds(UInt _index1, residue* _pRes2, UInt _index2);
	bool isSeparatedByOneOrTwoBackboneBonds(UInt _index1, residue* _pRes2, UInt _index2);
	bool isSeparatedByThreeBackboneBonds(UInt _index1, residue* _pRes2, UInt _index2);
	bool isClash(UInt _index1, UInt _index2);
	bool isClash(UInt _index1, residue* _other, UInt _index2);
	UInt getNumHardClashes(residue* _other);
	UInt getNumHardClashes();
	UInt getNumHardBackboneClashes(residue* _other);


	// Total Residues
	void queryChildren(const UInt start);

	// Chain Related
	void setNextRes(residue* _pnr) {pItsNextRes = _pnr;}
	void setPrevRes(residue* _ppr) {pItsPrevRes = _ppr;}
public:
	residue* getNextRes() {return pItsNextRes;}
	residue* getPrevRes() {return pItsPrevRes;}

// ***********************************************************************
// ***********************************************************************
// 	Geometry Related
// ***********************************************************************
// ***********************************************************************

	dblVec getCoords( const string _atomName);
    dblVec getCoords( const UInt _atomIndex){return itsAtoms[_atomIndex]->getCoords();}
	void setCoords(UInt _atomIndex, dblVec _coords){return itsAtoms[_atomIndex]->setCoords(_coords);}
	void rotate(const point& _point, const dblMat& _RMatrix);
	void rotate(UInt _first, UInt _second,const double _theta);
	void rotate(atom* _pAtom1, atom* _pAtom2, const double _theta, bool backboneRotation);
	void rotateDihedral(atom* _pAtom1, atom* _pAtom2, double _deltaTheta,  UInt _angleType, UInt _direction);
	void rotate(atom* _pAtom1, const dblVec& _R_axis, const double _theta);
	void rotate(const point& _point, const dblVec& _R_axis, const double _theta);
	void rotate_new(atom* _pivotAtom, const dblMat& _RMatrix);
	void rotate_new(atom* _pivotAtom, atom* _firstAtom, const dblMat& _RMatrix);
	string getTypeStringFromAtomNum(UInt _atomNum) { return itsAtoms[_atomNum]->getType(); }
	string getNameStringFromAtomNum(UInt _atomNum) { return itsAtoms[_atomNum]->getName(); }
	void translate(dblVec* _pDoubleVector);
	void translate(const dblVec& _dblVec);
	void recursiveTranslateWithDirection(dblVec& _dblVec, UInt _direction);
	void recursiveTranslateLocal(dblVec& _dblVec, int direction);
	void recursiveTranslate(dblVec& _dblVec);
	void transform(dblMat* _pDoubleMatrix);
	void transform(const dblMat& _dblMat);
	void recursiveTransformR(dblMat& _dblMat);
	void recursiveTransform(dblMat& _dblMat);
	void recursiveTransformLocal(dblVec& atomCoords, double _deltaTheta, UInt _direction);
    void listConnectivity();
	void coilcoil(const double _pitch);
	residue* superimposeGLY();
	dblVec getBackBoneCentroid();
	void alignAmideProtonToBackbone();

	// utilities
	void printCoords() const;
	double intraEnergy();
	double intraSoluteEnergy();
	void polarizability();
	void polarizability(residue* _other);
	void updateMovedDependence(residue* _other, UInt _EorC);
	void calculateDielectrics();
    double calculateSolvationEnergy(UInt _atomIndex);
    double getSolvationEnergy();
    double getDielectric();
    double maxwellGarnettApproximation(UInt _atomIndex1, UInt _atomIndex2, double _dielectric, double _distanceSquared);
    double maxwellGarnettApproximation(UInt _atomIndex1, residue* _other, UInt _atomIndex2, double _dielectric, double _distanceSquared);
    double approximateDipoleDipolePolarization(UInt _atomIndex1, UInt _atomIndex2);
    double approximateDipoleDipolePolarization(UInt _atomIndex1, residue *_other, UInt _atomIndex2);
	double getIntraEnergy(const UInt atom1, residue* _other, const UInt atom2);
	double interEnergy(residue* _other);
	double interSoluteEnergy(residue* _other);
	double getSelfEnergy(residue* _other);
	double calculateHCA_O_hBondEnergy(residue* _other);
	double getVolume(UInt _method);
    bool notHydrogen(UInt _atomIndex);
    double getTotalVolumeofBondedAtoms(UInt _atomIndex);

private:
	double wodakVolume();

// ***********************************************************************
// ***********************************************************************
// surface area related operations
// ***********************************************************************
// ***********************************************************************


public:

	void initializeSpherePoints();
	void removeInterResidueSpherePoints(residue* _other);
	void removeIntraResidueSpherePoints();
	double tabulateSurfaceArea();
	double tabulateSurfaceArea(UInt _atomIndex);
	double tabulateSolvationEnergy(UInt _param);

// ***********************************************************************
// ***********************************************************************
//	Residue Type Base
// ***********************************************************************
// ***********************************************************************

public:
	static UInt getDataBaseSize();
	static string getDataBaseItem(const UInt _itemIndex);
	static UInt getAtomNameBaseSize(const UInt _restype);
	static string getAtomNameBaseItem(const UInt _restype,
		const UInt _atomtype);

private:
	void buildDataBase();
	void buildResidueDataBaseAminoAcids();
	void buildDihedralDataBase();
	void buildDataBaseFromPrep();
	void buildRotamerLib();
	void interpretBondingPattern();

	// intended non-const reference, to be modified
	void defineType(const string& _aaType);
	void setType(const string& _aaType);

public:
	string getType() const;
	string getType(UInt resType);
	UInt getTypeIndex() const {return itsType;}
	void printMainChain() const;
	void printBranchPoints() const;
	int getResNum() const {return itsResNum;}
    double getRadius(UInt atomIndex) {return itsAtoms[atomIndex]->getRadius();}
    static UInt getHowMany() {return howMany;}
	void setResNum(const UInt _num) {itsResNum = _num;}
	UInt getNumAtoms() const {return itsAtoms.size();}

// ***********************************************************************
// ***********************************************************************
//	flag Related
// ***********************************************************************
// ***********************************************************************

public:
	bool getHydrogensOn() const {return hydrogensOn;}
	void setHydrogensOn(const bool _hydrogensOn) ;
	bool getPolarHydorgensOn() const {return polarHydrogensOn;}
	void setPolarHydrogensOn(const bool _polarHydrogensOn);
	bool getHasPolarHRotamers() const {return dataBase[itsType].getHasPolarHRotamers(); }
	void setMoved(bool _moved, UInt _EorC);
	void setMoved();
	void setCheckMovedDependence (bool _check, UInt _EorC);
	void clearEnvironment();
	bool getMoved(UInt EorC); //Energy 0 or clashes 1
	bool getCheckMovedDependence(UInt _EorC);
	void setClashes (UInt _clashes);
	void sumClashes (UInt _clashes);
	UInt getClashes() const {return clashes;}
	void setBackboneClashes (UInt _clashes);
	void sumBackboneClashes (UInt _clashes);
	UInt getBackboneClashes() const {return clashesB;}
	void setEnergy (double _energy);
	void sumEnergy (double _energy);
	double getEnergy() const {return Energy;}
	void setResiduesPerTurnType(UInt _RPT);
	UInt getResiduesPerTurn() const {return RPT;}

// ***********************************************************************
// ***********************************************************************
// 	Static Variable Accessors
// ***********************************************************************
// ***********************************************************************

public:
	// static function cannot have const modifier
	static void setupDataBase();
	static void setupDataBase(const bool _Hflag);
	static void setupDataBase(const bool _Hflag, const bool _HPflag);
	static double getCutoffDistance() {return cutoffDistance; }
	static void setCutoffDistance( const double _cutoff ) { cutoffDistance = _cutoff; cutoffDistanceSquared = _cutoff*_cutoff; }
	static void setTemperature( const double _temp ) { temperature = _temp; }
	static double getTemperature() { return temperature; }
	static double getKT() { return KT; }
	static void setPolarizableElec( bool _polElec ) { polarizableElec = _polElec; }
	static void setElectroSolvationScaleFactor( const double _Esolv ) { EsolvationFactor = _Esolv; }
	static double getElectroSolvationScaleFactor() { return EsolvationFactor; }
	static void setHydroSolvationScaleFactor( const double _Hsolv ) { HsolvationFactor = _Hsolv; }
	static double getHydroSolvationScaleFactor() { return HsolvationFactor; }


// ***********************************************************************
// ***********************************************************************
//	Variables
// ***********************************************************************
// ***********************************************************************Cutoff

private:
	//variable declarations
	vector<atom*> itsAtoms;
	UInt itsType;
	int itsResNum;
	residue* pItsNextRes;
	residue* pItsPrevRes;
	bool hydrogensOn;
	bool polarHydrogensOn;
	bool isArtificiallyBuilt;
	bool movedE = true;
	bool movedC = true;
	bool movedB = true;
	bool dependentMoveE = false;
	bool dependentMoveC = false;
	bool dependentMoveB = false;
	UInt clashes = 0;
	UInt clashesB = 0;
	double Energy = 0.0;
	UInt RPT = 0;


	//variables relating to the rotameric state
	//or lack thereof....
	vector< vector <double> > itsSidechainDihedralAngles;
	double itsPolarHDihedralAngle;

	public:
	static vector<residueTemplate> dataBase;
	static void printDataBaseData()
		{	for (UInt i=0; i<dataBase.size(); i++)
			{	cout << "database residue type " << i << endl;
				dataBase[i].printAtomEnergyTypeDefinitions();
			}
		}
	private:
	static UInt howMany;
	static bool dataBaseBuilt;
	static bool polarizableElec;
	static double temperature;
	static double EsolvationFactor;
	static double HsolvationFactor;
	static double cutoffDistance;
	static double cutoffDistanceSquared;
	static double cutoffCubeVolume;
	static double dielectricWidth;
	static double KT;
};

#endif

// ***********************************************************************
// ***********************************************************************
// END OF FILE
// ***********************************************************************
// ***********************************************************************
