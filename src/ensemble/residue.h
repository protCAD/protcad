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

#ifndef LIGAND_H
#include "ligand.h"
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
	residue* fixBrokenResidue();

public:
	// Atoms Related
	void addAtom(PDBAtomRecord& _theRecord);
	void addAtom(pdbAtom& _pdbAtom);
	void makeAtomSilent(const UInt _atomIndex);
        
        //******************testing junk************************
	void accessMe();

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
	double getPhi();
	double getPsi();
	double getAngle(UInt angleType);
	double getOmega();
	double getAmide();
	double getAtomCharge(UInt _atomNum) {return residueTemplate::itsAmberElec.getItsCharge(itsType, itsAtoms[_atomNum]->itsType); }
	int setPhi(double _phi);
	int setPsi(double _psi);
	int setAngleLocal(double _angle, double deltaAngle, UInt angleType, int distance, int direction);
	int setDihedralLocal(double _deltaTheta, UInt _angleType);
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
        void setPolarHChi(const UInt _rotamerIndex);
	void setPolarHChiByDelta(const UInt _atom1, const UInt _atom2, const double _angle);
	double getChi(const UInt _bpt, const UInt _index) const;
	double getChi(const UInt _index) const;
        double getPolarHChi() const;

public:
	atom* getMainChain(UInt _index);
	bool inCube(const residue* _other, double _cutoff);
	bool isBonded(atom* _pAtom1, atom* _pAtom2);
	bool isBonded(UInt _index1, UInt _index2);
	bool isSeparatedByFewBonds(UInt _index1, UInt _index2);
	bool isSeparatedByFewBonds(residue* _pRes1,UInt _index1,
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
	bool isClash(UInt _index1, UInt _index2);
	bool isClash(UInt _index1, residue* _other, UInt _index2);
	UInt getNumHardClashes(residue* _other);
	UInt getNumHardClashes();


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
        dblVec getCoords( const UInt _atomIndex)
        {return itsAtoms[_atomIndex]->getCoords(); }
	void setCoords(UInt _atomIndex, dblVec _coords)
		{ return itsAtoms[_atomIndex]->setCoords(_coords);}
	void rotate(const point& _point, const dblMat& _RMatrix);
	void rotate(UInt _first, UInt _second,const double _theta);
	void rotate(atom* _pAtom1, atom* _pAtom2, const double _theta, bool backboneRotation);
	void rotateLocal(atom* _pAtom1, atom* _pAtom2, const double _theta, double deltaTheta, int distance, int direction);
	void rotateDihedralLocal(atom* _pAtom1, atom* _pAtom2, double _deltaTheta, UInt _direction);
	void rotateDihedral(atom* _pAtom1, atom* _pAtom2, double _deltaTheta, UInt _direction);
	void rotate(atom* _pAtom1, const dblVec& _R_axis, const double _theta);
	void rotate(const point& _point, const dblVec& _R_axis, const double _theta);
	void rotate_new(atom* _pivotAtom, const dblMat& _RMatrix);
	void rotate_new(atom* _pivotAtom, atom* _firstAtom, const dblMat& _RMatrix);
	string getTypeStringFromAtomNum(UInt _atomNum) { return itsAtoms[_atomNum]->getType(); }
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

	void coilcoil(const double _pitch);
	residue* superimposeGLY();
	dblVec getBackBoneCentroid();

	// utilities
	void printCoords() const;
	double intraEnergy();
	double intraSoluteEnergy();
	vector <double> calculateDielectric(residue* _other, UInt _atomIndex);
	vector <double> calculateDielectric(residue* _other, atom* _atom);
    vector <double> calculateSolvationEnergy(UInt _atomIndex);
	double getIntraEnergy(const UInt atom1, residue* _other, const UInt atom2);
	double interEnergy(residue* _other);
	double interSoluteEnergy(residue* _other);
     double interEnergy(ligand* _other);
	double getSelfEnergy(residue* _other);
	double calculateHCA_O_hBondEnergy(residue* _other);
	double BBEnergy();
	double BBEnergy(residue* _other);
	/* 1-4 atom interaction variables and accessors
	static double oneFourVDWScaleFactor;
	static double oneFourAmberElecScaleFactor;
	static void setOneFourVDWScaleFactor(double scale) { oneFourVDWScaleFactor = scale; }
	static double getOneFourVDWScaleFactor() { return oneFourVDWScaleFactor; }
	static void setOneFourAmberElecScaleFactor(double scale) {oneFourAmberElecScaleFactor = scale; }
	static double getOneFourAmberElecScaleFactor() { return oneFourAmberElecScaleFactor; }
	*/
	double getVolume(UInt _method);

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
	void buildDataBaseAA();
	void buildResidueDataBaseAminoAcids();
	//void buildDataBaseAE();
	void buildDataBaseCC();
	void buildRotamerLib();
	void interpretBondingPattern();

	// intended non-const reference, to be modified
	void defineType(const string& _aaType);
	void setType(const string& _aaType);

public:
	string getType() const;
	UInt getTypeIndex() const {return itsType;}
	void printMainChain() const;
	void printBranchPoints() const;
	int getResNum() const {return itsResNum;}
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
	static UInt getHowMany() {return howMany;}
	static double getCutoffDistance() {return cutoffDistance; }
	static void setCutoffDistance( const double _cutoff ) { cutoffDistance = _cutoff; cutoffDistanceSquared = _cutoff*_cutoff; }

// ***********************************************************************
// ***********************************************************************
//	Variables
// ***********************************************************************
// ***********************************************************************

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
	static double cutoffDistance;
	static double cutoffDistanceSquared;
};

#endif

// ***********************************************************************
// ***********************************************************************
// END OF FILE
// ***********************************************************************
// ***********************************************************************
