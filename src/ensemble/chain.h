// filename: chain.h

#include "assert.h"
#include <string.h>
#include <vector>
#include "typedef.h"
#include "ran.h"
#include "enums.h"
//#include "CMath.h"
#include "chainPosition.h"
#include "residue.h"
#include "chainModBuffer.h"

#ifndef ATOMITERATOR_H
class atomIterator;
#endif

#ifndef CHAIN_H
#define CHAIN_H
//#warning "chain.h read in"
#include "secondaryStructure.h"

class chain
{
public:
	friend class atomIterator;
	friend class ruler;

	// Constructor and Destructor declaration
	chain();
	chain(char _id);
	chain(const chain& _rhs);
	~chain();
private:
	void initialize();

public:
	void add(residue* _pResidue);

	// Accessors
	bool isArtificiallyBuilt(UInt j) { return itsResidues[j]->isArtificiallyBuilt; }
    dblVec getCoords(const UInt _resIndex, const string _atomName)
        {return itsResidues[_resIndex]->getCoords(_atomName);}
    dblVec getCoords(const UInt _resIndex, const UInt _atomIndex)
        {return itsResidues[_resIndex]->getCoords(_atomIndex);}
    void setCoords(UInt _resIndex, UInt _atomIndex, dblVec _coords)
	{ return itsResidues[_resIndex]->setCoords(_atomIndex, _coords);}
    UInt getNumAtoms(const UInt _resIndex)
        {return itsResidues[_resIndex]->getNumAtoms();}
    UInt getNumResidues() {return itsResidues.size();}
	int mapResNumToChainPosition(const int _resNum);
private:
	void activateChainPosition(const UInt _index);
	void activateAllChainPositions();
public:
	bool activateForRepacking(const UInt _index);
	bool activateForRepacking(const UInt _start, const UInt _end);
	bool activateAllForRepacking();
	static UInt getHowMany() {return howMany; }
	void setChainID(char _id) {itsChainID = _id;}
	int getResNum(UInt _resIndex) {return itsResidues[_resIndex]->getResNum();}
	char getChainID() {return itsChainID;}
	residue* getResidue(UInt _resIndex) { return itsResidues[_resIndex]; }
	UInt getTypeFromResNum(UInt _resNum) { return itsResidues[_resNum]->getTypeIndex(); }
    double getRadius(UInt resIndex, UInt atomIndex) {return itsResidues[resIndex]->getRadius(atomIndex);}
	string getTypeStringFromAtomNum(UInt _resNum, UInt _atomNum) { return itsResidues[_resNum]->getTypeStringFromAtomNum( _atomNum); }
	string getTypeStringFromResNum(UInt _resNum) {return itsResidues[_resNum]->getType();}
	double getAtomCharge(UInt _resNum, UInt _atomNum) { return itsResidues[_resNum]->getAtomCharge(_atomNum); }
	void mutate(const UInt _indexInChain, const UInt _aaType);
	void mutateWithoutBuffering(const UInt _indexInChain, const UInt _aaType);
    void fixBrokenResidue(const UInt _indexInChain);
    void rebuildResidue(const UInt _indexInChain);
	void rebuildResiduesInChain();
	bool isDAminoAcid(residue* currentRes);
	void redoModification(chainModBuffer _redoBuffer);
	void makeAtomSilent(const UInt _resIndex, const UInt _atomIndex);
	void makeResidueSilent(const UInt _resIndex);
	int getLastModificationPosition() { return itsLastTargetResidue; }
	void randomizeSystem(ran& _ran);
	void makeAllAlanine();
    void removeResidue(UInt _resNum);
	bool isCofactor(UInt resIndex){return itsResidues[resIndex]->isCofactor();}   


	// single modification buffers
	vector<chainModBuffer> performRandomMutation(ran& _ran);
	vector<chainModBuffer> performRandomRotamerChange(ran& _ran);
	vector<chainModBuffer> performRandomRotamerRotation(ran& _ran);
	vector<chainModBuffer> performRandomMutation(ran& _ran, vector <int> _position);
	vector<chainModBuffer> performRandomRotamerChange(ran& _ran, vector <int> _position);
	vector<chainModBuffer> performRandomRotamerRotation(ran& _ran, vector <int> _position);

	// multiple modification buffers
	vector<chainModBuffer> saveCurrentState();
	vector<chainModBuffer> saveCurrentState(vector <int> _position);
	void commitState();
	void undoState();
	void copyState(vector <chainModBuffer> _externalBuffer);

	void commitLastMutation();
	void commitLastRotamerChange();
	void commitLastRotamerRotation();
	void undoLastMutation();
	void undoLastRotamerChange();
	void undoLastRotamerRotation();
	void repeatModification(const chainModBuffer& _redoBuffer);

	void finishChainBuild();
	void assignSecondaryStructure();

	void listSecondaryStructure() const;
	void listDihedrals();
	void listAllowedRotamers(UInt _indexInChain) const;
	void listChiDefinitions() const;
	UInt getNumChis(const UInt _resIndex, const UInt _bpt);
	double netCharge();
	vector <dblVec> saveCoords(UInt resIndex);
	void setAllCoords(UInt resIndex, vector<dblVec> allCoords);
	
	void setMoved (UInt resIndex, bool _moved, UInt _EorC) {itsResidues[resIndex]->setMoved(_moved, _EorC);}
	void setMoved (bool _moved, UInt _EorC);
	bool getMoved (UInt resIndex, UInt _EorC) {return itsResidues[resIndex]->getMoved(_EorC);}
	double getEnergy (UInt resIndex) {return itsResidues[resIndex]->getEnergy();}
	double getSolvationEnergy(const UInt _resIndex) {return itsResidues[_resIndex]->getSolvationEnergy();}
	double getDielectric(const UInt _resIndex) {return itsResidues[_resIndex]->getDielectric();}
	double getBetaChi(const UInt _resIndex) {return itsResidues[_resIndex]->getBetaChi();}
	void setBetaChi(const UInt _resIndex, double _chi) {return itsResidues[_resIndex]->setBetaChi(_chi);}
	UIntVec getActiveResidues() { return itsRepackActivePositionMap;}
	void setRotamerNotAllowed (const UInt _indexInChain, const UInt aaType, const UInt _bpt, const UInt _rotamer);
	UIntVec getAllowedRotamers ( const UInt _indexInChain, const UInt  _aaType, const UInt _bpt);
	vector <UIntVec> getAllowedRotamers ( const UInt _indexInChain, const UInt  _aaType);
	void setResNotAllowed(const UInt _indexInChain, const UInt _aaType);
	void setResAllowed(const UInt _indexInChain, const UInt _aaType);
	UIntVec getResAllowed (const UInt _indexInChain);
	void setOnlyNativeIdentity(const UInt _indexInChain);
	vector< vector< double > > randContinuousSidechainConformation(UInt _resIndex) {return itsResidues[_resIndex]->randContinuousSidechainConformation();}
	void setRotamerWithoutBuffering(const UInt _indexInChain, vector<UInt> _rotamer);
	void setRotamerWithoutBuffering(const UInt _indexInChain, vector<UInt> _rotamer, vector< vector< double> > dihedralAngles);
	void setRotamerWithoutBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _rotamer);
	void setRotamerWithBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _rotamer);
	void setRotamerWithBuffering(const UInt _indexInChain, vector<UInt> _rotamer);
	void setPolarHRotamerWithoutBuffering(const UInt _indexInChain, const UInt _rotamerIndex);
	void setAllHydrogensOn(const bool _hydrogensOn);
	void setAllPolarHydrogensOn(const bool _polarHydrogensOn);
	void setRelativeChi(const UInt _resIndex, const UInt _bpt, const UInt _chi, const double _angle);
	UIntVec getCurrentRotamer(UInt _resIndex) { return itsChainPositions[_resIndex]->getCurrentRotamerIndex(); }
	double getChi(const UInt _resIndex, const UInt _bpt, const UInt _chi)
		{ return itsResidues[_resIndex]->getChi(_bpt, _chi); }
	void setChi(const UInt _resIndex, const UInt _bpt, const UInt _chi, const double _angle)
		{ itsResidues[_resIndex]->setChi(_bpt, _chi, _angle); }
	void setDihedrals(UInt _resIndex, UInt _bpt, vector <double> _dihedrals);
	vector < vector <double> >  getDihedrals(UInt _resIndex) {return itsResidues[_resIndex]->getSidechainDihedralAngles() ;}
	void setDihedralWithBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _dihedralIndex, const double _dihedralAngle);
	void setDihedralWithoutBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _dihedralIndex, const double _dihedralAngle);
	vector< vector< double> > getSidechainDihedralAngles(UInt _indexInChain);
	void setSidechainDihedralAngles(UInt _indexInChain, vector< vector<double> > Angles);
	int chooseNextTargetPosition(ran& _ran);
	UInt chooseNextMutationIdentity(ran& _ran, vector <int> _position) { return chooseTargetIdentity(_position[2],_ran); }
	residue* superimposeGLY(const UInt _residue);
	double calculateHCA_O_hBondEnergy(chain* _other);
	UInt getBackboneClashes();
	void updateBackboneClashes(chain* _other);
	void updateBackboneClashes();
	UInt getClashes();
	void updateClashes(chain* _other);
	void updateClashes();
	UInt getClashes (UInt resIndex) {return itsResidues[resIndex]->getClashes();}
	void listConnectivity(UInt _resIndex) {return itsResidues[_resIndex]->listConnectivity();}
private:
	int chooseTargetResidue(ran& _ran);
	int chooseTargetIdentity(const UInt _indexInChain, ran& _ran);
	int chooseTargetBranchpoint(const UInt _indexInChain, ran& _ran);
	int chooseTargetRotamer(const UInt _indexInChain, const UInt _bpt, ran& _ran);
	int chooseTargetDihedralIndex(const UInt _indexInChain, const UInt _bpt, ran& _ran);
	double chooseTargetDihedralAngle(const UInt _indes, const UInt _bpt, const UInt _angleIndex, ran& _ran);

	void addResidue(residue* _pResidue);
	void addChainPosition(chainPosition* _pChainPosition);
	void addSecondaryStructure(secondaryStructure* _pSecondaryStructure);

public:

	dblVec getSpaceLink() { return itsSpaceLink; }
	void setSpaceLink(const dblVec& _relation) { itsSpaceLink = _relation; }
	dblVec getSpinLink() {return itsSpinLink; }
	void setSpinLink(const dblVec& _relation) {itsSpaceLink = _relation; }

	void translate(const dblVec& _dblVec);
	void translate(const double _x, const double _y, const double _z);
	void translateR(const double _x, const double _y, const double _z);
	void transform(const dblMat& _dblMat);
	void coilcoil(const double _pitch);
	dblVec getBackBoneCentroid();

	void rotate(const axis _axis,const double _theta);
	void rotateRelative(const axis _axis,const double _theta);
	void rotate(const point& _point, const dblVec& _R_axis, const double _theta);
	double getResiduesPerTurn(const UInt _resIndex);
	double getResiduesPerTurn(double phi, double psi);
	UInt getBackboneSequenceType(const UInt _resIndex);
	UInt getBackboneSequenceType(double RPT, double phi);
	void updateResiduesPerTurnType();
	double getPhi(const UInt _indexInChain);
	double getPsi(const UInt _indexInChain);
	double getAngle(const UInt _indexInChain, UInt angleType);
	double getAmide(const UInt _indexInChain);
	int setPhi(const UInt _indexInChain, double _phi);
	int setPsi(const UInt _indexInChain, double _psi);
	int setDihedral(const UInt _resIndex, double _dihedral, UInt _angleType, UInt _direction);
	double getDielectric(UInt _resIndex, UInt _atomIndex) {return itsResidues[_resIndex]->itsAtoms[_atomIndex]->getDielectric();}
	double intraEnergy();
	void updateEnergy();
	void updateEnergy(chain* _other);
	void updateMovedDependence(UInt _EorC);
	void updateMovedDependence(chain* _other, UInt _EorC);
	double getEnergy();
	void polarizability();
	void polarizability(chain* _other);
	void calculateDielectrics();
	double interEnergy(chain* _other);
	double getInterEnergy(const UInt _res1, chain* _other, const UInt _res2);
	double getInterEnergy(const UInt _residue1, const UInt _atom1, chain* _other, const UInt _residue2, const UInt _atom2);
	double getSelfEnergy(UInt _residueIndex);

	double getPositionIntraEnergy(vector<int> _position);
	double getPositionInterEnergy(vector<int> _position, chain* _other);
	double getPositionIntraSoluteEnergy(vector<int> _position);
	double getPositionIntraSoluteEnergy(UInt _residueIndex);
	double getPositionInterSoluteEnergy(vector<int> _position, chain* _other);

	double getVolume(UInt _method);

	double rotamerEnergy();
	double rotamerEnergy(const UInt _index);

	vector<chainPosition*> getChainPositionVector() {return itsChainPositions;}
	vector<UInt> getRepackActivePositionMap() {return itsRepackActivePositionMap;}

	void symmetryLinkResidueAtoB(UInt _masterResIndex, UInt _slaveResIndex);

// surface area

public:

	void initializeSpherePoints();
	void initializeSpherePoints(UInt _residue);
	double tabulateSurfaceArea();
	double tabulateSurfaceArea(UInt _residue);
	double tabulateSurfaceArea(UInt _residueIndex, UInt _atomIndex);
	double tabulateSolvationEnergy(UInt _param);
	double tabulateSolvationEnergy(UInt _residue, UInt _param);
	void removeIntraChainSpherePoints();
	void removeIntraChainSpherePoints(UInt _residue);
	void removeInterChainSpherePoints(chain* _other);
	void removeInterChainSpherePoints(UInt _residue, chain* _other);

private:
//	private transformation variables

	dblVec itsSpaceLink;
	dblVec itsSpinLink;
private:
//	Undo and buffering accessors
	void bufferResidueIntoUndoBuffer(UInt indexInChain);
	void resetUndoBuffer();
	void bufferResidueIntoRedoBuffer(UInt indexInChain);
	void resetRedoBuffer();

//  residue level symmetry linking for periodic boundary conditions, etc ...
    UIntVec itsIndependentPositions;
    vector <vector <int> > itsResidueLinkageMap;
public:
	//variable declarations
	static UInt howMany;
	vector<residue*> itsResidues;
	vector<chainPosition*> itsChainPositions;
	vector<UInt> itsRepackActivePositionMap;
	vector<secondaryStructure*> itsSecondaryStructures;
	char itsChainID;
	int itsLastTargetResidue;

private:
	// variables to buffer undo info
	vector<chainModBuffer> itsBuffers; // single point changes 0 -> undo buffer 1 -> redobuffer
	vector<chainModBuffer> stateBuffers; // array of buffers containing list of residues
};
#endif
