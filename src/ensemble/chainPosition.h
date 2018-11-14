// filename: chainPosition.h

#include "assert.h"
#include <vector>
#include "typedef.h"
#include "ran.h"
#include "residue.h"
#include "allowedResidue.h"
#include "rotamerLib.h"

#ifndef CHAINPOSITION_H
#define CHAINPOSITION_H
//#warning "chain position read in"
class chainPosition
{
public:
	friend class chain;

	// Constructor and Destructor declaration
	chainPosition(UInt _index, UInt _aatype);
	chainPosition(const chainPosition& _rhs);
	~chainPosition();

	// Accessors
	void setResNotAllowed(UInt _aaType);
	void setResAllowed(UInt _aaType);
	void setOnlyAllowedIdentity(UInt _aaType);
	void setRotamerNotAllowed(UInt _aaType, UInt _bpt, UInt _rotamer);
	void setRotamerNotAllowed(UInt _aaType, UInt _rotamer);
	UIntVec getAllowedRotamers(UInt _resType, UInt _bpt);
    vector <UIntVec> getAllowedRotamers(UInt _resType);
	bool residueIsAllowed(UInt _aaType);
	UInt getResNum() {return itsResNum;}
	UInt getNumAllowedRes() {return itsAllowedResidues.size();}
	UInt getCurrentAllowedResIdentity() const;
	UIntVec getResAllowed();

	int chooseResidueFromLibrary(ran& _ran);
	int chooseRotamerFromLibrary(UInt _bpt, ran& _ran);
	int chooseBranchpoint(ran& _ran);
	
	vector<UInt> getCurrentRotamerIndex() const;
	void setCurrentRotamerIndex(const UInt _index);
	void setCurrentRotamerIndex(vector<UInt> _rotamer);
	void setCurrentPolarHRotamerIndex(const UInt _rotamer);
	void listAllowedRotamers();

//	void listAllowedConformations( UInt _aaType );  CONFORMATION CODE NOT USED
//	void setConformationNotAllowed ( UInt _aaType, vector <UInt> _conformation);


	void  setRotamerLibIndex(const UInt _index);
	UInt  getRotamerLibIndex() const {return itsRotamerLibIndex;}

	static UInt getHowMany() {return howMany;}

	void setSecondaryStructureIndex(UInt _index) {itsSecondaryStructureIndex = _index;}
	UInt getSecondaryStructureIndex() {return itsSecondaryStructureIndex;}
	int getCurrentResIndex() {return itsCurrentAllowedResIndex;}

	/*this function returns a three dimensional vector with
	 the following properties:  the first dimension is the set of allowed residues,
	 with each UInt being the index of an amino acid.  the second dimension, if NULL,
	 indicates that there are no occurrences of that amino acid.  if not null, then 
	 the 0th element of that vector is a vector of UInts which indicate allowed 
	 rotamer conformations for that amino acid. This is called by the frontend
	 graphics routines */

	vector < vector <UIntVec> > getAllowedDB() const;

private:

	void setCurrentResIndex(const UInt _aaType);
	void removeResidue(UInt _index);
	void addResidue(UInt _index);

	//variable declarations
	UInt itsResNum;
	vector<allowedResidue> itsAllowedResidues;
	int itsCurrentAllowedResIndex;
	static UInt howMany;
	UInt itsSecondaryStructureIndex;
	UInt itsRotamerLibIndex;
};

#endif
