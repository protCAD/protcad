// filename allowedResidue.h

#include "assert.h"
#include <string.h>
#include <vector>
#include "typedef.h"
#include "ran.h"
#include "residue.h"
#include "residueTemplate.h"

#ifndef ALLOWEDRESIDUE_H
#define ALLOWEDRESIDUE_H

class allowedResidue 
{
public:
	friend class chainPosition;
	friend class chain;

	allowedResidue();
	allowedResidue(const UInt _aaType);
	allowedResidue(const UInt _aaType, const UInt _lib);
	allowedResidue(const allowedResidue& _rhs);
	~allowedResidue();

	void buildAllowedRotamers(const UInt _lib);
	UInt getIdentity() const {return itsIdentity;}
	void setRotamerNotAllowed(const UInt _bpt, const UInt _rotamer);
	UInt getRotamerLibSize(const UInt _bpt) const;
	int chooseRotamerFromLibrary(const UInt _bpt, ran& _ran) const;
	int chooseBranchpoint(ran& _ran) const;
	void listAllowedRotamers() const;

	vector<UIntVec> getAllowedRotamers() const {return itsAllowedRotamers;}


	vector<UInt> getCurrentRotamerIndex() const { return itsCurrentRotamerIndex;}
	UInt getCurrentPolarHRotamerIndex() const { return itsCurrentPolarHRotamerIndex;}
	void setCurrentRotamerIndex(const UInt _rotamerIndex);
	void setCurrentRotamerIndex(vector<UInt> _rotamerIndex);
	void setCurrentPolarHRotamerIndex(const UInt _rotamerIndex);

private:
	int getRotamerFromIndex(const UInt _bpt, const UInt _index) const;
	void deleteRotamer(const UInt _bpt, const UInt _index);
	UInt itsIdentity;
	vector<UInt> itsCurrentRotamerIndex;
	vector<UIntVec> itsAllowedRotamers;
	UInt itsCurrentPolarHRotamerIndex;

/*   CONFORMATION CODE NOT USED
public:
	void buildAllowedConformations(UInt _rotamerModulus); 
	void listAllowedConformations();
	void setConformationNotAllowed(UIntVec _conformation);
private:
	UIntVec incrementConf(UIntVec _tmpVec, UInt _position, UInt _rotamerModulus);
	vector <UIntVec> itsAllowedConformations;
*/
};
#endif
