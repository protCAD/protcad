#include "assert.h"
#include "typedef.h"
#include "rotamer.h"
#include <iostream>
#include <fstream>

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef ROTAMERLIB_H
#define ROTAMERLIB_H

class rotamerLib
{	
public:
	friend class allowedResidue;
	rotamerLib(const UInt _numOfBpt);
	rotamerLib(const rotamerLib& _otherRotamerLib);  //deep copy
	~rotamerLib();

	void addRotamer(const StrVec& _strVec);
	void addRotamer(int _bpt, UIntVec _angles);
	
	bool rotamersExist(const UInt _bpt, const UInt _rotamer) const;
	DouVec getAngles(const UInt _bpt, const UInt _rotamer) const;
	double getEnergy(const UInt _bpt, const UInt _rotamer) const;

	void print() const;
	
	static UInt howMany;
	static UInt getHowMany() {return howMany;}

private:
	vector< vector<rotamer> > itsRotamers;
};

#endif
