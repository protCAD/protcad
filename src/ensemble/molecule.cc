// filename: molecule.cc
// contents: class molecule implementation

#include "molecule.h"

void molecule::initialize()
{
}

molecule::molecule()
{
#ifdef MOLECULE_DEBUG
	cout<< "default molecule constructor called" << endl;
#endif
	initialize();
	setMoleculeType(0);
	itsName = "UNK";
}

molecule::molecule(const string& _name)
{	initialize();
	setMoleculeType(0);
	itsName = _name;
}

molecule::molecule(const molecule& _rhs)
{	setMoleculeType(_rhs.itsMolType);
	itsName = _rhs.itsName;
}

molecule::~molecule()
{
#ifdef MOLECULE_DEBUG
		cout<< "molecule destructor called " << endl;
#endif
}
