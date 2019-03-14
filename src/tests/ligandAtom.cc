//
// Start ligandAtom code
//

#include "ligandAtom.h"

UInt ligandAtom::howManyLigAtoms=0;

ligandAtom::ligandAtom() : atom::atom()
{
    //ligAtomInitialize();
     itsAmberAllType="";
    itsAmberUnitedType="";
    itsAmberAllCharge=0.0;
    itsAmberUnitedCharge=0.0;
    itsSummaAtomType="";
    itsSummaEnvType="";
    itsSixAtomSolvationType="";
    
    howManyLigAtoms++;
    
}

ligandAtom::ligandAtom(bool _hetflag, const PDBAtomRecord& _theRecord) : atom::atom(_theRecord, _hetflag)
{
        //ligandAtom-specific stuff...
   // ligAtomInitialize();
    itsAmberAllType="";
    itsAmberUnitedType="";
    itsAmberAllCharge=0.0;
    itsAmberUnitedCharge=0.0;
    itsSummaAtomType="";
    itsSummaEnvType="";
    itsSixAtomSolvationType="";
    
    howManyLigAtoms++;
}

// deep-copy constructor
ligandAtom::ligandAtom(const ligandAtom& _rhs) : atom::atom(_rhs)
{
    
    // ligandAtom stuff...
        itsAmberAllType=_rhs.itsAmberAllType;
        itsAmberUnitedType=_rhs.itsAmberUnitedType;
        itsAmberAllCharge=_rhs.itsAmberAllCharge;
        itsAmberUnitedCharge=_rhs.itsAmberUnitedCharge;
        itsSummaAtomType=_rhs.itsSummaAtomType;
        itsSummaEnvType=_rhs.itsSummaEnvType;
        itsSixAtomSolvationType=_rhs.itsSixAtomSolvationType;
        
        howManyLigAtoms++;
        
    //cout << "Done with deep copy constructor" << endl;
}

void ligandAtom::ligAtomInitialize()
{
    itsAmberAllType="";
    itsAmberUnitedType="";
    itsAmberAllCharge=0.0;
    itsAmberUnitedCharge=0.0;
    itsSummaAtomType="";
    itsSummaEnvType="";
    itsSixAtomSolvationType="";
}
