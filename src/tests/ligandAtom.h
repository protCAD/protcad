//
// LigandAtom.h
//

#ifndef ASSERT_H
#include "assert.h"
#endif

#include <string.h>
#include <vector.h>

#ifndef TYPEDEF_H
#include "typedef.h"
#endif

#ifndef GENERALIO_H
#include "generalio.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIGANDATOM_H
#define LIGANDATOM_H
class ligandAtom : private atom
{
public:
    
    // Constructor and Destructor declaration
    ligandAtom();
    ligandAtom(bool _hetflag, const PDBAtomRecord& _theRecord);
    ligandAtom(const ligandAtom& _rhs); //Deep copy constructor
    ~ligandAtom();
    
    // Initialize and setup
    void ligAtomInitialize(); //set everything to defaults   
    
    void setAmberUnited(string _type){itsAmberUnitedType=_type;}
    void setAmberAllCharge(double _charge){itsAmberAllCharge=_charge;}
    void setAmberUnitedCharge(double _charge){itsAmberUnitedCharge=_charge;}
    void setSummaAtom(string _type){itsSummaAtomType=_type;}
    void setSummaEnv(string _type){itsSummaEnvType=_type;}
    void setSixAtomSolvation(string _type){itsSixAtomSolvationType=_type;}
    
    // Accessor Functions
    string getAmberAll(){return itsAmberAllType;}
    string getAmberUnited(){return itsAmberUnitedType;}
    double getAmberAllCharge(){return itsAmberAllCharge;}
    double getAmberUnitedCharge(){return itsAmberUnitedCharge;}
    string getSummaAtom(){return itsSummaAtomType;}
    string getSummaEnv(){return itsSummaEnvType;}
    string getSixAtomSolvation(){return itsSixAtomSolvationType;}
    
private:
    string itsAmberAllType;
    string itsAmberUnitedType;
    double itsAmberAllCharge;
    double itsAmberUnitedCharge;
    string itsSummaAtomType;
    string itsSummaEnvType;
    string itsSixAtomSolvationType;
    
    static UInt howManyLigAtoms;

};
#endif
