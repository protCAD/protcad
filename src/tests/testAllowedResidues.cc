#include "pdbReader.h"
#include "atomIterator.h"
#include "ran.h"

int main(void)
{
    string fileName = "1UBI.pdb";
    molecule* pTheProtein = pdbReader(fileName);
    if ( pTheProtein == 0) return 1;   // exit if null

    UIntVec allowedResidues;

    static_cast<protein*>(pTheProtein)->activateForRepacking(0,2);

    allowedResidues = static_cast<protein*>(pTheProtein)->getResAllowed(0,2);
    for (UInt i = 0; i < allowedResidues.size(); i ++)
    {
        cout << residue::getDataBaseItem(allowedResidues[i]) << " ";
    }
    cout << endl;


    static_cast<protein*>(pTheProtein)->setOnlyROCHydrophobic(0,2);
    cout << "Set to hydrophobic only" << endl;

    allowedResidues = static_cast<protein*>(pTheProtein)->getResAllowed(0,2);
    for (UInt i = 0; i < allowedResidues.size(); i ++)
    {
        cout << residue::getDataBaseItem(allowedResidues[i]) << " ";
    }
    cout << endl;
}

