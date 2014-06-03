// test conformation data structure  3/22/01

#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"

int main (void)
{
	string fileName = "1HOL.pdb";
	molecule* pMol = pdbReader(fileName);
	protein* pProt = static_cast<protein*>(pMol);

	pProt->symmetryLinkChainAtoB(1,0);
	
	pProt->activateForRepacking(0,pProt->getIndexFromResNum(0,8));

    vector <UInt> bannedResList;
    bannedResList.resize(0);

    // keep ALA
    bannedResList.push_back(1); // ban ARG
    bannedResList.push_back(2); // ban ASN
    bannedResList.push_back(3); // ban ASP
    bannedResList.push_back(4); // ban CYS
    bannedResList.push_back(5); // ban GLN
    bannedResList.push_back(6); // ban GLU
    bannedResList.push_back(7); // ban GLY
    bannedResList.push_back(8); // ban HIS
    // keep ILE
    // keep LEU
    bannedResList.push_back(11); // ban LYS 
    bannedResList.push_back(12); // ban MET
    bannedResList.push_back(13); // ban PHE
    bannedResList.push_back(14); // ban PRO
    bannedResList.push_back(15); // ban SER
    // keep THR
    bannedResList.push_back(17); // ban TRP
    bannedResList.push_back(18); // ban TYR
    // keep VAL

	pProt->setListNotAllowed(0, pProt->getIndexFromResNum(0,8), bannedResList);

	cout << "Listing allowed conformations for position (0,8)\n";
	pProt->listAllowedConformations(0, pProt->getIndexFromResNum(0,8));
	cout << "Listing allowed conformations for position (1,8)\n";
	pProt->listAllowedConformations(1, pProt->getIndexFromResNum(1,8));

	cout << "Setting canonical helix rotamers for position (0,8)\n";
	pProt->setCanonicalHelixRotamersOnly(0, pProt->getIndexFromResNum(0,8));

	cout << "Listing allowed conformations for position (0,8)\n";
	pProt->listAllowedConformations(0, pProt->getIndexFromResNum(0,8));
	cout << "Listing allowed conformations for position (1,8)\n";
	pProt->listAllowedConformations(1, pProt->getIndexFromResNum(1,8));

	return 0;
}	
