
#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"

int main (int argc, char* argv[])
{
    string fileName = "apoBPB.pdb";
    molecule* pMol = pdbReader(fileName);
    protein* pProt = static_cast<protein*>(pMol);



	pProt->symmetryLinkChainAtoB(1,0);
	pProt->symmetryLinkChainAtoB(2,0);
	pProt->symmetryLinkChainAtoB(3,0);


	pProt->setSpaceLink(0, 0.5,0.5,0.5);	
	pProt->setSpaceLink(1, -0.5,0.5,-0.5);
	pProt->setSpaceLink(2, 0.5,-0.5,-0.5);
	pProt->setSpaceLink(3, -0.5,-0.5,0.5);

	pProt->translateWithLinkage(0, 10,10,10);
	pdbWriter(pProt, "out1.pdb");

	pProt->translateWithLinkage(0, -4,7,-15);
	pdbWriter(pProt, "out2.pdb");

	pProt->translateWithLinkage(0, -6,-17,5);
	pdbWriter(pProt, "out3.pdb");

	return 0;
}
