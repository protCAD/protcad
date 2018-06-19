//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protOpt        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
	if (argc !=3)
	{   cout << "protOpt <inFile.pdb> <outFile.pdb>" << endl;
		exit(1); }

	string infile = argv[1];
	string outFile = argv[2];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
    bool homosymmetric = true;
    bool backbone = false;

    UInt _frozenResidues[] = {216};
	UInt _activeChains[] = {0};
	UInt activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]), frozenResiduesSize = sizeof(_frozenResidues)/sizeof(_frozenResidues[0]);
	UIntVec activeChains, frozenResidues;
	for (UInt i = 0; i < activeChainsSize; i++)
	{
		activeChains.push_back(_activeChains[i]);
	}
	for (UInt i = 0; i < frozenResiduesSize; i++)
	{
		frozenResidues.push_back(_frozenResidues[i]);
    }
    if (homosymmetric)
    {
        _prot->symmetryLinkChainAtoB(activeChains[0],1);
    }

	double startEnergy = _prot->protEnergy();
	time_t start,end;
	time (&start);
    //_prot->protOpt(backbone);
    _prot->protOpt(backbone,frozenResidues,activeChains);
	time (&end);
	double endEnergy = _prot->protEnergy();
	cout << _prot->protEnergy() << " " << (endEnergy-startEnergy)/difftime(end,start) << " ";
	pdbWriter(_prot, outFile);

	return 0;
}

