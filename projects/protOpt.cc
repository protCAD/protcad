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
    bool homosymmetric = false;
    bool backbone = true;
    clock_t start, end;
	double cpu_time_used;
    
    residue::setElectroSolvationScaleFactor(1.0);
    residue::setHydroSolvationScaleFactor(1.0);
    amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);

    UInt _frozenResidues[] = {};
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
        for (UInt i = 0; i < _prot->getNumChains(); i++)
        {_prot->symmetryLinkChainAtoB(activeChains[0],i);}
    }
    cout << "start Energy: " << _prot->protEnergy() << endl;
    //_prot->protOpt(backbone,frozenResidues,activeChains);
    start = clock();
    _prot->protOpt(backbone);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    _prot->setMoved(true);
    cout << "end Energy: "  << _prot->protEnergy() << " time: " << cpu_time_used << endl;
	pdbWriter(_prot, outFile);

	return 0;
}

