#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "annealer.h"



int main(int argc, char* argv[])
{
   residue::setCutoffDistance(6.0);
    pmf::setScaleFactor(0);
    rotamer::setScaleFactor(0);
    microEnvironment::setScaleFactor(0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.9);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0);

    string inputFileName = argv[1];
    string outputFileName= argv[2];
    enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    PDBInterface* thePDB = new PDBInterface(inputFileName);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* prot  = static_cast<protein*>(pMol);

	ran ranNumber;
	ranNumber.setSeed(132);


	for (UInt chain = 0; chain < prot->getNumChains(); chain ++)
	{
		for (UInt res = 0; res < prot->getNumResidues(chain); res ++)
		{	
			prot->activateForRepacking(chain,res);
			prot->mutate(chain,res,A);
		}
	}

	double oldEnergy = 10E5;
	for (UInt i = 0; i < 10000; i ++)
	{
		UInt chain = (UInt)(ranNumber.getNext()*prot->getNumChains());
		if (chain == prot->getNumChains()) chain = chain -1;
		UInt res = (UInt)(ranNumber.getNext()*prot->getNumResidues(chain));
		if (res == prot->getNumResidues(chain)) res = res -1;
		UInt choice = (UInt)(ranNumber.getNext()*20);
		prot->saveCurrentState();
		
		double beta = 1E-4 - i*(1E-4 - 1)/10000;
		cout << "mutating " << chain << " " << res << " to " << choice << endl;
		prot->mutate(chain, res, choice);
		prot->optimizeRotamers();
		double energy = prot->intraEnergy();
		if (energy < oldEnergy)
		{
			cout << "ACCEPTED" << endl;
			pdbWriter(prot, "current.pdb");
			prot->commitState();
			oldEnergy = energy;
		}
		else
		{
			double deltaE = energy - oldEnergy;
			double prob = pow(2.718, -beta*deltaE);
			if (prob > ranNumber.getNext())
			{
				cout << "ACCEPTED" << endl;
				pdbWriter(prot, "current.pdb");
				prot->commitState();
				oldEnergy = energy;
			}
			else
			{
				cout << "REJECT" << endl;
				prot->undoState();
			}
		}
		pdbWriter(prot, "current.pdb");
	}

	pdbWriter(prot,outputFileName);
	return 0;
}
