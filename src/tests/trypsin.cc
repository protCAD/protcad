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

	UIntVec active;

	active.push_back(prot->getIndexFromResNum(0,116));
	active.push_back(prot->getIndexFromResNum(0,134));
	active.push_back(prot->getIndexFromResNum(0,136));
	active.push_back(prot->getIndexFromResNum(0,170));
	active.push_back(prot->getIndexFromResNum(0,174));

	UIntVec allowed;
	allowed.push_back(A);
	allowed.push_back(I);
	allowed.push_back(L);
	allowed.push_back(F);
	allowed.push_back(V);
	allowed.push_back(W);

	ran ranNumber;
	ranNumber.setSeed(132);

	for (UInt i = 0; i < active.size(); i ++)
	{	
		prot->activateForRepacking(0,active[i]);
		prot->mutate(0,active[i],A);
	}

	double oldEnergy = 10E5;
	for (UInt i = 0; i < 1000; i ++)
	{
		UInt pos = (UInt)(ranNumber.getNext()*active.size());
		if (pos >= active.size()) pos = 0;
		UInt choice = (UInt)(ranNumber.getNext()*allowed.size());
		if (pos >= allowed.size()) choice = 0;
		prot->saveCurrentState();
		
		double beta = 1E-4 - i*(1E-4 - 1)/1000;
		cout << "mutating " << active[pos] << " to " << allowed[choice] << endl;
		prot->mutate(0, active[pos], allowed[choice]);
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
	}
	pdbWriter(prot,outputFileName);
	return 0;
}
