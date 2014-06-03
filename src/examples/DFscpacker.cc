#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "annealer.h"
#include <sstream>


int main(int argc, char* argv[])
{
   residue::setCutoffDistance(6.0);
    pmf::setScaleFactor(0);
    rotamer::setScaleFactor(1.0);
    microEnvironment::setScaleFactor(0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
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

	active.push_back(prot->getIndexFromResNum(0,20));
	active.push_back(prot->getIndexFromResNum(0,23));
	active.push_back(prot->getIndexFromResNum(0,24));
	active.push_back(prot->getIndexFromResNum(0,31));
	active.push_back(prot->getIndexFromResNum(0,32));
	active.push_back(prot->getIndexFromResNum(0,29));
	active.push_back(prot->getIndexFromResNum(0,33));

	UIntVec allowed;
	allowed.push_back(A);
	allowed.push_back(I);
	allowed.push_back(L);
	allowed.push_back(F);
	allowed.push_back(V);
	allowed.push_back(M);

	prot->silenceMessages();
	
	ran ranNumber;
	ranNumber.setSeed(132);

	for (UInt i = (UInt)prot->getIndexFromResNum(0,20); i < (UInt)prot->getIndexFromResNum(0,38); i ++)
	{
		prot->activateForRepacking(0,i);
		if (i != (UInt)prot->getIndexFromResNum(0,27) && i != (UInt)prot->getIndexFromResNum(0,28))
		{
			prot->setCanonicalHelixRotamersOnly(0,i);
		}
	}
	
	for (UInt i = 0; i < active.size(); i ++)
	{
		prot->mutate(0,active[i], A);
	}
	prot->optimizeRotamers();
	double alaE = prot->intraEnergy();
	
	vector < vector < bool > > selfGood(0);
	
	for (UInt i = 0; i < active.size(); i ++)
	{
		vector < bool > temp(0);
		for (UInt j = 0; j < allowed.size(); j ++)
		{
			prot->mutate(0, active[i], allowed[j]);
			prot->optimizeRotamers();
			double energy = prot->intraEnergy() - alaE;
			if (energy > 50.0)
			{
				temp.push_back(false);
			}
			else temp.push_back(true);
			prot->mutate(0, active[i], A);
		}
		selfGood.push_back(temp);
	}	
	UInt allowedsize = allowed.size();
	UInt activesize = active.size();
	double lowestE = 5000.0;

	for (UInt i = 0; i < active.size(); i ++)
	{
		cout << "truth table " << i << " ";
		for (UInt j = 0; j < allowed.size(); j ++)
		{
			if (selfGood[i][j]) cout << "T ";
			else cout << "F ";
		}
		cout << endl;
	}
	
	ofstream fout ("energylog.out", ios::app);
	for (UInt i = 0; i < (UInt)pow((double)allowedsize,(double)activesize); i ++)
	{
		UInt temp = i;
		bool calculate = true;
		for (UInt j = 0; j < activesize; j ++)
		{
			if (!selfGood[j][temp%allowedsize]) calculate = false;
			temp = (UInt)((float)temp/(float)allowedsize);
		}		
		if (calculate)			
		{
			temp = i;	
			for (UInt j = 0; j < activesize; j ++)
			{
				UInt resID = allowed[temp%allowedsize];
				prot->mutate(0, active[j], resID);
				temp = (UInt)((float)temp/(float)allowedsize);
			}
		
			prot->optimizeRotamers();
			double energy = prot->intraEnergy();
			if (energy < lowestE)
			{
				pdbWriter(prot, "lowest.pdb");
				cout << "NEW LOW ";
				lowestE = energy;
			}
		
		
			cout << i << ":  ";
			fout << i << ":  ";
			for (UInt j = 0; j < activesize; j ++)
			{
				fout << prot->getTypeStringFromResNum(0, active[j]) << " ";
				cout << prot->getTypeStringFromResNum(0, active[j]) << " ";
			}
			fout << energy << endl;
			cout << energy << endl;
		}
	}
	pdbWriter(prot,outputFileName);
	return 0;
}
