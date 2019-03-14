#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include "microEnvDB.h"
#include "microEnvironment.h"
#include <time.h>
#include <fstream>
#include <string>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);

	molecule* pProt = pdbReader("1UBI.pdb");
	delete pProt;
	microEnvironment* pBrian = new microEnvironment(4.0,2);

				string filename = "/home/summa/misfolds/alpha2D/r1_noh.pdb";
				pProt = pdbReader(filename);
				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}
				double josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << josh << endl;
				delete pProt;
				pProt = 0;

				filename = "/home/summa/misfolds/alpha2D/r2_noh.pdb";
				pProt = pdbReader(filename);
				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}
				josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << josh << endl;
				delete pProt;
				pProt = 0;

				filename = "/home/summa/misfolds/alpha2D/s1_noh.pdb";
				pProt = pdbReader(filename);
				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}
				josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << josh << endl;
				delete pProt;
				pProt = 0;

				filename = "/home/summa/misfolds/alpha2D/s2_noh.pdb";
				pProt = pdbReader(filename);
				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}
				josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << josh << endl;
				delete pProt;
				pProt = 0;

	time(&endTime);
	cout << "startTime = " << startTime << endl;
	cout << "endTime = " << endTime << endl;
	runTime = endTime - startTime;
	cout << "runTime = " << runTime << endl;

	return 0;
}

void printPerAtomEnergy(microEnvironment* pBrian)
{
		vector<double> perAtomEnergy = pBrian->getEnergyPerAtom();
		vector< vector<UInt> > perAtomEnv = pBrian->getEnvPerAtom();
		vector<UInt> atomTypes = pBrian->getEnvAssignedAtomTypes();
		for (UInt i=0; i< perAtomEnergy.size(); i++)
		{
			cout << i << "  - " << atomTypes[i] << " | ";
			for (UInt j=0; j<perAtomEnv[i].size(); j++)
			{	cout << perAtomEnv[i][j] << " ";
			}
			cout << perAtomEnergy[i] << endl;
		}
}
