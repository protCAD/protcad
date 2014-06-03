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
	microEnvironment::setScaleFactor(2.0);

	//pBrian->setDistanceDependenceOn();
	pBrian->setDistanceDependenceOff();

	string path_level1 = "/home/summa/baker/1gb1-Rosetta";
	string listname = path_level1 + "/list";
	ifstream listfile;
	listfile.open(listname.c_str());
	string linebuffer;
	while (listfile >> linebuffer)
	{
				string filename = path_level1 + "/" + linebuffer;
				pProt = pdbReader(filename);

				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}

				double josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << atom::getHowMany() << "  " << josh << endl;
				delete pProt;
				pProt = 0;
	}
	listfile.close();
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
