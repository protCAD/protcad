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

	string path_level1 = "/usr4/people/summa4/dd/multiple";
	string listname = path_level1 + "/list";
	ifstream listfile;
	listfile.open(listname.c_str());
	//cout << listname << endl;
	string linebuffer;
	while (listfile >> linebuffer)
	{
		string path_level2 = path_level1 + "/" + linebuffer;
		string listname2 = path_level2 + "/list";
		//cout << listname2 << endl;
		ifstream listfile2;
		listfile2.open(listname2.c_str());
		string linebuffer2;
		while (listfile2 >> linebuffer2)
		{	
			string path_level3 = path_level2 + "/" + linebuffer2;
			string listname3 = path_level3 + "/list";
			//cout << listname3 << endl;
			ifstream listfile3;
			listfile3.open(listname3.c_str());
			string linebuffer3;
			while (listfile3 >> linebuffer3)
			{
				string filename = path_level3 + "/" + linebuffer3;
				pProt = pdbReader(filename);
				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}
				double josh = pBrian->calculateEnergy(pProt);
				cout << filename << "  " << josh << endl;
				//printPerAtomEnergy(pBrian);
				delete pProt;
				pProt = 0;
			}
			listfile3.close();
		}
		listfile2.close(); 
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
