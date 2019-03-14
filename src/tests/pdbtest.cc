#include "pdbReader.h"
#include "pdbWriter.h"
#include "ran.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	string filename = "temp.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	static_cast<protein*>(pTheProtein1)->activateAllChainPosition(0);
	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();
	static_cast<protein*>(pTheProtein1)->listSecondaryStructure();

	ran rnum(53490);

	double oldEnergy;
	oldEnergy = pTheProtein1->intraEnergy();
	double newEnergy;

	time(&startTime);
	for (unsigned int i=0; i<10; i++)
	{	int success = pTheProtein1->modify(rnum);
		if (success != -1)
		{	cout << "Before change : " << oldEnergy;
			newEnergy = pTheProtein1->intraEnergy();
			cout << " After change : " << newEnergy;
			if (newEnergy < oldEnergy)
			{	pTheProtein1->acceptModification();
				oldEnergy = newEnergy;
				cout << " Change accepted : " << oldEnergy << endl;
			}
			else
			{	pTheProtein1->rejectModification();
				//oldEnergy = pTheProtein1->intraEnergy();
				cout << " Change rejected " << endl;
				cout << pTheProtein1->intraEnergy() << endl;
			}
		}
	}

	string outfilename = "jason.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	delete pTheProtein1;

	time(&endTime);
	outfilename = "time.dat";
	ofstream oFile;
	oFile.open(outfilename.c_str());
	oFile << "startTime = " << startTime << "\n";
	oFile << "endTime = " << endTime << "\n";
	runTime = endTime - startTime;
	oFile << "runTime = " << runTime << "\n";
	oFile.close();

	return 0;
}
