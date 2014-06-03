#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "vanDerWaals.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "temp.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	//static_cast<protein*>(pTheProtein1)->activateChainPosition(0,26);
	//static_cast<protein*>(pTheProtein1)->activateChainPosition(0,40);
	static_cast<protein*>(pTheProtein1)->activateAllChainPosition(0);

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	
	cout << " energy : " << pEnsemble->energy() << endl;
	annealer* pAnneal = new annealer(pEnsemble);
	pAnneal->run(5.0, 5.0, 1000, 32442);

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
