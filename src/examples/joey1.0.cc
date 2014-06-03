#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "mmo_core.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,251-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,251-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,138-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,138-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,142-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,142-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,115-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,115-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,125-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,125-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,202-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,202-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,212-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,212-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,216-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,216-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,223-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,223-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,211-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,211-95);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,177-95);
	static_cast<protein*>(pTheProtein1)->setOnlyHydrophilic(0,177-95);
	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();

	ensemble* pEnsemble = new ensemble;
    pEnsemble->add(pTheProtein1);
    annealer* pAnneal = new annealer(pEnsemble);
    string outfilename;
    pAnneal->run(600.0, 10.0, 20000, 1197);
    pAnneal->run(10.0, 10.0, 2000, 1347);
    outfilename = "joeyout1.pdb";
    pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);

    pAnneal->run(600.0, 10.0, 20000, 12338);
    pAnneal->run(10.0,10.0,2000,9348);
    outfilename = "joeyout2.pdb";
    pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);

    pAnneal->run(600.0, 10.0, 20000, 21994);
    pAnneal->run(10.0, 10.0, 2000, 44278);
    outfilename = "joeyout3.pdb";
    pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);

    pAnneal->run(600.0, 10.0, 20000, 63420);
    pAnneal->run(10.0, 10.0, 2000, 3989);
    outfilename = "joeyout4.pdb";
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
