#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "DF4_no_H.pdb";

	molecule* pTheProtein1 = pdbReader(filename);

	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,2);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,2);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,9);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,9);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,12);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,12);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,56);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,56);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,59);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,59);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,63);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,63);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,66);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,66);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,94);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,94);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,101);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,101);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,106);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,106);

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,110);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,110);

	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;

	outfilename = "DF4start.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	pAnneal->run(600.0, 10.0, 20000, 111197);
	pAnneal->run(10.0, 10.0, 2000, 13417);
	outfilename = "DF4out1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 12398);
	pAnneal->run(10.0,10.0,2000,9348);
	outfilename = "DF4out2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 1994);
	pAnneal->run(10.0, 10.0, 2000, 4478);
	outfilename = "DF4out3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 634200);
	pAnneal->run(10.0, 10.0, 2000, 39894);
	outfilename = "DF4out4.pdb";
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
