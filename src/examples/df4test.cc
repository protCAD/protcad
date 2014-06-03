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
	string filename = "repackme.pdb";

	molecule* pTheProtein1 = pdbReader(filename);

	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

//	residue 3: anything
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,2);

//  residue 67: anything
    static_cast<protein*>(pTheProtein1)->activateChainPosition(0,66);

//  residue 134: anything
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,100);

//  residue 136: anything
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,102);

//  residue 6: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,5);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,5);

//  residue 64: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,63);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,63);

//  residue 139 : hydrophobic
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,105);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,105);

	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();
	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;

	outfilename = "DF4start.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	pAnneal->run(600.0, 10.0, 20000, 104368);
	pAnneal->run(10.0, 10.0, 2000, 73462);
	outfilename = "DF4out1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 95122);
	pAnneal->run(10.0,10.0,2000, 3300);
	outfilename = "DF4out2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 52569);
	pAnneal->run(10.0, 10.0, 2000, 60270);
	outfilename = "DF4out3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 4465);
	pAnneal->run(10.0, 10.0, 2000, 202444);
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
