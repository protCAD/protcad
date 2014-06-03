#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	//string filename = "pdb1dhp.ent";
	//string filename = "dsdhalren.pdb";
	string filename = "forchris.pdb";

	molecule* pTheProtein1 = pdbReader(filename);

	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

//	residue 2: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,1);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,1);

//  residue 5: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,4);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,4);

//  residue 9: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,8);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,8);

//  residue 12: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,11);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,11);

//  residue 22: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,21);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,21);

//  residue 25: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,24);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,24);

//  residue 29: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,28);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,28);

//  residue 32: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,31);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,31);

//  residue 38 Chain B: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,2);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,2);

//  residue 41 Chain B: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,5);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,5);

//  residue 45 Chain B: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,9);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,9);

//  residue 48 Chain B: hydrophobe
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,12);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,12);

	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();
	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;

	outfilename = "GIOstart.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	pAnneal->run(600.0, 10.0, 20000, 104368);
	pAnneal->run(10.0, 10.0, 2000, 73462);
	outfilename = "GIOout1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 95122);
	pAnneal->run(10.0,10.0,2000, 3300);
	outfilename = "GIOout2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 52569);
	pAnneal->run(10.0, 10.0, 2000, 60270);
	outfilename = "GIOout3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 4465);
	pAnneal->run(10.0, 10.0, 2000, 202444);
	outfilename = "GIOout4.pdb";
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
