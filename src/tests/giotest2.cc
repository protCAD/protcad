#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "dsd1.1.pdb";

	molecule* pTheProtein1 = pdbReader(filename);

	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	int tempint;
//  residue 46 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,46);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 43 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,43);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 36 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,36);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 12 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,12);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 22 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,22);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 9 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,9);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 25 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,25);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 29 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,29);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 5 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,5);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 39 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,39);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 32 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,32);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 2 chain B: anything
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,2);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);

//  residue 3 chain B: anything
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,3);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);

//  residue 36 chain B: anything
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,36);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);

//-------start-------------
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
