#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "joslyn_start.pdb";

	molecule* pTheProtein1 = pdbReader(filename);

	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	int tempint;

//-------you can edit from here---------------------

//  residue 46 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,46);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 46 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,46);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 43 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,43);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 43 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,43);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 12 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,12);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 12 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,12);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 22 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,22);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 22 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,22);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 9 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,9);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 9 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,9);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 25 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,25);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 25 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,25);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 29 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,29);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 29 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,29);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 5 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,5);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 5 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,5);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 39 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,39);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 39 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,39);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 32 chain A: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,32);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,tempint);

//  residue 32 chain B: hydrophobe
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,32);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(1,tempint);

//  residue 36 chain A: anything
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,36);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,tempint);

//  residue 36 chain B: anything
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,36);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(1,tempint);

//-------to here !!!!-------------

	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();

	cout << "Allowed Rotamers" << endl << endl;

	for (UInt i=0;i< 40; i++)
	{ cout << " residue " << i << endl;
	  static_cast<protein*>(pTheProtein1)->printAllowedRotamers(0,i);
	 }

	cout << endl << endl;

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;

	pAnneal->run(600.0, 10.0, 10000, 43725);
	pAnneal->run(10.0, 10.0, 4000, 57264);
	outfilename = "joslynout1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 10000, 17982);
	pAnneal->run(10.0,10.0,4000, 5427);
	outfilename = "joslynout2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 10000, 7269);
	pAnneal->run(10.0, 10.0, 4000, 6300);
	outfilename = "joslynout3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 10000, 9376);
	pAnneal->run(10.0, 10.0, 4000, 2444);
	outfilename = "joslynout4.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);

	pAnneal->run(600.0, 10.0, 10000, 3376);
	pAnneal->run(10.0, 10.0, 4000, 24834);
	outfilename = "joslynout5.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);

	pAnneal->run(600.0, 10.0, 10000, 39376);
	pAnneal->run(10.0, 10.0, 4000, 23334);
	outfilename = "joslynout6.pdb";
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
