#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "rub_dimer_noh.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;
	cout << "number of res templates generated:" << residueTemplate::getHowMany() << endl;

	UInt tempint;


    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,35);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,36);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,37);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,39);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,40);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,42);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,43);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,44);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,45);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,46);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,47);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,48);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,49);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,50);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,51);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);

    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,35);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,36);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,37);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,39);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,40);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,42);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,43);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,44);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,45);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,46);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,47);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,48);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,49);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,50);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);
    tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(1,51);
    static_cast<protein*>(pTheProtein1)->activateForRepacking(1,tempint);

	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;
	pAnneal->run(600.0, 10.0, 20000, 111197);
	pAnneal->run(10.0, 10.0, 2000, 13417);
	outfilename = "zk1out1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
/*	
	pAnneal->run(600.0, 10.0, 20000, 12398);
	pAnneal->run(10.0,10.0,2000,9348);
	outfilename = "zk1out2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 1994);
	pAnneal->run(10.0, 10.0, 2000, 4478);
	outfilename = "zk1out3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(600.0, 10.0, 20000, 634200);
	pAnneal->run(10.0, 10.0, 2000, 39894);
	outfilename = "zk1out4.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
*/
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
