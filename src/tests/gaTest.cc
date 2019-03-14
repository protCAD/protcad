#include <pdbReader.h>
#include <pdbWriter.h>
#include <atomIterator.h>
#include <ran.h>
#include <time.h>
#include <myGenome.h>


int main()
{
	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "1UBI.pdb";
  
//	exit(1);
	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;
  
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,2);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,2);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,4);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,4);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,12);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,12);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,14);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,14);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,16);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,16);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,22);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,22);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,25);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,25);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,29);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,29);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,42);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,42);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,49);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,49);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,55);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,55);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,60);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,60);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,66);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,66);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,68);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,68);
	static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();
  
	myGenome theGA = myGenome(pTheProtein1);
	theGA.setNumGen(10);
	theGA.setPopSize(10);
	theGA.run();

	string outfilename = "JoshsBaby.pdb";
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
