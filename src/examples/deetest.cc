#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "microEnvDB.h"
#include "microEnvironment.h"
#include "deadEndEliminator.h"
#include <time.h>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{	
	if (argc != 3)
	{	cout << "Usage: deetest pmfScale microEnvScale" << endl;
		exit (1);
	}	
		
	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "1UBI.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;
	cout << "number of res templates generated:" << residueTemplate::getHowMany() << endl;

	double tempDouble;
	sscanf(argv[1],"%lf", &tempDouble);
	double pmfScale = tempDouble;
	sscanf(argv[2],"%lf", &tempDouble);
	double microScale = tempDouble;
	pmf::setScaleFactor(pmfScale);
	microEnvironment::setScaleFactor(microScale);
	cout << "pmf scaling factor = " << pmf::getScaleFactor() << endl;
	cout << "microEnvironent scaling factor = " << microEnvironment::getScaleFactor() << endl;

	string outfilename;

	UInt tempint;
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,27);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
/*
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
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,52);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,53);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,54);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,55);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,56);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,57);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,58);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,59);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,61);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,62);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,64);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,65);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,66);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,67);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,68);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,69);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,70);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,71);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
	tempint = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,72);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,tempint);
*/
	
	deadEndEliminator jason(static_cast<protein*>(pTheProtein1));
	jason.run(0,100.0);

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
