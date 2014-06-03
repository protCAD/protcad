#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argc, char* argv[])
{	
	if (argc != 6)
	{	cout << "Usage: ubitest pmfScale microEnvScale VDWscale rotamerScale randomSeed" << endl;
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
	int tempInt;
	sscanf(argv[1],"%lf", &tempDouble);
	double pmfScale = tempDouble;
	sscanf(argv[2],"%lf", &tempDouble);
	double microScale = tempDouble;
	sscanf(argv[3],"%lf", &tempDouble);
	double vdwScale = tempDouble;
	sscanf(argv[4],"%lf", &tempDouble);
	double rotamerScale = tempDouble;
	sscanf(argv[5],"%i", &tempInt);
	int randomSeed = tempInt;

	pmf::setScaleFactor(pmfScale);
	microEnvironment::setScaleFactor(microScale);
	amberVDW::setScaleFactor(vdwScale);
	rotamer::setScaleFactor(rotamerScale);
	cout << "pmf scaling factor = " << pmf::getScaleFactor() << endl;
	cout << "microEnvironent scaling factor = " << microEnvironment::getScaleFactor() << endl;
	cout << "vanDerWaals scaling factor = " << amberVDW::getScaleFactor() << endl;
	cout << "rotamer scaling factor = " << rotamer::getScaleFactor() << endl;
	cout << "randomSeed = " << randomSeed << endl;

	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,2);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,2);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,4);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,4);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,12);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,12);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,14);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,14);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,16);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,16);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,22);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,22);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,25);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,25);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,29);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,29);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,42);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,42);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,49);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,49);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,55);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,55);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,60);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,60);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,66);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,66);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,68);
	static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,68);

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;
	pAnneal->run(400.0, 10.0, 7000, randomSeed);
	pAnneal->run(10.0, 10.0, 2000);
	outfilename = "ubiout1.pdb";
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
