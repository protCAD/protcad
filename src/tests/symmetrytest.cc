#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argc, char* argv[])
{	
	if (argc != 7)
	{	cout << "Usage: symmetrytest pmfScale microEnvScale VDWscale rotamerScale randomSeed elecScale" << endl;
		exit (1);
	}	
		
	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "newsymmetric.pdb";
	string outfilename;

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
	sscanf(argv[6],"%lf", &tempDouble);
	double elecScale = tempDouble;

	pmf::setScaleFactor(pmfScale);
	microEnvironment::setScaleFactor(microScale);
	amberVDW::setScaleFactor(vdwScale);
	rotamer::setScaleFactor(rotamerScale);
	amberElec::setScaleFactor(elecScale);
	cout << "pmf scaling factor = " << pmf::getScaleFactor() << endl;
	cout << "microEnvironent scaling factor = " << microEnvironment::getScaleFactor() << endl;
	cout << "vanDerWaals scaling factor = " << amberVDW::getScaleFactor() << endl;
	cout << "rotamer scaling factor = " << rotamer::getScaleFactor() << endl;
	cout << "randomSeed = " << randomSeed << endl;


	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(1,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(2,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(3,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(4,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(5,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(6,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(7,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(8,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(9,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(10,0);
	static_cast<protein*>(pTheProtein1)->symmetryLinkChainAtoB(11,0);

	static_cast<protein*>(pTheProtein1)->activateAllForRepacking(0);

	//static_cast<protein*>(pTheProtein1)->mutateWBC(0,5,3);

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	pAnneal->run(400.0, 10.0, 7000, randomSeed);
	pAnneal->run(10.0, 10.0, 2000);
	/*
	outfilename = "ubiout1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(300.0, 10.0, 7000);
	pAnneal->run(10.0,10.0,2000);
	outfilename = "ubiout2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(200.0, 10.0, 7000);
	pAnneal->run(10.0, 10.0, 2000);
	outfilename = "ubiout3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(100.0, 10.0, 7000);
	pAnneal->run(10.0, 10.0, 2000);
	*/
	outfilename = "symmetryout.pdb";
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
