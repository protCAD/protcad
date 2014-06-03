#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argc, char* argv[])
{	
	if (argc != 5)
	{	cout << "Usage: ubitest pmfScale microEnvScale VDWscale randomSeed" << endl;
		exit (1);
	}	
		
	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "bhelix.pdb";

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
	sscanf(argv[4],"%i", &tempInt);
	int randomSeed = tempInt;

    pmf::setScaleFactor(pmfScale);
    microEnvironment::setScaleFactor(microScale);
	amberVDW::setScaleFactor(vdwScale);
    cout << "pmf scaling factor = " << pmf::getScaleFactor() << endl;
    cout << "microEnvironent scaling factor = " << microEnvironment::getScaleFactor() << endl;
    cout << "vanDerWaals scaling factor = " << amberVDW::getScaleFactor() << endl;
	cout << "randomSeed = " << randomSeed << endl;

	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,0);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,1);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,2);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,3);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,4);
	static_cast<protein*>(pTheProtein1)->activateChainPosition(0,5);
	static_cast<protein*>(pTheProtein1)->setAllBetaAminoAcids();

	static_cast<protein*>(pTheProtein1)->listAllowedRotamers(0,0);
	static_cast<protein*>(pTheProtein1)->listAllowedRotamers(0,1);
	static_cast<protein*>(pTheProtein1)->listAllowedRotamers(0,2);
	static_cast<protein*>(pTheProtein1)->listAllowedRotamers(0,3);
	static_cast<protein*>(pTheProtein1)->listAllowedRotamers(0,4);

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	annealer* pAnneal = new annealer(pEnsemble);
	string outfilename;
	pAnneal->run(400.0, 10.0, 7000, randomSeed);
	pAnneal->run(10.0, 10.0, 2000);
	outfilename = "betaout1.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(300.0, 10.0, 7000);
	pAnneal->run(10.0,10.0,2000);
	outfilename = "betaout2.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(200.0, 10.0, 7000);
	pAnneal->run(10.0, 10.0, 2000);
	outfilename = "betaout3.pdb";
	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	pAnneal->run(100.0, 10.0, 7000);
	pAnneal->run(10.0, 10.0, 2000);
	outfilename = "betaout4.pdb";
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
