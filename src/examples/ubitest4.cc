#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argc, char* argv[])
{	
	if (argc != 7)
	{	cout << "Usage: ubitest pmfScale microEnvScale VDWscale VDWRadiusScale rotamerScale randomSeed" << endl;
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
	double vdwRadiusScale = tempDouble;
	sscanf(argv[5],"%lf", &tempDouble);
	double rotamerScale = tempDouble;
	sscanf(argv[6],"%i", &tempInt);
	int randomSeed = tempInt;

	pmf::setScaleFactor(pmfScale);
	microEnvironment::setScaleFactor(microScale);
	amberVDW::setScaleFactor(vdwScale);
	amberVDW::setRadiusScaleFactor(vdwRadiusScale);
	rotamer::setScaleFactor(rotamerScale);
	//cout << "pmf scaling factor = " << pmf::getScaleFactor() << endl;
	//cout << "microEnvironent scaling factor = " << microEnvironment::getScaleFactor() << endl;
	//cout << "vanDerWaals scaling factor = " << amberVDW::getScaleFactor() << endl;
	//cout << "vanDerWaals radius scaling factor = " << amberVDW::getRadiusScaleFactor() << endl;
	//cout << "rotamer scaling factor = " << rotamer::getScaleFactor() << endl;
	//cout << "randomSeed = " << randomSeed << endl;

	int temp = 0;
	temp = static_cast<protein*>(pTheProtein1)->getIndexFromResNum(0,27);
	static_cast<protein*>(pTheProtein1)->activateForRepacking(0,temp);
	static_cast<protein*>(pTheProtein1)->setOnlyNativeIdentity(0,temp);

	double lowest = 500.0;
	int lowestrotamer = 0;
	string outfile = "testit.out";
	char fmodchar[2];

	microEnvDB::setUseSingleBodyEnergyOn();

	for (UInt i=0; i<81; i++)
	{
		//outfile = "UBI_out";
		static_cast<protein*>(pTheProtein1)->setRotamerWBC(0,temp,0,i);
		double theEnergy = static_cast<protein*>(pTheProtein1)->intraEnergy();
		//sprintf(fmodchar,"%i",i);
		//string fmodstring = fmodchar;
                //outfile += fmodstring;
                //outfile += ".pdb";

		//cout << outfile << "  ";
		//pdbWriter(static_cast<protein*>(pTheProtein1),outfile);
		if (theEnergy < lowest)
		{
			lowest = theEnergy;
			lowestrotamer = i;
		}
	//	cout << i << "  " << theEnergy << endl;
	}
	if (lowestrotamer == 67)
	{	
		ofstream fout;
		fout.open(outfile.c_str(),ios::app);
		fout << "SUCCESS " <<  pmf::getScaleFactor() << " " ;
		fout <<  microEnvironment::getScaleFactor() << " ";
		fout << amberVDW::getScaleFactor() << " ";
		fout << amberVDW::getRadiusScaleFactor() << " ";
		fout <<  rotamer::getScaleFactor() << endl;
		fout.close();
	}


	//pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	delete pTheProtein1;
/*
	time(&endTime);
	outfilename = "time.dat";
	ofstream oFile;
	oFile.open(outfilename.c_str());
	oFile << "startTime = " << startTime << "\n";
	oFile << "endTime = " << endTime << "\n";
	runTime = endTime - startTime;
	oFile << "runTime = " << runTime << "\n";
	oFile.close();
*/
	return 0;
}
