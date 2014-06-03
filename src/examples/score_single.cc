#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include "microEnvDB.h"
#include "microEnvironment.h"
#include <time.h>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{

        if (argc != 4)
        {   cout << "Usage: score_multiple radius skip 1:calculate_all" << endl;
                exit (1);
        }

        double tempDouble;
        int tempInt;
        sscanf(argv[1],"%lf", &tempDouble);
        double radius = tempDouble;
        sscanf(argv[2],"%i", &tempInt);
        int skip = tempInt;
        sscanf(argv[3],"%i", &tempInt);
        int flag = tempInt;

        cout << "Radius = " << radius << endl;
        cout << "Skip = " << skip << endl;


	molecule* pProt = pdbReader("1UBI.pdb");
	delete pProt;
	microEnvironment* pBrian = new microEnvironment(radius,skip);

	string path_level1 = "/usr4/people/summa4/dd/single";
	string listname = path_level1 + "/list";
	ifstream listfile;
	listfile.open(listname.c_str());
	//cout << listname << endl;
	string linebuffer;
	while (listfile >> linebuffer)
	{
		string path_level2 = path_level1 + "/" + linebuffer;
		string listname2 = path_level2 + "/list";
		//cout << listname2 << endl;
		ifstream listfile2;
		listfile2.open(listname2.c_str());
		string linebuffer2;
		while (listfile2 >> linebuffer2)
		{	
			string path_level3 = path_level2 + "/" + linebuffer2;
			string listname3 = path_level3 + "/list";
			//cout << listname3 << endl;
			ifstream listfile3;
			listfile3.open(listname3.c_str());
			string linebuffer3;
			while (listfile3 >> linebuffer3)
			{
				string filename = path_level3 + "/" + linebuffer3;
				//cout << "About to read " << filename << endl;
				pProt = pdbReader(filename);

				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}

				pBrian->initialize();

				double microE = pBrian->calculateEnergy(pProt);

                                double rocco = 0.0;
                                double isabella=0.0;
                                double slovic=0.0;

                                if (flag==1)
                                {
                                amberVDW::setScaleFactor(1.0);
                                amberElec::setScaleFactor(0.0);
                                pmf::setScaleFactor(0.0);
                                rotamer::setScaleFactor(0.0);

                                rocco = static_cast<protein*>(pProt)->intraEnergy();

                                amberVDW::setScaleFactor(0.0);
                                amberElec::setScaleFactor(1.0);
                                pmf::setScaleFactor(0.0);
                                rotamer::setScaleFactor(0.0);

				isabella = static_cast<protein*>(pProt)->intraEnergy();

                                amberVDW::setScaleFactor(0.0);
                                amberElec::setScaleFactor(0.0);
                                pmf::setScaleFactor(1.0);
                                rotamer::setScaleFactor(0.0);

                                slovic = static_cast<protein*>(pProt)->intraEnergy();

                                }

                                cout << filename << "  " << atom::getHowMany() ;
				cout << "  " << microE <<  "  " << rocco ;
                                cout << "  " << isabella << "  " << rocco + isabella ;
				cout << "  " << slovic << endl;

				//printPerAtomEnergy(pBrian);
				delete pProt;
				pProt = 0;
			}
			listfile3.close();
		}
		listfile2.close(); 
	}
	listfile.close();
	return 0;
}

void printPerAtomEnergy(microEnvironment* pBrian)
{
		vector<double> perAtomEnergy = pBrian->getEnergyPerAtom();
		vector< vector<UInt> > perAtomEnv = pBrian->getEnvPerAtom();
		vector<UInt> atomTypes = pBrian->getEnvAssignedAtomTypes();
		for (UInt i=0; i< perAtomEnergy.size(); i++)
		{
			cout << i << "  - " << atomTypes[i] << " | ";
			for (UInt j=0; j<perAtomEnv[i].size(); j++)
			{	cout << perAtomEnv[i][j] << " ";
			}
			cout << perAtomEnergy[i] << endl;
		}
}
