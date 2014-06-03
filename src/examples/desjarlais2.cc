#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include "deadEndEliminator.h"

int main (int argc, char* argv[])
{
	if (argc < 9)
	{
		cout << "Usage:  desjarlais2 (1-8) \n\t(1)pmfScale \n\t(2)microEnvScale \n\t(3)vdWScale \n\t(4)amberElecScale \n\t(5)distance-dielectric-on (0 off, 1 on) \n\t(6) radius scaling factor \n\t(7) vdw linear repulsion dampening on (0 off, 1 on) \n\t(8) # of cycles \n\t(9)randomseed" << endl;
	exit(1);
	}
	string fileName = "1UBI.pdb"; // ubiquitin pdb file
	molecule* pTheMolecule = pdbReader(fileName);
	cout << "protein read in" << endl;
	protein* pTheProtein = static_cast<protein*>(pTheMolecule);
	
	if (pTheProtein == 0) return 1;

    cout << "number of residues generated:  " << residue::getHowMany() << endl;
    cout << "number of chains generated:  " << chain::getHowMany() << endl;
    cout << "number of residue templates generated: " << residueTemplate::getHowMany() << endl;

	double tmpDbl; int tmpInt;
	sscanf(argv[1], "%lf", &tmpDbl);
	pmf::setScaleFactor(tmpDbl);
	rotamer::setScaleFactor(1.0);
	sscanf(argv[2], "%lf", &tmpDbl);
	microEnvironment::setScaleFactor(tmpDbl);
	sscanf(argv[3], "%lf", &tmpDbl);
	amberVDW::setScaleFactor(tmpDbl);
	sscanf(argv[4], "%lf", &tmpDbl);
	amberElec::setScaleFactor(0.0);
	sscanf(argv[5], "%lf", &tmpDbl);
	if (tmpDbl == 0) amberElec::distanceDependanceOff();
	if (tmpDbl == 1) amberElec::distanceDependanceOn();
	sscanf(argv[6], "%lf", &tmpDbl);
	amberVDW::setRadiusScaleFactor(tmpDbl);
	sscanf(argv[7], "%lf", &tmpDbl);
	if (tmpDbl == 1) amberVDW::setLinearRepulsionDampeningOn();
	if (tmpDbl == 0) amberVDW::setLinearRepulsionDampeningOff();
	sscanf(argv[9], "%i", &tmpInt);
	int randomSeed = tmpInt; cout << "random seed: " << randomSeed << endl;
	sscanf(argv[8], "%i", &tmpInt);
	int cycles = tmpInt; cout << "number of cycles: " << cycles << endl;


	pTheProtein->activateAllForRepacking(0);


    ensemble *pEnsemble = new ensemble;
    pEnsemble->add(pTheMolecule);
    annealer* pAnneal = new annealer(pEnsemble);

    for (int i = 0; i < cycles; i ++)
    {
        pAnneal->run(500.0, 100.0, 2000, randomSeed+(i));
        pAnneal->run(100.0, 10.0, 2000, randomSeed+(i*i));
        pAnneal->run(10.0, 10.0, 1000, randomSeed+(10*i));
        string outFileName = "UBIout_ " ;
        outFileName[7] = 'a' + i;
        outFileName[8] = '\0';
        outFileName += ".pdb";

        pdbWriter(pTheProtein, outFileName);
  	}

	delete pAnneal;
	delete pTheProtein;
	delete pTheMolecule;	
	delete pEnsemble;

	return 0;
}
