#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"

int main (int argc, char* argv[])
{
	if (argc < 7)
	{
		cout << "Usage:  desjarlais (1)pmfScale\n\t(2)microEnvScale \n\t(3)vdWScale \n\t(4)amberElecScale \n\t(5)distance-dielectric-on (0 off, 1 on) \n\t(6) # of cycles \n\t(7)randomseed" << endl;
	exit(1);
	}

	string fileName = "1UBI_h.pdb"; // ubiquitin pdb file
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
	sscanf(argv[7], "%i", &tmpInt);
	int randomSeed = tmpInt; cout << "random seed: " << randomSeed << endl;
	sscanf(argv[6], "%i", &tmpInt);
	int cycles = tmpInt; cout << "number of cycles: " << cycles << endl;

	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,3));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,5));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,13));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,15));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,17));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,23));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,26));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,30));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,43));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,50));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,56));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,61));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,67));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,69));

    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,3));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,5));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,13));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,15));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,17));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,23));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,26));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,30));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,43));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,50));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,56));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,61));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,67));
    pTheProtein->setOnlyROCHydrophobic(0, pTheProtein->getIndexFromResNum(0,69));

    ensemble *pEnsemble = new ensemble;
    pEnsemble->add(pTheMolecule);
    annealer* pAnneal = new annealer(pEnsemble);

    for (UInt i = 0; i < cycles; i ++)
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
