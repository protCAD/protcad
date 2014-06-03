#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
int main (int argc, char* argv[])
{
    if (argc < 10)
    {
        cout << "Usage:  desjarlais2 (1-10) \n\t(1)filename \n\t(2)pmfScale \n\t(3)microEnvScale \n\t(4)vdWScale \n\t(5)amberElecScale \n\t(6)distance-dielectric-on (0 off, 1 on) \n\t(7) radius scaling factor \n\t(8) vdw linear repulsion dampening on (0 off, 1 on) \n\t(9) # of cycles \n\t(10)randomseed" << endl;
    exit(1);
    }
    string fileName = argv[1];
	    molecule* pTheMolecule = pdbReader(fileName);
    cout << fileName << " protein read in" << endl;

    protein* pTheProtein = static_cast<protein*>(pTheMolecule);
	
    if (pTheProtein == 0) return 1;
    cout << "number of residues generated:  " << residue::getHowMany() << endl;
    cout << "number of chains generated:  " << chain::getHowMany() << endl;
    cout << "number of residue templates generated: " << residueTemplate::getHowMany() << endl;

    double tmpDbl; int tmpInt;
    sscanf(argv[2], "%lf", &tmpDbl);
    pmf::setScaleFactor(tmpDbl);
    rotamer::setScaleFactor(1.0);
    sscanf(argv[3], "%lf", &tmpDbl);
    microEnvironment::setScaleFactor(tmpDbl);
    sscanf(argv[4], "%lf", &tmpDbl);
    amberVDW::setScaleFactor(tmpDbl);
    sscanf(argv[5], "%lf", &tmpDbl);
    amberElec::setScaleFactor(tmpDbl);
    sscanf(argv[6], "%lf", &tmpDbl);
    if (tmpDbl == 0) amberElec::distanceDependanceOff();
    if (tmpDbl == 1) amberElec::distanceDependanceOn();
    sscanf(argv[7], "%lf", &tmpDbl);
    amberVDW::setRadiusScaleFactor(tmpDbl);
    sscanf(argv[8], "%lf", &tmpDbl);
    if (tmpDbl == 1) amberVDW::setLinearRepulsionDampeningOn();
    if (tmpDbl == 0) amberVDW::setLinearRepulsionDampeningOff();
    sscanf(argv[10], "%i", &tmpInt);
    int randomSeed = tmpInt; cout << "random seed: " << randomSeed << endl;
    sscanf(argv[9], "%i", &tmpInt);
    int cycles = tmpInt; cout << "number of cycles: " << cycles << endl;

    ensemble *pEnsemble = new ensemble;
    pEnsemble->add(pTheMolecule);
    annealer* pAnneal = new annealer(pEnsemble);

	amberVDW::setAttractionScaleFactor(0.0);

	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,26), pTheProtein->getIndexFromResNum(0,41)); 

	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,26), 4); // disallow cysteine
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,27), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,28), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,29), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,30), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,31), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,32), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,33), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,34), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,35), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,36), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,37), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,38), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,39), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,40), 4);
	pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,41), 4);

        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,26), 7); // disallow glycine
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,27), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,28), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,29), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,30), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,31), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,32), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,33), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,34), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,35), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,36), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,37), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,38), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,39), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,40), 7);
        pTheProtein->setResNotAllowed(0, pTheProtein->getIndexFromResNum(0,41), 7);

    for (int i = 0; i < cycles; i ++)
    {
        pAnneal->run(500.0, 100.0, 2000, randomSeed+(i));
        pAnneal->run(100.0, 10.0, 2000, randomSeed+(i*i));
        pAnneal->run(10.0, 10.0, 1000, randomSeed+(10*i));
        string outFileName = "repack_ " ;
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
