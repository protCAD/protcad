#include "pdbReader.h"
#include "atomIterator.h"
#include "annealer.h"
#include "ran.h"

int main(void)
{
	string fileName = "coil.pdb";
	molecule *pTheMolecule  = pdbReader(fileName);
	protein* pTheProtein = static_cast<protein*>(pTheMolecule);

	if (pTheProtein == 0) return 1;
	cout << "------- exercising the energy functions -------\n";
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(0.0);
	cout << "default distance dependance " << amberElec::isDistanceDependanceOn() << endl;
	amberElec::distanceDependanceOn();
	cout << "default amberElec::getScaleFactor() -> " << amberElec::getScaleFactor() << " energy -> " << pTheProtein->intraEnergy() << endl;
	amberElec::setScaleFactor(0.5);
	cout << "amberElec::setScaleFactor(0.5) -> " << amberElec::getScaleFactor() << " energy -> " << pTheProtein->intraEnergy() << endl;
	amberElec::distanceDependanceOff();
	amberElec::setScaleFactor(1.0);
    cout << "amberElec::getScaleFactor() -> " << amberElec::getScaleFactor() << " energy -> " << pTheProtein->intraEnergy() << endl;
    amberElec::setScaleFactor(0.5);
    cout << "amberElec::setScaleFactor(0.5) -> " << amberElec::getScaleFactor() << " energy -> " << pTheProtein->intraEnergy() << endl;
	cout << "default amberElec::getDielectricConstant() ->" << amberElec::getDielectricConstant() << " energy -> " << pTheProtein->intraEnergy() << endl;
	amberElec::setDielectricConstant(80.0);
	cout << "amberElec::setDielectricConstant(80.0) ->" << amberElec::getDielectricConstant() << " energy -> " << pTheProtein->intraEnergy() << endl;
	cout << "--------  exercise over -------" << endl;


	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	
	amberVDW::setScaleFactor(1.0);

	amberElec::setScaleFactor(1.0);
	amberElec::distanceDependanceOn();

	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,5));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,9));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,12));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,16));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,19));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,23));

	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,5), 11);
	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,9), 11);
	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,12), 11);
	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,16), 11);
	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,19), 11);
	pTheProtein->mutate(0, pTheProtein->getIndexFromResNum(0,23), 11);

	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,5));
	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,9));
	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,12));
	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,16));
	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,19));
	pTheProtein->setOnlyCharged(0, pTheProtein->getIndexFromResNum(0,23));

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheMolecule);
	annealer* pAnneal = new annealer(pEnsemble);
	UInt randomSeed = 10;

    for (int i = 0; i < 20; i ++)
    {
		// BROILING
		amberVDW::setRadiusScaleFactor(1.0);
		amberVDW::setLinearRepulsionDampeningOn();
        pAnneal->run(300.0, 100.0, 2000, randomSeed+(i));

		// SIMMERING
		amberVDW::setLinearRepulsionDampeningOff();
        pAnneal->run(100.0, 10.0, 2000, randomSeed+(i*i));

		// COOLING
		amberVDW::setRadiusScaleFactor(0.95);
        pAnneal->run(10.0, 1.0, 1000, randomSeed+(10*i));

        string outFileName = "coiles_ " ;
        outFileName[7] = 'a' + i;
        outFileName[8] = '\0';
        outFileName += ".pdb";

        pdbWriter(pTheProtein, outFileName);
    }


	
	return 0;
}
 
