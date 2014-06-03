#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"

int main ( void )
{
	string fileName = "1ENH.pdb";
	molecule* pTheMolecule1 = pdbReader(fileName);
	protein* pTheProtein1 = static_cast<protein*>(pTheMolecule1);

	if (pTheProtein1 == 0) return 1;
	
	cout << "number of residues generated:  " << residue::getHowMany() << endl;
	cout << "number of chains generated:  " << chain::getHowMany() << endl;
	cout << "number of residue templates generated: " << residueTemplate::getHowMany() << endl;


	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(1.0);
	amberVDW::setScaleFactor(0.0);
	
  	pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,5));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,8));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,13));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,16));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,17));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,20));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,26));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,31));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,45));
	pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,48));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,49));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,52));
    pTheProtein1->activateForRepacking(0, pTheProtein1->getIndexFromResNum(0,55));

    ensemble *pEnsemble = new ensemble;
    pEnsemble->add(pTheMolecule1);
    annealer* pAnneal = new annealer(pEnsemble);

	for (UInt i = 0; i < 20; i ++)
	{
   		pAnneal->run(400.0, 100.0, 10000, 1776+(i*i));
		pAnneal->run(100.0, 10.0, 10000, 1234+(i*i));
   		pAnneal->run(10.0, 10.0, 1000, 1812+(i*i));
		string outFileName = "hdTest_ " ;
		outFileName[7] = 'a' + i;
		outFileName[8] = '\0';
		outFileName += ".out.pdb";

		pdbWriter(pTheProtein1, outFileName);
	}
	delete pTheProtein1;
	delete pTheMolecule1;

	return 0;
}


