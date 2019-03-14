// repacking of transmembrane dimer with coil face-face interactions 
// made of glycine pairs
// 3/09/2001

#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"

int main (int argc, char* argv[])
{

	string fileName = argv[1];
	molecule* pMol = pdbReader(fileName);
	protein* pTheProtein = static_cast<protein*>(pMol);

	// use only vdW energies
	residue::setCutoffDistance(8.0);
	//residue::setOneFourVDWScaleFactor(0.0);
	//residue::setOneFourAmberElecScaleFactor(0.0);
	pmf::setScaleFactor(0.0);	
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberVDW::setRadiusScaleFactor(1.0);
	amberElec::setScaleFactor(0.0);
	cout << "*****************************************" << endl;
	cout << "           energy parameters             " << endl;
	cout << endl;
	cout << "pmf scale: " << pmf::getScaleFactor() << endl;
	cout << "rotamer scale: " << rotamer::getScaleFactor() << endl;
	cout << "microEnvironment scale: " << microEnvironment::getScaleFactor() << endl;
	cout << "AMBER vdW scale: " << amberVDW::getScaleFactor() << endl;
	cout << "\tradius scale factor: " << amberVDW::getRadiusScaleFactor() << endl;
	cout << "\tlinear repulsion dampening: "; if (amberVDW::linearRepulsionDampening) cout << "on"; else cout << "off"; cout << endl;
 	cout << "AMBER electrostatics scale: " << amberElec::getScaleFactor() << endl;
	cout << "\tdistance dependent dielectric: "; if (amberElec::isDistanceDependanceOn()) cout << "on"; else cout << "off"; cout << endl;
	cout << endl;
	cout << "*****************************************" << endl;

	pTheProtein->symmetryLinkChainAtoB(1,0);
	//pTheProtein->symmetryLinkChainAtoB(2,0);


	// ACTIVATE ALL RESIDUES FOR REPACKING FROM 6-34 EXCEPT FOR GLYCINE
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,6));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,7));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,8));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,9));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,10));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,11));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,13));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,14));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,16));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,17));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,18));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,20));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,21));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,23));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,24));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,25));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,26));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,27));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,28));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,29));
    pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,30));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,31));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,32));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,33));
	pTheProtein->activateForRepacking(0, pTheProtein->getIndexFromResNum(0,34));



	vector <UInt> bannedResList;
	bannedResList.resize(0);

	// BUILD BANNED LIST FOR FLOATED POSITIONS
	// keep ALA
	bannedResList.push_back(1); // ban ARG
	bannedResList.push_back(2); // ban ASN
	bannedResList.push_back(3); // ban ASP
	bannedResList.push_back(4); // ban CYS
	bannedResList.push_back(5); // ban GLN
	bannedResList.push_back(6); // ban GLU
	bannedResList.push_back(7); // ban GLY
	bannedResList.push_back(8); // ban HIS
	// keep ILE
	// keep LEU
	bannedResList.push_back(11); // ban LYS 
	bannedResList.push_back(12); // ban MET
	bannedResList.push_back(13); // ban PHE
	bannedResList.push_back(14); // ban PRO
	bannedResList.push_back(15); // ban SER
	// keep THR
	bannedResList.push_back(17); // ban TRP
	bannedResList.push_back(18); // ban TYR
	// keep VAL

	
	// SET IDENTITES FOR FLOATED POSTIONS
	pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,8), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,9), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,11), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,16), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,18), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,23), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,25), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,26), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,29), bannedResList);
    pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,30), bannedResList);
	pTheProtein->setListNotAllowed(0, pTheProtein->getIndexFromResNum(0,33), bannedResList);

	// MUTATE FIXED POSITIONS TO MS1 SEQUENCE AND HOLD
	pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,6),0);		
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,6));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,7),6);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,7));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,10),9);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,10));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,13),10);  	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,13));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,14),10);  	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,14));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,17),0);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,17));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,20),10);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,20));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,21),9);  	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,21));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,24),9);  	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,24));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,27),0);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,27));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,28),1);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,28));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,31),18);   	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,31));
    pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,32),10);  	
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,32));
	pTheProtein->mutate(0,pTheProtein->getIndexFromResNum(0,34),7);		
	pTheProtein->setOnlyNativeIdentity(0,pTheProtein->getIndexFromResNum(0,34));

	pTheProtein->setCanonicalHelixRotamersOnly(0);

	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,8), 16, 0, 0); // remove rotamer 0 of THR
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,9), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,11), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,16), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,18), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,23), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,25), 16, 0, 0);	
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,26), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,29), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,30), 16, 0, 0);
	pTheProtein->setRotamerNotAllowed(0, pTheProtein->getIndexFromResNum(0,33), 16, 0, 0);

	cout << " rotamer and identities set ... proceeding with annealing." << endl;

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pMol);
	annealer* pAnneal = new annealer(pEnsemble);
    for (int i = 0; i <= 25; i ++)
    {
		cout << "***************** ITERATION -" << i << "- *****************" << endl;
		pAnneal->run(600.0,600.0, 500, i*5); // DEEP FRY
		pAnneal->run(50.0, 10.0, 2000, i);   // BAKE
		pAnneal->run(10.0, 2.0, 3000, 2*i);  // SIMMER
		pAnneal->run(2.0, 1.0, 6000, 3*i);   // COOL
        string outFileName = "repack_ " ;
        outFileName[7] = 'a' + i;
        outFileName[8] = '\0';
        outFileName += ".pdb";

        pdbWriter(pTheProtein, outFileName);
    }
	delete pTheProtein;

	return 0;
}



