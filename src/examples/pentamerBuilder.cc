#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"


int main (int argc, char* argv[])
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);
    	protein* bundle = static_cast<protein*>(pMol);

	pmf::setScaleFactor(0.0);
    	rotamer::setScaleFactor(1.0);
    	microEnvironment::setScaleFactor(0.0);
    	amberVDW::setScaleFactor(1.0);
    	amberVDW::setRadiusScaleFactor(1.0);
    	amberVDW::setLinearRepulsionDampeningOn();
    	amberElec::setScaleFactor(0.0);
    	solvation::setItsScaleFactor(0.0);

//    	int resID[] = {A,A,A,A,A,I,M,I,I,I,C,C,V,I,L,G,I,I,I,A,C,T,I,G,G,I,F,G,A,A,A,A,A};
    //	int arraySize = sizeof(resID)/sizeof(resID[0]);

	for (UInt i = 1 ; i < 5; i ++)
	{
		bundle->symmetryLinkChainAtoB(i,0);
	}

   // 	for (int i = 0; i < arraySize; i ++)
    //	{
   //     	bundle->activateForRepacking(0,i);
  //      	bundle->setCanonicalHelixRotamersOnly(0,i);
  //      	bundle->mutate(0, i, (UInt)resID[i]);
    //	}

	double radius, phase, coil;

	sscanf(argv[3], "%lf", &radius);
	sscanf(argv[4], "%lf", &phase);
	sscanf(argv[5], "%lf", &coil);

	bundle->silenceMessages();
//	bundle->optimizeRotamers();

	bundle->rotate(Z_axis, phase);	
	bundle->translate(radius, 0.0, 0.0);
	for (UInt i = 1; i < 5; i ++)
	{
		double angle = 72.0 * (double)i;
		bundle->rotate(i, Z_axis, angle); 
	}
	bundle->coilcoil(coil);

	pdbWriter(bundle, argv[2]);
}




