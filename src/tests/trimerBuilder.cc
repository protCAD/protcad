#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{

//  STEP 1:  generate bundle
    if (argc != 7)
    {
        cout << "trimerBuilder  <infilename.pdb>  <bundle radius>  <zrot>  <supercoil pitch>  <outfilename.pdb> <yes/no optimize rots>" << endl;
    exit(1);
    }
	 enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    string infile = argv[1];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);

    protein* bundle = static_cast<protein*>(pMol);
    residue::setCutoffDistance(5.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(1.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

    // get parameters from input
    double radius;
    double face;
    double pitch;
    sscanf(argv[2], "%lf", &radius);
	sscanf(argv[3], "%lf", &face);
    sscanf(argv[4], "%lf", &pitch);

	bundle->rotate(Z_axis, face);
	bundle->translate(0,radius,0);
	bundle->rotate(1,Z_axis, 120);
	bundle->rotate(2,Z_axis, 240);


    bundle->coilcoil(pitch);

    string outfile = argv[5];
    bundle->symmetryLinkChainAtoB(1,0);
	bundle->symmetryLinkChainAtoB(2,0);

    int resID[] = {A,A,A,A,A,V,E,A,L,E,K,K,V,A,L,L,E,S,K,V,Q,A,L,E,K,K,V,K,Q,L,L,I,A,V,L,L,L,I,A,V,N,L,I,L,L,I,A,V,A,R,L,R,Y,L,V,G,A,A,A,A};

    int arraySize = sizeof(resID)/sizeof(resID[0]);

    for (int i = 0; i < arraySize; i ++)
    {
        bundle->activateForRepacking(0,i);
        bundle->setCanonicalHelixRotamersOnly(0,i);
        bundle->mutate(0, i, (UInt)resID[i]);
        cout << resID[i] << endl;
    }
	string answer = argv[6];
	if (answer == "yes") bundle->optimizeRotamers();


    pdbWriter(bundle, outfile);

    return 0;
}

