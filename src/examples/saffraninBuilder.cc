#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

VecDblVec getBackbone(string _backboneFile);
VecDblVec getMCParams(string _mcParmasFile);
VecDblVec monteCarlo(VecDblVec _mcParams, ran _ranNumber, VecDblVec _backboneParams);


enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
protein* prot;


int main (int argc, char* argv[])
{

        // read in PDB file
	string protFile = argv[1];
	PDBInterface* thePDB = new PDBInterface(protFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	prot = static_cast<protein*>(theMol);

        // read in backbone params
        string backboneFile = argv[2];
        VecDblVec backboneParams = getBackbone(backboneFile);

        string mcFile = argv[3];
        VecDblVec mcParams = getMCParams(mcFile);


	// mutate i, i+4 to serine and i+1 to leucine, set rotamers

	// align chains so that center of serine CBs point towards -x and are at z = 0

