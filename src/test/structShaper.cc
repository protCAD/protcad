//*******************************************************************************************************
//*******************************************************************************************************
//************************************                    ***********************************************
//************************************  structShaper 1.0  ***********************************************
//************************************                    ***********************************************
//*******************************************************************************************************
//**********************  -modify phi psi according to desired structure- *******************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will output modified structure.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=3)
	{
		cout << "structShaper <inFile.pdb> <outFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));

	//--Initialize variables for loop
	stringstream convertphi, convertpsi;
	string inFile = argv[1], outFile = argv[2];
	string phistr, psistr;
    UInt chainNum, resNum, restype;
    double phi, psi;


	//secondary structure library////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	////Variables////////////////////////////
	vector <double> phisL(3), phisD(3), phisLD(6), psisL(4), psisD(4), psisLD(8);
	vector <double> randCoilL(2), randCoilD(2), coilL(2), coilD(2), alphaL(2), alphaD(2), three10L(2), three10D(2), aTurn1L(2), aTurn2L(2), aTurn3L(2), aTurn1D(2), aTurn2D(2), aTurn3D(2);
	vector <double> ppL(2), ppD(2), pBetaL(2), pBetaD(2), aBetaL(2), aBetaD(2), bTurn1L(2), bTurn2L(2), bTurn3L(2), bTurn1D(2), bTurn2D(2), bTurn3D(2);
	vector < vector <double> > Lfolds(12);
	vector < vector <double> > Dfolds(12);

	////Angles///////////////////////////////
	phisL[0] = -63.68, phisL[1] = -99.18, phisL[2] = -139.84;
	psisL[0] = -41.19, psisL[1] = 4.09, psisL[2] = 117.23, psisL[3] = 155.22;
	phisD[0] = (phisL[0]*-1), phisD[1] = (phisL[1]*-1), phisD[2] = (phisL[2]*-1);
	psisD[0] = (psisL[0]*-1), psisD[1] = (psisL[1]*-1), psisD[2] = (psisL[2]*-1), psisD[3] = (psisL[3]*-1);
    phisLD[0] = phisL[0], phisLD[1] = phisD[0], phisLD[2] = phisL[1], phisLD[3] = phisD[1], phisLD[4] = phisL[2], phisLD[5] = phisD[2];
	psisLD[0] = psisL[0], psisLD[1] = psisD[0], psisLD[2] = psisL[1], psisLD[3] = psisD[1], psisLD[4] = psisL[2], psisLD[5] = psisD[2], psisLD[6] = psisL[3], psisLD[7] = psisD[3];

	////L Motifs////////////////////////////
	//Random Coil
	randCoilL[0] = -81.43, randCoilL[1] = 79.65;

	//Alpha-type
	//unfolded
	coilL[0] = phisL[1], coilL[1] = psisL[1];
	//folded
	alphaL[0] = phisL[0], alphaL[1] = psisL[0];
	three10L[0] = phisL[0], three10L[1] = psisL[1];
	//turns
	aTurn1L[0] = phisL[1], aTurn1L[1] = psisL[0];
	aTurn2L[0] = phisL[2], aTurn2L[1] = psisL[0];
	aTurn3L[0] = phisL[2], aTurn3L[1] = psisL[1];

	////Beta-type motifs
	//unfolded
	ppL[0] = phisL[0], ppL[1] = psisL[3];
	//folded
	pBetaL[0] = phisL[1], pBetaL[1] = psisL[2];
	aBetaL[0] = phisL[2], aBetaL[1] = psisL[3];
	//turns
	bTurn1L[0] = phisL[0], bTurn1L[1] = psisL[2];
	bTurn2L[0] = phisL[1], bTurn2L[1] = psisL[3];
	bTurn3L[0] = phisL[2], bTurn3L[1] = psisL[2];

	////D Motifs//////////////////////////////
	//Random Coil
	randCoilD[0] = 81.43, randCoilL[1] = -79.65;

	//Alpha-type
	//unfolded
	coilD[0] = phisD[1], coilD[1] = psisD[1];
	//folded
	alphaD[0] = phisD[0], alphaD[1] = psisD[0];
	three10D[0] = phisD[0], three10D[1] = psisD[1];
	//turns
	aTurn1D[0] = phisD[1], aTurn1D[1] = psisD[0];
	aTurn2D[0] = phisD[2], aTurn2D[1] = psisD[0];
	aTurn3D[0] = phisD[2], aTurn3D[1] = psisD[1];

	//// Beta-type motifs
	//unfolded
	ppD[0] = phisD[0], ppD[1] = psisD[3];
	//folded
	pBetaD[0] = phisD[1], pBetaD[1] = psisD[2];
	aBetaD[0] = phisD[2], aBetaD[1] = psisD[3];
	//turns
	bTurn1D[0] = phisD[0], bTurn1D[1] = psisD[2];
	bTurn2D[0] = phisD[1], bTurn2D[1] = psisD[3];
	bTurn3D[0] = phisD[2], bTurn3D[1] = psisD[2];

	////L folds/////////////////////////////////
	Lfolds[0] = alphaL, Lfolds[1] = ppL, Lfolds[2] = pBetaL, Lfolds[3] = aBetaL, Lfolds[4] = coilL, Lfolds[5] = bTurn2L;
	Lfolds[6] = bTurn3L, Lfolds[7] = bTurn1L, Lfolds[8] = aTurn1L, Lfolds[9] = three10L, Lfolds[10] = aTurn2L, Lfolds[11] = aTurn2L;

	////D folds/////////////////////////////////
	Dfolds[0] = alphaD, Dfolds[1] = ppD, Dfolds[2] = pBetaD, Dfolds[3] = aBetaD, Dfolds[4] = coilD, Dfolds[5] = bTurn2D;
	Dfolds[6] = bTurn3D, Dfolds[7] = bTurn1D, Dfolds[8] = aTurn1D, Dfolds[9] = three10D, Dfolds[10] = aTurn2D, Dfolds[11] = aTurn2D;

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	//--loop
	PDBInterface* theFramePDB = new PDBInterface(inFile);
	ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
	molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
	protein* frame = static_cast<protein*>(frameMol);
	frame->silenceMessages();
	chainNum = frame->getNumChains();

	for (UInt i = 0; i < chainNum; i++)
	{
		resNum = frame->getNumResidues(i);
		for (UInt j = 0; j < resNum; j++)
		{

            if (j >= 22)
			{

                phi = frame->getPhi(i,j);
                frame->setDihedral(i, j, (phi*-1), 0, 0);
                psi = frame->getPsi(i,j);
                frame->setDihedral(i, j, (psi*-1), 1, 0);
                restype = frame->getTypeFromResNum(i,j);
                frame->activateForRepacking(i, j);
                frame->mutateWBC(i, j, (restype+27));
                //randomizeSideChain(frame, i, j);
			}

		}
	}
	pdbWriter(frame, outFile);

//--Print end and write a pdb file--------------------------------------------------------------
	cout << endl << "Structure reshaped!!" << endl << endl;
	return 0;
}

