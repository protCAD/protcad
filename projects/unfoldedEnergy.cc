//*******************************************************************************************************
//*******************************************************************************************************
//************************************                      *********************************************
//************************************  unfoldedEnergy 1.0  *********************************************
//************************************                      *********************************************
//*******************************************************************************************************
//**************************** -find weighted unfolded state energy- ************************************
//*******************************************************************************************************

/////// Just specify a infile and it will result in monomeric folds

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

UInt getFoldingPosition(protein* _prot);

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "unfoldedEnergy <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
    residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

    ///secondary structure library///////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

    //Variables////////////////////////////
	vector <double> phisL(3), phisD(3), phisLD(6), psisL(4), psisD(4), psisLD(8);
	vector <double> randCoilL(2), randCoilD(2), coilL(2), coilD(2), alphaL(2), alphaD(2), three10L(2), three10D(2), aTurn1L(2), aTurn2L(2), aTurn3L(2), aTurn1D(2), aTurn2D(2), aTurn3D(2);
	vector <double> ppL(2), ppD(2), pBetaL(2), pBetaD(2), aBetaL(2), aBetaD(2), bTurn1L(2), bTurn2L(2), bTurn3L(2), bTurn1D(2), bTurn2D(2), bTurn3D(2);
    vector < vector <double> > Lfolds(5);
    vector < vector <double> > Dfolds(5);

    //Angles///////////////////////////////
	phisL[0] = -63.68, phisL[1] = -99.18, phisL[2] = -139.84;
    psisL[0] = -41.19, psisL[1] = 4.09, psisL[2] = 117.23, psisL[3] = 155.22;
	phisD[0] = (phisL[0]*-1), phisD[1] = (phisL[1]*-1), phisD[2] = (phisL[2]*-1);
	psisD[0] = (psisL[0]*-1), psisD[1] = (psisL[1]*-1), psisD[2] = (psisL[2]*-1), psisD[3] = (psisL[3]*-1);
	phisLD[0] = phisL[0], phisLD[1] = phisD[0], phisLD[2] = phisL[1], phisLD[3] = phisD[1], phisLD[4] = phisL[2], phisLD[5] = phisD[2];
	psisLD[0] = psisL[0], psisLD[1] = psisD[0], psisLD[2] = psisL[1], psisLD[3] = psisD[1], psisLD[4] = psisL[2], psisLD[5] = psisD[2], psisLD[6] = psisL[3], psisLD[7] = psisD[3];

    //L Motifs////////////////////////////
	//Random Coil
	randCoilL[0] = -100.9, randCoilL[1] = 58.8;

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

    //Beta-type motifs
	//unfolded
	ppL[0] = phisL[0], ppL[1] = psisL[3];
	//folded
	pBetaL[0] = phisL[1], pBetaL[1] = psisL[2];
	aBetaL[0] = phisL[2], aBetaL[1] = psisL[3];
	//turns
	bTurn1L[0] = phisL[0], bTurn1L[1] = psisL[2];
	bTurn2L[0] = phisL[1], bTurn2L[1] = psisL[3];
	bTurn3L[0] = phisL[2], bTurn3L[1] = psisL[2];

    //D Motifs//////////////////////////////
	//Random Coil
	randCoilD[0] = 100.9, randCoilL[1] = -58.8;

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

    //Beta-type motifs
	//unfolded
	ppD[0] = phisD[0], ppD[1] = psisD[3];
	//folded
	pBetaD[0] = phisD[1], pBetaD[1] = psisD[2];
	aBetaD[0] = phisD[2], aBetaD[1] = psisD[3];
	//turns
	bTurn1D[0] = phisD[0], bTurn1D[1] = psisD[2];
	bTurn2D[0] = phisD[1], bTurn2D[1] = psisD[3];
	bTurn3D[0] = phisD[2], bTurn3D[1] = psisD[2];

    //L folds/////////////////////////////////
    Lfolds[0] = alphaL, Lfolds[1] = ppL, Lfolds[2] = pBetaL, Lfolds[3] = aBetaL, Lfolds[4] = coilL;

    //D folds/////////////////////////////////
    Dfolds[0] = alphaD, Dfolds[1] = ppD, Dfolds[2] = pBetaD, Dfolds[3] = aBetaD, Dfolds[4] = coilD;

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	//--Initialize variables for loop
    string inFile = argv[1], startstr, endstr, outFile;
    UInt resNum, foldPosition, restype, DorL, NCorCN, plateau, nobetter, buffer = 0;
    string tempModel = "temp.pdb";
    double newphi, newpsi, randangle, pastEnergy, Energy, bestEnergy, finalEnergy, foldedEnergy, worstFinalEnergy = 1E10;

	//--loop
    for (int a = 0; a < 10; a++)
	{
		PDBInterface* thePDB = new PDBInterface(inFile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
        resNum = bundle->getNumResidues(0);

        if (a == 0)
        {
            foldedEnergy = bundle->intraSoluteEnergy(true);
        }

        // initial unfold
        for (UInt i = 0; i < resNum; i++)
        {
            restype = bundle->getTypeFromResNum(0, i);
            if (restype < 26)
            {
                newphi = randCoilL[0];
                newpsi = randCoilL[1];
            }
            if (restype > 26)
            {
                newphi = randCoilD[0];
                newpsi = randCoilD[1];
            }
            if (restype == 26)
            {
                DorL = rand() % 1;
                if (DorL == 0)
                {
                    newphi = randCoilL[0];
                    newpsi = randCoilL[1];
                }
                else
                {
                    newphi = randCoilD[0];
                    newpsi = randCoilD[1];
                }
            }
            if (i == 0)
            {
                bundle->setDihedral(0, i, newphi, 0, 0);
            }
            if (i == (resNum-1))
            {
                bundle->setDihedral(0, i, newpsi, 1, 0);
            }
            if (i != 0 && i != (resNum-1))
            {
                bundle->setDihedral(0, i, newphi, 0, 0);
                bundle->setDihedral(0, i, newpsi, 1, 0);
            }
        }
        pdbWriter(bundle, tempModel);

        //starting point
        pastEnergy = bundle->intraSoluteEnergy(true), resNum = bundle->getNumResidues(0), nobetter = 0, bestEnergy = 1E10, plateau = resNum*20;
        foldPosition = getFoldingPosition(bundle);
		dblVec phis(resNum), psis(resNum);
		for (UInt i = 0; i < resNum; i++)
		{
			phis[i] = bundle->getPhi(0,i);
			psis[i] = bundle->getPsi(0,i);
		}
		
		//--Run simulation loop till grand minimum---------------------------------------------------
		do
		{  
			PDBInterface* thePDB = new PDBInterface(inFile);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* bundle = static_cast<protein*>(pMol);
            nobetter++, buffer++;

			//generate random angles---------------------------
            randangle = rand() % Lfolds.size();

			restype = bundle->getTypeFromResNum(0, foldPosition);
            if (restype < 26)
            {
                newphi = Lfolds[randangle][0];
                newpsi = Lfolds[randangle][1];
            }
            if (restype > 26)
            {
                newphi = Dfolds[randangle][0];
                newpsi = Dfolds[randangle][1];
            }
            if (restype == 26)
            {
                DorL = rand() % 1;
                if (DorL == 0)
                {
                    newphi = Lfolds[randangle][0];
                    newpsi = Lfolds[randangle][1];
                }
                else
                {
                    newphi = Dfolds[randangle][0];
                    newpsi = Dfolds[randangle][1];
                }
            }

            //unfold up down the chain
            NCorCN = rand() % 1;
            if (NCorCN == 0)
            {
                //NtoC fold
                for (UInt i = 0; i < resNum; i++)
                {
                    //psi only (n-terminus)
                    if (i == 0 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newpsi, 1, 0);
                    }
                    if (i == 0 && i != foldPosition)
                    {
                        bundle->setDihedral(0, i, psis[i], 1, 0);
                    }
                    //phi only (c-terminus)
                    if (i == resNum-1 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newphi, 0, 0);
                    }
                    if (i == resNum-1 && i != foldPosition)
                    {
                        bundle->setDihedral(0, i, phis[i], 0, 0);
                    }
                    //phi and psi (non-termini)
                    if (i != 0 && i != resNum-1 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newphi, 0, 0);
                        bundle->setDihedral(0, i, newpsi, 1, 0);
                    }
                    else
                    {
                        bundle->setDihedral(0, i, phis[i], 0, 0);
                        bundle->setDihedral(0, i, psis[i], 1, 0);
                    }
                }
            }
            else
            {
                for (UInt i = resNum-1; i >= 0; i--)
                {
                    //psi only (n-terminus)
                    if (i == 0 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newpsi, 1, 1);
                    }
                    if (i == 0 && i != foldPosition)
                    {
                        bundle->setDihedral(0, i, psis[i], 1, 1);
                    }
                    //phi only (c-terminus)
                    if (i == resNum-1 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newphi, 0, 1);
                    }
                    if (i == resNum-1 && i != foldPosition)
                    {
                        bundle->setDihedral(0, i, phis[i], 0, 1);
                    }
                    //phi and psi (non-termini)
                    if (i != 0 && i != resNum-1 && i == foldPosition)
                    {
                        bundle->setDihedral(0, i, newphi, 0, 1);
                        bundle->setDihedral(0, i, newpsi, 1, 1);
                    }
                    else
                    {
                        bundle->setDihedral(0, i, phis[i], 0, 1);
                        bundle->setDihedral(0, i, psis[i], 1, 1);
                    }
                }
            }

            //optimize and check Energy------------------------------
            Energy = bundle->intraSoluteEnergy(true);
            if (Energy < (pastEnergy+buffer))
            {
                pastEnergy = Energy, nobetter--, buffer = 0;
                phis[foldPosition] = newphi, psis[foldPosition] = newpsi;

                if (Energy < bestEnergy)
                {
                    bestEnergy = Energy;
                    pdbWriter(bundle, tempModel);
                }
            }
            if (buffer > 0)
            {
                foldPosition = getFoldingPosition(bundle);
            }
            else
            {
                if (NCorCN == 0 && foldPosition < (resNum-1))
                {
                    foldPosition = foldPosition+1;
                }
                else if (NCorCN == 1 && foldPosition > (0))
                {
                    foldPosition = foldPosition-1;
                }
                else
                {
                    foldPosition = getFoldingPosition(bundle);
                }
            }
			delete thePDB;
        }while (nobetter < plateau*1.2);
		
		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
        model->protOptSolvent(50);
        finalEnergy = model->intraSoluteEnergy(true);
        cout << foldedEnergy-finalEnergy << endl;
        if (finalEnergy > worstFinalEnergy)
        {
            worstFinalEnergy = finalEnergy;
            outFile = "unfolded.pdb";
            pdbWriter(model, outFile);
        }
		delete theModelPDB;
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


UInt getFoldingPosition(protein* _prot)
{
	//--get average position energy
    UInt randres, resNum = _prot->getNumResidues(0), position, rescount = 0;
    double posE, totalposE = 0, aveposE;
    vector <double> posEs(resNum);
	_prot->updateDielectrics();

	for (UInt i = 0; i < resNum; i++)
	{
        rescount++;
        posE = _prot->getPositionSoluteEnergy(0, i, false);
        posEs[i] = posE;
        totalposE = totalposE + posE;
	}
    aveposE = totalposE/rescount;

    //--find position with worse than average energy
    posE = -1E10;
    do
    {
        randres = rand() % resNum;
    }while (posEs[randres] <= aveposE);
    position = randres;
	return position;
}
