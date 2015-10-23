//*******************************************************************************************************
//*******************************************************************************************************
//************************************                    ***********************************************
//************************************   protFolder 2.0   ***********************************************
//************************************                    ***********************************************
//*******************************************************************************************************
//******************** -sample phi psi space according to energy and converge- **************************
//*******************************************************************************************************

/////// Just specify a infile and it will result in monomeric folds

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

UInt getFoldingPosition(protein* _prot, UInt _nobetter);

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protFolder <inFile.pdb>" << endl;
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

    //--Initialize variables for loop and generate unique seed per thread
    srand (getpid());
    string inFile = argv[1], startstr, endstr, outFile;
    UInt resNum, foldPosition, name, restype, DorL, test = 0, count, plateau, nobetter;
    int direction;
    bool gtest;
	stringstream convert;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
    double newphi = 0.0, newpsi = 0.0, randphi, randpsi, pastEnergy, Energy, bestEnergy;

	//--loop
    for (int a = 1; a < 1000; a++)
	{
		PDBInterface* thePDB = new PDBInterface(inFile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
        pastEnergy = bundle->protEnergy(), resNum = bundle->getNumResidues(0), count = 0, nobetter = 0, bestEnergy = 1E10, plateau = resNum;
        foldPosition = getFoldingPosition(bundle,nobetter+1);
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
            nobetter++, gtest = true;

			//generate random angles---------------------------
			restype = bundle->getTypeFromResNum(0, foldPosition);
			randphi = rand() % phisL.size();
			randpsi = rand() % psisL.size();
            if (restype < G)
			{
				newphi = phisL[randphi];
				newpsi = psisL[randpsi];
			}
            if (restype > G)
			{
				newphi = phisD[randphi];
				newpsi = psisD[randpsi];
			}
            if (restype == G)
			{
				DorL = rand() % 1;
				if (DorL == 0)
				{
					newphi = phisL[randphi];
					newpsi = psisL[randpsi];
				}
				else
				{
					newphi = phisD[randphi];
					newpsi = psisD[randpsi];
				}
			}

			//fold and continuously fold loop----------------------
			do
			{
				//NtoC fold
				test++;
				for (UInt i = 0; i < resNum; i++)
				{
                    bundle->setMoved(0,i,1);
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

				//optimize and check Energy------------------------------
                bundle->protOpt(false);
                Energy = bundle->protEnergy();
                if (Energy < (pastEnergy + (nobetter*nobetter)) && (Energy < (pastEnergy-.5) || Energy > (pastEnergy+.5)))
				{
                    //cout << Energy << " " << nobetter << endl;
                    pastEnergy = Energy, test = 0, count++;
                    if (nobetter > 0) { nobetter--;
                    }
                    else{ nobetter = 0;
                    }
					phis[foldPosition] = newphi, psis[foldPosition] = newpsi;

					if (Energy < bestEnergy)
					{
						bestEnergy = Energy;
						pdbWriter(bundle, tempModel);
					}

					//check folding status---------------------------------
					restype = bundle->getTypeFromResNum(0, foldPosition);
                    do
                    { direction = ((rand() % 3) -1);
                    } while (direction == 0);

					//L chirality
                    if (newphi < 0 && foldPosition != 0 && foldPosition != resNum-1 && foldPosition != 0)
					{
						//alpha
						if (newpsi < phisD[0])
						{
							if (newphi == phisL[0] && newpsi == psisL[0])
							{
                                foldPosition = foldPosition+direction;
							}
							else
							{
								newphi = phisL[0];
								newpsi = psisL[0];
							}
						}
						//beta
						if (newpsi > phisD[0])
						{
							if ((newphi == phisL[1] && newpsi == psisL[2]) || (newphi == phisL[2] && newpsi == psisL[3]))
							{
                                foldPosition = foldPosition+direction;
							}
							else
							{
								newphi = phisL[2];
								newpsi = psisL[3];
							}
						}
					}
					//D chirality
                    if (newphi > 0 && foldPosition != 0 && foldPosition != resNum-1 && foldPosition != 0)
					{
						//alpha
						if (newpsi > phisL[0])
						{
							if (newphi == phisD[0] && newpsi == psisD[0])
							{
                                foldPosition = foldPosition+direction;
							}
							else
							{
								newphi = phisD[0];
								newpsi = psisD[0];
							}
						}
						//beta
						if (newpsi < phisL[0])
						{
							if ((newphi == phisD[1] && newpsi == psisD[2]) || (newphi == phisD[2] && newpsi == psisD[3]))
							{
                                foldPosition = foldPosition+direction;
							}
							else
							{
								newphi = phisD[2];
								newpsi = psisD[3];
							}
						}
					}
                }
                else if(restype == G && gtest)
                {
                    test = 0;
                    gtest = false;
                    if (DorL == 1)
                    {
                        newphi = phisL[randphi];
                        newpsi = psisL[randpsi];
                    }
                    else
                    {
                        newphi = phisD[randphi];
                        newpsi = psisD[randpsi];
                    }
                }
			} while (test == 0);
            foldPosition = getFoldingPosition(bundle,nobetter);
			delete thePDB;
		}while (nobetter < plateau);
		
		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		name = rand() % 1000000;
        cout << name << " " << model->protEnergy() << endl;
		stringstream convert; 
		convert << name, endstr = convert.str();
		outFile = endstr + "_model.pdb";
		pdbWriter(model, outFile);
		delete theModelPDB;
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


UInt getFoldingPosition(protein* _prot, UInt _nobetter)
{
    //--get median residue energy
    UInt randres, resNum = _prot->getNumResidues(0);
    double posE, medE;
    medE = _prot->getMedianResEnergy();

    //--find random position with worse than median energy
	posE = -1E10;
	do
	{
		randres = rand() % resNum;
        posE = _prot->resEnergy(0,randres);
    }while (posE < (medE/(_nobetter+1)));
    return randres;
}
