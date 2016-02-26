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

void buildAntiParallelBetaBarrel (protein* _prot, double _pitch);
void buildAntiParallelHelixDimer (protein* _prot, double _radius, double _phase, double _coil);
void buildAlaCoilHexamer (protein* _prot, double _radius, double _phase);

int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
        cout << "structShaper <inFile.pdb>" << endl;
		exit(1);
	}
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hcd,Pch};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
    residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
    amberElec::setScaleFactor(1.0);
    srand (time(NULL));

	//--Initialize variables for loop
	stringstream convertphi, convertpsi;
    string outFile = "out.pdb", inFile = argv[1];
	string phistr, psistr;


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
    //turns2.5
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
    double radius = 9.9;
    //double coilstart = 90;
    //double coil = coilstart;
    double best = 1E10;
    double phase = 0;
    //double angle = 0;
    int count = 0;//, rotamer;

    //rotamer optimizations
    /*for (UInt k = 0; k < 82; k++)
    {
        //for (UInt m = 0; m < 36; m++)
        //{
            count++;
            PDBInterface* theFramePDB = new PDBInterface(inFile);
            ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
            molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
            protein* frame 6= static_cast<protein*>(frameMol);
            frame->silenceMessages();
            //buildAlaCoilHexamer(frame, radius, phase);

            for (UInt l = 0; l < frame->getNumChains(); l++)
            {
                UIntVec allowedRots = frame->getAllowedRotamers(0, 23, R, 0);
                frame->activateForRepacking(l, 0);
                frame->mutateWBC(l, 0, R);
                frame->setRotamerWBC(l, 0, 0, allowedRots[14]);
                frame->activateForRepacking(l, 35);
                frame->mutateWBC(l, 35, R);
                frame->setRotamerWBC(l, 35, 0, allowedRots[14]);
                //angle = 10*m;
                //frame->setChi(l, 5, 1, 0, angle);
                //frame->setChi(l, 13, 1, 0, angle);
                //frame->setChi(l, 40, 1, 0, angle);
                //frame->setChi(l, 48, 1, 0, angle);
            }
            double Energy = frame->protEnergy();
            cout << Energy << " " << k;
            if (Energy < 0)
            {
                cout << " hit!" << endl;
                best = Energy;
                stringstream convert;
                string countstr;
                convert << count, countstr = convert.str();
                outFile = countstr + ".good.pdb";
                pdbWriter(frame, outFile);
            }
            else
            {
                cout << endl;
            }
            delete theFramePDB;
        //}
    }*/

    //bundle optimizations
    UIntVec allowedRots;
    cout << "iteration Energy radius rotamer chi" << endl;
    for (UInt i=0; i < 50; i++)
    {
        radius = radius + 0.1;
        //for (UInt j=10; j < 30; j++)
        //{
            //for (UInt k = 0; k < 9; k++)
            //{
                //for (UInt m = 0; m < 36; m++)
                //{
                    count++;
                    phase = 22;//, angle = 10*m;//, rotamer = k;
                    PDBInterface* theFramePDB = new PDBInterface(inFile);
                    ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
                    molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
                    protein* frame = static_cast<protein*>(frameMol);
                    frame->silenceMessages();
                    buildAlaCoilHexamer(frame, radius, phase);
                    for (UInt l = 0; l < frame->getNumChains(); l++)
                    {
                        allowedRots = frame->getAllowedRotamers(0, 0, R, 0);
                        frame->activateForRepacking(l, 0);
                        frame->mutateWBC(l, 0, R);
                        frame->setRotamerWBC(l, 0, 0, allowedRots[14]);
                        frame->activateForRepacking(l, 35);
                        frame->mutateWBC(l, 35, R);
                        frame->setRotamerWBC(l, 35, 0, allowedRots[14]);
                        frame->activateForRepacking(l, 23);
                        frame->mutateWBC(l, 23, R);
                        frame->setRotamerWBC(l, 23, 0, allowedRots[41]);
                        frame->activateForRepacking(l, 58);
                        frame->mutateWBC(l, 58, R);
                        frame->setRotamerWBC(l, 58, 0, allowedRots[41]);

                        allowedRots = frame->getAllowedRotamers(0, 5, Hcd, 0);
                        frame->activateForRepacking(l, 5);
                        frame->mutateWBC(l, 5, Hcd);
                        frame->setRotamerWBC(l, 5, 0, allowedRots[6]);
                        frame->activateForRepacking(l, 13);
                        frame->mutateWBC(l, 13, Hcd);
                        frame->setRotamerWBC(l, 13, 0, allowedRots[7]);
                        frame->activateForRepacking(l, 40);
                        frame->mutateWBC(l, 40, Hcd);
                        frame->setRotamerWBC(l, 40, 0, allowedRots[6]);
                        frame->activateForRepacking(l, 48);
                        frame->mutateWBC(l, 48, Hcd);
                        frame->setRotamerWBC(l, 48, 0, allowedRots[7]);
                        frame->setChi(l, 5, 1, 0, 300);//300
                        frame->setChi(l, 13, 1, 0, 170);
                        frame->setChi(l, 40, 1, 0, 300);
                        frame->setChi(l, 48, 1, 0, 170);
                    }
                    double Energy = frame->protEnergy();
                    cout << count << " " << Energy << " " << radius;
                    if (Energy < best)
                    {
                        cout << " hit!" << endl;
                        best = Energy;
                        stringstream convert;
                        string countstr;
                        convert << count, countstr = convert.str();
                        outFile = countstr + ".good.pdb";
                        pdbWriter(frame, outFile);
                    }
                    else
                    {
                        cout << endl;
                    }
                    delete theFramePDB;
                //}
            //}
        //}
    }


//--Print end and write a pdb file--------------------------------------------------------------
    cout << endl << "Structure reshaped!! best= " << best << endl << endl;
	return 0;
}
void buildAlaCoilHexamer (protein* _prot, double _radius, double _phase)
{
    _prot->rotate(Z_axis, _phase);
    _prot->translate(0.0, _radius, 0.0);
    _prot->rotate(1, Z_axis, 120);
    _prot->rotate(2, Z_axis, 240);
    return;
}

void buildAntiParallelHelixDimer (protein* _prot, double _radius, double _phase, double _coil)
{ 
    _prot->rotate(Z_axis, _phase);
    _prot->rotate(1, Y_axis, 180);
    _prot->translate(0.0, _radius, 0.0);
    _prot->rotate(1, Z_axis, 180);
    _prot->coilcoil(_coil);
    return;
}

void buildAntiParallelBetaBarrel (protein* _prot, double _pitch)
{
    UInt chainNum = _prot->getNumChains();
    UInt resNum;
    //double angle = 0.0, rot1 = 0.0, rot2 = 0.0;
    double dist = 0.0;
    //_prot->coilcoil(_pitch);
    for (UInt i = 0; i < chainNum; i++)
    {
        //--mod structure
        resNum = _prot->getNumResidues(i);
        for (UInt j = 0; j < resNum; j++)
        {
            _prot->setDihedral(i, j, -135, 0, 0);
            _prot->setDihedral(i, j, 135, 1, 0);
        }

        _prot->translateChain(i, 5, 2, dist);
        //_prot->rotateChain(i, Z_axis, 60);
        //_prot->rotateChain(i, Y_axis, 30);
        //_prot->rotateChain(i, X_axis, 20);
        //_prot->translateChain(i, cos(angle)*12, sin(angle)*12, dist);

        /*--set params for next chain
        if (rot1 == 0.0)
        {
            rot1 = 0.0;
            rot2 = 0.0;
        }
        else
        {
            rot1 = 0.0;
            rot2 = 0.0;
        }
        angle = angle+29.5;*/

    }

    return;
}
// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
