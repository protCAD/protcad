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
void buildHelixOligamer (protein* _prot, UInt numChains, bool antiParallel, double _radius, double _phase1, double _phase2, double _coil, double _offset);
double getBackboneHBondEnergy (protein* _prot);

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
    amberElec::setScaleFactor(0.0);
    //amberVDW::setLinearRepulsionDampeningOn();
    //amberVDW::itsRepulsionScaleFactor = 0.9;
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
    vector <double> hetPhis(7);
    vector <double> hetPsis(7);

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
    //unfolded36
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

    ////////////////////////////////////////////
    ////denovo heterochiral dihedral set //////
    hetPhis[0] = -60.0, hetPhis[1] = -100.0, hetPhis[2] = -140.0, hetPhis[3] = -180.0, hetPhis[4] = 140.0, hetPhis[5] = 100.0, hetPhis[6] = 60.0;
    hetPsis[0] = -40.0, hetPsis[1] = 0.0, hetPsis[2] = 40.0, hetPsis[3] = 120.0, hetPsis[4] = 160.0, hetPsis[5] = -160.0, hetPsis[6] = -120.0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

    double radius;
   double best = 1E100;
    double phase1, phase2;
    UInt count = 0;
    double coil;
    double offset = 0.0;
    //rotamer optimizations
    /*PDBInterface* theFramePDB = new PDBInterface(inFile);
    ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
    molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(frameMol);
    double startE = bundle->protEnergy();
    for (UInt f = 33; f < 37; f++)
    {
        for (UInt g = 33; g < 37; g++)
        {
            for (UInt h = 33; h < 37; h++)
            {
                for (UInt k = 33; k < 37; k++)
                {
                    for (UInt i = 0; i < phisLD.size(); i++)
                    {
                        for (UInt j = 0; j < psisLD.size(); j++)
                        {
                            count++;
                            protein* frame = new protein(*bundle);
                            frame->setDihedral(1,f,phisLD[i],0,0);
                            frame->setDihedral(1,f,psisLD[j],1,0);
                            frame->setDihedral(1,g,phisLD[i],0,0);
                            frame->setDihedral(1,g,psisLD[j],1,0);
                            frame->setDihedral(1,h,phisLD[i],0,0);
                            frame->setDihedral(1,h,psisLD[j],1,0);
                            frame->setDihedral(1,k,phisLD[i],0,0);
                            frame->setDihedral(1,k,psisLD[j],1,0);
                            dblVec Ccoords = frame->getCoords(1, 36, "C");
                            dblVec Ncoords = frame->getCoords(2, 0, "N");
                            double dist = CMath::distance(Ccoords, Ncoords);
                            double Energy = 0.0;
                            Energy += -1000/dist;
                            Energy += frame->protEnergy();
                            cout << count << " " << dist << " " << Energy << " " << f << " " << phisLD[i] << " " << psisLD[j] << " " << g << " " << phisLD[i] << " " << psisLD[j] << h << " " << phisLD[i] << " " << psisLD[j];
                            if (Energy < best)
                            {
                                cout << " hit!!!!!!!!" << endl;
                                best = Energy;
                                stringstream convert;
                                string countstr;
                                convert << count, countstr = convert.str();
                                outFile = countstr + ".turn.pdb";
                                pdbWriter(frame, outFile);
                            }
                            else
                            {
                                cout << endl;
                            }
                            delete frame;
                       }
                    }
                }
            }
        }
    }

    ///bundle optimizations/
    //UIntVec allowedRots;*/
    PDBInterface* theFramePDB = new PDBInterface(inFile);
    ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
    molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(frameMol);
     for (int k = 0; k < 10; k++)
     {
        //for (int l = 0; l < 10; l++)
        //{
            for (int n = 0; n < 20; n++)
            {
                for (int m = 0; m < 10; m++)
                {
                    protein* frame = new protein(*bundle);
                    count++;
                    radius = 7.0+ (k*0.1), coil = 0, offset = -10.0+(0.1*m), phase1 = 0, phase2 = 150+n;
                    buildHelixOligamer(frame, 4, true, radius, phase1, phase2, coil, offset);
                    /*dblVec HNcoords = frame->getCoords(0, 18, "NE2");
                    dblVec YOcoords = frame->getCoords(2, 18, "OH");
                    double dist = CMath::distance(HNcoords, YOcoords);
                    double Energy = 0.0;
                    if (dist <= 2.7 && dist >= 2.5)
                    {
                        Energy = -500;
                    }*/
                    double Energy = frame->intraSoluteEnergy(false);
                    cout << count << " " << Energy << " " << radius << " " << offset << " " << phase1 << " " << phase2;
                    if (Energy < best)
                    {
                        cout << " hit!!!!!!!!" << endl;
                        best = Energy;
                        stringstream convert;
                        string countstr;
                        convert << count, countstr = convert.str();
                        outFile = countstr + ".olig.pdb";
                        pdbWriter(frame, outFile);
                    }
                    else
                    {
                        cout << endl;
                    }
                    delete frame;
               }
            }
        //}
    }
    /*UInt res1, res2, rot;
    double chi1, chi2;
    double angle1, angle2;
    cout << "iteration Energy radius angle" << endl;/
    double phase, radius;
            for (int m = 0; m < 100; m++)
            {
               for (int k = 80; k < 110; k++)
               {
                    count++;
                    phase = k, radius = 4+ (m*0.1);
                    PDBInterface* theFramePDB = new PDBInterface(inFile);
                    ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
                    molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
                    protein* frame = static_cast<protein*>(frameMol);
                    frame->silenceMessages();
                    buildAlaCoilTetramer(frame, radius, phase);
                    /for (UInt l = 0; l < frame->getNumChains(); l++)
                    {
                        allowedRots = frame->getAllowedRotamers(0, res1, Hcd, 0);
                        frame->activateForRepacking(l, res1);
                        frame->mutateWBC(l, res1, Hcd);
                        frame->setRotamerWBC(l,frame res1, 0, allowedRots[rot]);
                        frame->setChi(l, res1, 0, 0, chi1);
                        frame->setChi(l, res1, 0, 1, chi2);
                        frame->setChi(l, res1, 1, 0, angle1);
                        frame->activateForRepacking(l, res2);
                        frame->mutateWBC(l, res2, Hcd);
                        frame->setRotamerWBC(l, res2, 0, allowedRots[rot]);
                        frame->setChi(l, res2, 0, 0, chi1);
                        frame->setChi(l, res2, 0, 1, chi2);
                        frame->setChi(l, res2, 1, 0, angle2);
                    }/
                    double Energy = frame->protEnergy();                  
                    cout << count << " " << Energy << " " << radius << " " << phase;
                    if (Energy < best)
                    {
                        best = Energy;
                        cout << " hit!!!!!!!!" << endl;
                        stringstream convert;
                        string countstr;
                        convert << count, countstr = convert.str();
                        outFile = countstr + ".tetramer.pdb";
                        pdbWriter(frame, outFile);
                    }
                    else
                    {
                        cout << endl;
                    }

                    //else
                    //{
                    //    cout << endl;
                    //}
                    delete theFramePDB;
                }
            }
    //}
    //getdihed0rals
    PDBInterface* theFramePDB = new PDBInterface(inFile);
    ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
    molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(frameMol);
    for (UInt a = 0; a < hetPhis.size(); a++)
    {
        for (UInt b = 0; b < hetPsis.size(); b++)
        {
            for (UInt c = 0; c < hetPhis.size(); c++)
            {
                for (UInt d = 0; d < hetPsis.size(); d++)
                {
                    count++;
                    protein* frame = new protein(*bundle);
                    for (UInt i = 0; i < frame->getNumChains(); i++)
                    {
                        bool even = true;
                        for (UInt j = 0; j < frame->getNumResidues(i); j++)
                        {
                            if (j != 0 && j != frame->getNumResidues(i)-1)
                            {
                                if (even)
                                {
                                    frame->setDihedral(i, j, hetPhis[a], 0, 0);
                                    frame->setDihedral(i, j, hetPsis[b], 1, 0);
                                    even = false;
                                }
                                else
                                {
                                    frame->setDihedral(i, j, hetPhis[c], 0, 0);
                                    frame->setDihedral(i, j, hetPsis[d], 1, 0);
                                    even = true;
                                }
                            }
                        }
                        double intraEnergy = frame->protEnergy();
                        double HBondEnergy = getBackboneHBondEnergy(frame);
                        double Energy = intraEnergy+HBondEnergy;
                        if (Energy < -300 && Energy > -500)
                        {
                            stringstream convert;
                            string countstr;
                            convert << count, countstr = convert.str();
                            outFile = countstr + ".struct.pdb";
                            pdbWriter(frame, outFile);
                        }
                        cout << count << " " << Energy << " " << hetPhis[a] << " " << hetPsis[b] << " " << hetPhis[c] << " " << hetPsis[d] <<endl;
                        delete frame;
                    }
                }
            }
        }
    }*/

//--Print end and write a pdb file--------------------------------------------------------------
	return 0;
}
void buildHelixOligamer (protein* _prot, UInt numChains, bool antiParallel, double _radius, double _phase1, double _phase2, double _coil, double _offset)
{
    double rotationInterval = 360/numChains;
    double rotation = 0.0;
    bool odd = false;
    _prot->rotate(0,Z_axis,_phase1);
    _prot->rotate(2,Z_axis,_phase1);
    _prot->rotate(1,Z_axis,_phase2);
    _prot->rotate(3,Z_axis,_phase2);


    if (antiParallel)
    {
        for (UInt i = 0; i < numChains; i++)
        {
            if (odd)
            {
                _prot->rotate(i, Y_axis, 180);
                odd = false;
            }
            else
            {
                odd = true;
            }
        }
        odd = false;
    }
    _prot->translate(0.0, _radius, 0.0);
    for (UInt i = 0; i < numChains; i++)
    {
        _prot->rotate(i, Z_axis, rotation);
        if (odd)
        {
            _prot->translate(i, 0.0, 0.0, _offset);
            odd = false;
        }
        else
        {
            odd = true;
        }
        rotation += rotationInterval;
    }
    //_prot->coilcoil(_coil);
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

double getBackboneHBondEnergy (protein* _prot)
{
    double Energy = 0;
    UInt size = _prot->getNumResidues(0);
    UInt numHbonds = 0;
    vector < dblVec > Ncoords;
    vector < dblVec > Ocoords;
    vector < dblVec > Ccoords;
    vector < dblVec > CAcoords;

    for (UInt i = 0; i < size; i ++)
    {
        dblVec tempcoord =_prot->getCoords(0, i, 0);
        Ncoords.push_back(tempcoord);
        tempcoord =_prot->getCoords(0, i, 3);
        Ocoords.push_back(tempcoord);
        tempcoord =_prot->getCoords(0, i, 1);
        CAcoords.push_back(tempcoord);
        tempcoord = _prot->getCoords(0, i, 2);
        Ccoords.push_back(tempcoord);
    }

    // for first residue = no angular constraint
    for (UInt j = 0; j < Ocoords.size(); j ++)
    {
        dblVec temp = Ncoords[0] - Ocoords[j];
        double distance = sqrt(CMath::dotProduct(temp,temp));
        int resSpace = (0 - j);
        if (distance < 3.2 && abs(resSpace) > 2 )
        {
            numHbonds ++;
        }
    }


    for (UInt i = 1; i < Ncoords.size(); i ++)
    {
        for (UInt j = 0; j < Ocoords.size(); j ++)
        {
            int resSpace = (i - j);
            if ( abs(resSpace) > 2)
            {
                dblVec NO = Ncoords[i] - Ocoords[j];
                double distance = sqrt(CMath::dotProduct(NO,NO));

                dblVec pseudoAtom = (Ccoords[i-1] + CAcoords[i])/2.0;
                dblVec NH = Ncoords[i] - pseudoAtom;

                double magNH = sqrt(CMath::dotProduct(NH,NH));
                double angle = acos( CMath::dotProduct(NO,NH) / (magNH * distance) );
                angle = angle * 180 / 3.14159;
                if (distance < 3.6 && angle > 140)
                {
                    Energy += -40.0/distance;
                }
            }
        }
    }
    return Energy;
}
// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
