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
void setDihedrals (UInt res, UInt maxres, UInt angletype, UInt angletypemax);

int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
        cout << "structShaper <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};

	//--running parameters
	residue::setCutoffDistance(8.0);
	residue::setTemperature(300);
	residue::setElectroSolvationScaleFactor(0.0);
	residue::setHydroSolvationScaleFactor(0.0);
	amberElec::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	srand (time(NULL));

	string outFile, inFile = argv[1];


	///database secondary structure library////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	////Variables////////////////////////////
	vector <double> phisL(3), phisD(3), phisLD(6), psisL(4), psisD(4), psisLD(8);
	vector <double> randCoilL(2), randCoilD(2), coilL(2), coilD(2), alphaL(2), alphaD(2), three10L(2), three10D(2), aTurn1L(2), aTurn2L(2), aTurn3L(2), aTurn1D(2), aTurn2D(2), aTurn3D(2);
	vector <double> ppL(2), ppD(2), pBetaL(2), pBetaD(2), aBetaL(2), aBetaD(2), bTurn1L(2), bTurn2L(2), bTurn3L(2), bTurn1D(2), bTurn2D(2), bTurn3D(2);
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
	//Random Coilbundle->mutateWBC(i, j, P);
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

	////////////////////////////////////////////
	/// denovo idealized heterochiral dihedral set //////
	hetPhis[0] = -60.0, hetPhis[1] = -100.0, hetPhis[2] = -140.0, hetPhis[3] = -180.0, hetPhis[4] = 140.0, hetPhis[5] = 100.0, hetPhis[6] = 60.0;
	hetPsis[0] = -40.0, hetPsis[1] = 0.0, hetPsis[2] = 40.0, hetPsis[3] = 120.0, hetPsis[4] = 160.0, hetPsis[5] = -160.0, hetPsis[6] = -120.0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	UInt interval = 60, count = 0;
	for (int h2 = -120; h2 < 180; h2+=interval)///2
	{
		for (int s2 = -120; s2 < 180; s2+=interval)
		{
			for (int h3 = -120; h3 < 180; h3+=interval) ///3
			{
				for (int s3 = -120; s3 < 180; s3+=interval)
				{
					for (int h4 = -120; h4 < 180; h4+=interval)///4
					{
						for (int s4 = -120; s4 < 180; s4+=interval)
						{
							for (int h5 = -120; h5 < 180; h5+=interval) ///5
							{
								for (int s5 = -120; s5 < 180; s5+=interval)
								{
									/*for (int h6 = -180; h6 < 180; h6+=interval)///6
									{
										for (int s6 = -180; s6 < 180; s6+=interval)
										{
											for (int h7 = -180; h7 < 180; h7+=interval) ///7
											{
												for (int s7 = -180; s7 < 180; s7+=interval)
												{
													for (int h8 = -180; h8 < 180; h8+=interval)///8
													{
														for (int s8 = -180; s8 < 180; s8+=interval)
														{
															for (int h9 = -180; h9 < 180; h9+=interval) ///9
															{
																for (int s9 = -180; s9 < 180; s9+=interval)
																{
																	for (int h10 = -180; h10 < 180; h10+=interval)///10
																	{
																		for (int s10 = -180; s10 < 180; s10+=interval)
																		{
																			for (int h11 = -180; h11 < 180; h11+=interval) ///11
																			{
																				for (int s11 = -180; s11 < 180; s11+=interval)
																				{
																					for (int h12 = -180; h12 < 180; h12+=interval)///12
																					{
																						for (int s12 = -180; s12 < 180; s12+=interval)
																						{*/
																						if (h2 != 0 && h3 != 0 && h4 != 0 && h5 != 0)
																						{
																							count++;
																							PDBInterface* thePDB = new PDBInterface(inFile);
																							ensemble* theEnsemble = thePDB->getEnsemblePointer();
																							molecule* pMol = theEnsemble->getMoleculePointer(0);
																							protein* frame = static_cast<protein*>(pMol);

																							frame->makeResidueSilent(0,0);
																							frame->makeResidueSilent(0,5);

																							frame->setDihedral(0,1,h2,0,0);
																							frame->setDihedral(0,1,s2,1,0);
																							if (h2 > 0){ frame->mutateWBC(0,1,dA);}
																							else { frame->mutateWBC(0,1,A);}

																							frame->setDihedral(0,2,h3,0,0);
																							frame->setDihedral(0,2,s3,1,0);
																							if (h3 > 0){ frame->mutateWBC(0,2,dA);}
																							else { frame->mutateWBC(0,2,A);}

																							frame->setDihedral(0,3,h4,0,0);
																							frame->setDihedral(0,3,s4,1,0);
																							if (h4 > 0){ frame->mutateWBC(0,3,dA);}
																							else { frame->mutateWBC(0,3,A);}

																							frame->setDihedral(0,4,h5,0,0);
																							frame->setDihedral(0,4,s5,1,0);
																							if (h5 > 0){ frame->mutateWBC(0,4,dA);}
																							else { frame->mutateWBC(0,4,A);}

																							/*frame->setDihedral(0,5,h6,0,0);
																							frame->setDihedral(0,5,s6,1,0);
																							if (h6 > 0){ frame->mutateWBC(0,5,dA);}
																							else { frame->mutateWBC(0,5,A);}

																							frame->setDihedral(0,6,h7,0,0);
																							frame->setDihedral(0,6,s7,1,0);
																							if (h7 > 0){ frame->mutateWBC(0,6,dA);}
																							else { frame->mutateWBC(0,6,A);}

																							frame->setDihedral(0,7,h8,0,0);
																							frame->setDihedral(0,7,s8,1,0);
																							if (h8 > 0){ frame->mutateWBC(0,7,dA);}
																							else { frame->mutateWBC(0,7,A);}

																							frame->setDihedral(0,8,h9,0,0);
																							frame->setDihedral(0,8,s9,1,0);
																							if (h9 > 0){ frame->mutateWBC(0,8,dA);}
																							else { frame->mutateWBC(0,8,A);}

																							frame->setDihedral(0,9,h10,0,0);
																							frame->setDihedral(0,9,s10,1,0);
																							if (h10 > 0){ frame->mutateWBC(0,9,dA);}
																							else { frame->mutateWBC(0,9,A);}

																							frame->setDihedral(0,10,h11,0,0);
																							frame->setDihedral(0,10,s11,1,0);
																							if (h11 > 0){ frame->mutateWBC(0,10,dA);}
																							else { frame->mutateWBC(0,10,A);}

																							frame->setDihedral(0,11,h12,0,0);
																							if (h12 > 0){ frame->mutateWBC(0,11,dA);}
																							else { frame->mutateWBC(0,11,A);}*/

																							dblVec C1coords = frame->getCoords(0, 1, "CB");
																							dblVec C2coords = frame->getCoords(0, 4, "CB");
																							double dist1 = CMath::distance(C1coords, C2coords);
																							if (dist1 < 10.6)
																							{
																								/*dblVec C3coords = frame->getCoords(0, 7, "CB");
																								dblVec C4coords = frame->getCoords(0, 10, "CB");
																								double dist2 = CMath::distance(C3coords, C4coords);
																								if (dist2 < 11)
																								{
																									double dist3 = CMath::distance(C1coords, C3coords);
																									if(dist3 < 11)
																									{
																										double dist4 = CMath::distance(C2coords, C4coords);
																										if (dist4 < 11)
																										{
																											double dist5 = CMath::distance(C1coords, C4coords);
																											if (dist5 < 11)
																											{
																												double dist6 = CMath::distance(C2coords, C3coords);
																												if (dist6 < 11)
																												{*/
																													double Energy = frame->protEnergy();
																													if (Energy < 45)
																													{
																														cout << count << " " << Energy << " " << h2 << " " << s2 << " " << h3 << " " << s3 << " " << h4 << " " << s4 << " " << h5 << " " << s5 << endl;//" " << h6 << " " << s6 << " " << h7 << " " << s7 << " " << h8 << " " << s8 << " " << h9 << " " << s9 << " " << h10 << " " << s10 << " " << h11 << " " << s11 << " " << h12 << " " << s12;
																														stringstream convert;
																														string countstr;
																														convert << count, countstr = convert.str();
																														outFile = countstr + ".motif.pdb";
																														pdbWriter(frame, outFile);
																													}
																													//else{ cout << endl; }
																												/*}
																											}
																										}
																									}
																								}*/
																							}
																							delete thePDB;
																						}
																						/*}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}*/
								}
							}
						}
					}
				}
			}
		}
	}
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

void setDihedrals (UInt res, UInt maxres, UInt angletype, UInt angletypemax)
{

	if (res == maxres)
	{
		res = 0;
	}
	if (angletype == 2)
	{
		angletype = 0;
	}
	cout << res << " " << angletype << " || ";
	setDihedrals(res+1, maxres, angletype+1,angletypemax);
}

// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
