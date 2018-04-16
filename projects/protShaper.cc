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
#include <dirent.h>
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
	/*DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");

	vector < double > angles3(6);
	vector < double > angles6(12);
	vector < vector < double > > XCXmotifs;
	vector < vector < double > > XCXXCXmotifs;
	UInt count = 0, interval = 15;
	for (int h1 = -150; h1 < -15; h1+=interval)///1
	{
		for (int s1 = -150; s1 < 165; s1+=interval)
		{
			for (int h2 = 30; h2 < 165; h2+=interval)///2
			{
				for (int s2 = -150; s2 < 165; s2+=interval)
				{
					for (int h3 = -150; h3 < -15; h3+=interval)///3
					{
						for (int s3 = -150; s3 < 165; s3+=interval)
						{
							count++;
							PDBInterface* thePDB = new PDBInterface(inFile);
							ensemble* theEnsemble = thePDB->getEnsemblePointer();
							molecule* pMol = theEnsemble->getMoleculePointer(0);
							protein* frame = static_cast<protein*>(pMol);

							frame->makeResidueSilent(0,0);
							frame->makeResidueSilent(0,4);
							frame->makeResidueSilent(0,5);
							frame->makeResidueSilent(0,6);
							frame->makeResidueSilent(0,7);
							frame->makeResidueSilent(0,8);
							frame->makeResidueSilent(0,9);
							frame->makeResidueSilent(0,10);
							frame->makeResidueSilent(0,11);
							frame->makeResidueSilent(0,12);
							frame->makeResidueSilent(0,13);

							frame->setDihedral(0,1,h1,0,0);
							frame->setDihedral(0,1,s1,1,0);
							if (h1 > 0){ frame->mutateWBC(0,1,dA);}
							else { frame->mutateWBC(0,1,A);}

							frame->setDihedral(0,2,h2,0,0);
							frame->setDihedral(0,2,s2,1,0);
							if (h2 > 0){ frame->mutateWBC(0,2,dA);}
							else { frame->mutateWBC(0,2,A);}

							frame->setDihedral(0,3,h3,0,0);
							frame->setDihedral(0,3,s3,1,0);
							if (h3 > 0){ frame->mutateWBC(0,3,dA);}
							else { frame->mutateWBC(0,3,A);}

							double Energy = frame->protEnergy();
							if (Energy < 3.8)
							{
								cout << count << " " << Energy << " " << h1 << " " << s1 << " " << h2 << " " << s2 << " " << h3 << " " << s3 << endl;
								angles3[0] = h1*-1,angles3[1] = s1*-1,angles3[2] = h2*-1,angles3[3] = s2*-1,angles3[4] = h3*-1,angles3[5] = s3*-1;
								XCXmotifs.push_back(angles3);
								stringstream convert;
								string countstr;
								convert << count, countstr = convert.str();
								outFile = countstr + ".XCX.pdb";
								pdbWriter(frame, outFile);
							}
							delete thePDB;
						}
					}
				}
			}2
		}
	}
	while ((pent=readdir(pdir)))
	{
		string inFrame = pent->d_name;
		if (inFrame.find(".XCX.pdb") != std::string::npos)
		{
			for (UInt i = 0; i < XCXmotifs.size(); i++)
			{
				count++;
				PDBInterface* thePDB = new PDBInterface(inFrame);
				ensemble* theEnsemble = thePDB->getEnsemblePointer();
				molecule* pMol = theEnsemble->getMoleculePointer(0);
				protein* frame = static_cast<protein*>(pMol);

				frame->makeResidueSilent(0,0);
				frame->makeResidueSilent(0,7);
				frame->makeResidueSilent(0,8);
				frame->makeResidueSilent(0,9);
				frame->makeResidueSilent(0,10);
				frame->makeResidueSilent(0,11);
				frame->makeResidueSilent(0,12);
				frame->makeResidueSilent(0,13);

				frame->setDihedral(0,4,XCXmotifs[i][0],0,0);
				frame->setDihedral(0,4,XCXmotifs[i][1],1,0);
				if (XCXmotifs[i][0] > 0){ frame->mutateWBC(0,4,dA);}
				else { frame->mutateWBC(0,4,A);}

				frame->setDihedral(0,5,XCXmotifs[i][2],0,0);
				frame->setDihedral(0,5,XCXmotifs[i][3],1,0);
				if (XCXmotifs[i][2] > 0){ frame->mutateWBC(0,5,dA);}
				else { frame->mutateWBC(0,5,A);}

				frame->setDihedral(0,6,XCXmotifs[i][4],0,0);
				frame->setDihedral(0,6,XCXmotifs[i][5],1,0);
				if (XCXmotifs[i][4] > 0){ frame->mutateWBC(0,6,dA);}
				else { frame->mutateWBC(0,6,A);}

				dblVec C1coords = frame->getCoords(0, 2, "CB");
				dblVec C2coords = frame->getCoords(0, 5, "CB");
				double dist = CMath::distance(C1coords, C2coords);
				if (dist < 10.6)
				{
					double Energy = frame->protEnergy();
					if (Energy < 7.6)
					{
						cout << count << " " << Energy << " " << frame->getAngle(0,1,0) << " " << frame->getAngle(0,1,1) << " " << frame->getAngle(0,2,0) << " " << frame->getAngle(0,2,1) << " " << frame->getAngle(0,3,0) << " " << frame->getAngle(0,3,1) << " " << XCXmotifs[i][0] << " " << XCXmotifs[i][1] << " " << XCXmotifs[i][2] << " " << XCXmotifs[i][3] << " " << XCXmotifs[i][4] << " " << XCXmotifs[i][5] << endl;
						angles6[0] = frame->getAngle(0,1,0),angles6[1] = frame->getAngle(0,1,1),angles6[2] = frame->getAngle(0,2,0),angles6[3] = frame->getAngle(0,2,1),angles6[4] = frame->getAngle(0,3,0),angles6[5] = frame->getAngle(0,3,1),angles6[6] = frame->getAngle(0,4,0),angles6[7] = frame->getAngle(0,4,1),angles6[8] = frame->getAngle(0,5,0),angles6[9] = frame->getAngle(0,5,1),angles6[10] = frame->getAngle(0,6,0),angles6[11] = frame->getAngle(0,6,1);
						XCXXCXmotifs.push_back(angles6);
						stringstream convert;
						string countstr;
						convert << count, countstr = convert.str();
						outFile = countstr + ".XCXXCX.pdb";
						pdbWriter(frame, outFile);
					}
				}
				delete thePDB;
			}
		}
	}
	closedir(pdir);
	while ((pent=readdir(pdir)))
	{
		string inFrame = pent->d_name;
		if (inFrame.find(".XCXXCX.pdh") != std::string::npos)
		{
			PDBInterface* thePDB = new PDBInterface(inFrame);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(pMol);
			angles[0] = frame->getAngle(0,1,0),angles[1] = frame->getAngle(0,1,1),angles[2] = frame->getAngle(0,2,0),angles[3] = frame->getAngle(0,2,1),angles[4] = frame->getAngle(0,3,0),angles[5] = frame->getAngle(0,3,1),angles[6] = frame->getAngle(0,4,0),angles[7] = frame->getAngle(0,4,1),angles[8] = frame->getAngle(0,5,0),angles[9] = frame->getAngle(0,5,1),angles[10] = frame->getAngle(0,6,0),angles[11] = frame->getAngle(0,6,1);
			XCXXCXmotifs.push_back(angles);
		}
	}
	closedir(pdir);
	DIR *pdir2;
	struct dirent *pent2;
	pdir2=opendir(".");
	while ((pent2=readdir(pdir2)))
	{
		string inFrame = pent2->d_name;
		if (inFrame.find(".XCXXCX.pdb") != std::string::npos)
		{
			for (UInt i = 0; i < XCXXCXmotifs.size(); i++)
			{
				count++;
				PDBInterface* thePDB = new PDBInterface(inFrame);
				ensemble* theEnsemble = thePDB->getEnsemblePointer();
				molecule* pMol = theEnsemble->getMoleculePointer(0);
				protein* frame = static_cast<protein*>(pMol);

				frame->makeResidueSilent(0,0);
				frame->makeResidueSilent(0,13);

				frame->setDihedral(0,7,XCXXCXmotifs[i][0],0,0);
				frame->setDihedral(0,7,XCXXCXmotifs[i][1],1,0);
				if (XCXmotifs[i][0] > 0){ frame->mutateWBC(0,7,dA);}
				else { frame->mutateWBC(0,7,A);}

				frame->setDihedral(0,8,XCXXCXmotifs[i][2],0,0);
				frame->setDihedral(0,8,XCXXCXmotifs[i][3],1,0);
				if (XCXmotifs[i][2] > 0){ frame->mutateWBC(0,8,dA);}
				else { frame->mutateWBC(0,8,A);}

				frame->setDihedral(0,9,XCXXCXmotifs[i][4],0,0);
				frame->setDihedral(0,9,XCXXCXmotifs[i][5],1,0);
				if (XCXmotifs[i][4] > 0){ frame->mutateWBC(0,9,dA);}
				else { frame->mutateWBC(0,9,A);}

				frame->setDihedral(0,10,XCXXCXmotifs[i][6],0,0);
				frame->setDihedral(0,10,XCXXCXmotifs[i][7],1,0);
				if (XCXmotifs[i][6] > 0){ frame->mutateWBC(0,10,dA);}
				else { frame->mutateWBC(0,10,A);}

				frame->setDihedral(0,11,XCXXCXmotifs[i][8],0,0);
				frame->setDihedral(0,11,XCXXCXmotifs[i][9],1,0);
				if (XCXmotifs[i][8] > 0){ frame->mutateWBC(0,11,dA);}
				else { frame->mutateWBC(0,11,A);}

				frame->setDihedral(0,12,XCXXCXmotifs[i][10],0,0);
				frame->setDihedral(0,12,XCXXCXmotifs[i][11],1,0);
				if (XCXmotifs[i][10] > 0){ frame->mutateWBC(0,12,dA);}
				else { frame->mutateWBC(0,12,A);}

				dblVec C1coords = frame->getCoords(0, 2, "CB");
				dblVec C2coords = frame->getCoords(0, 5, "CB");
				dblVec C3coords = frame->getCoords(0, 8, "CB");
				dblVec C4coords = frame->getCoords(0, 11, "CB");
				double dist1 = CMath::distance(C1coords, C4coords);
				if (dist1 < 10.6)
				{
					double dist2 = CMath::distance(C2coords, C4coords);
					if (dist2 < 10.6)
					{
						double dist1 = CMath::distance(C3coords, C4coords);
						if (dist1 < 10.6)
						{
							double dist2 = CMath::distance(C1coords, C3coords);
							if (dist2 < 10.6)
							{
								double dist3 = CMath::distance(C2coords, C3coords);
								if (dist3 < 10.6)
								{
									double Energy = frame->protEnergy();
									if (Energy < 20)
									{
										cout << count << " " << Energy << " " << frame->getAngle(0,1,0) << " " << frame->getAngle(0,1,1) << " " << frame->getAngle(0,2,0) << " " << frame->getAngle(0,2,1) << " " << frame->getAngle(0,3,0) << " " << frame->getAngle(0,3,1) << " " << frame->getAngle(0,4,0) << " " << frame->getAngle(0,4,1) << " " << frame->getAngle(0,5,0) << " " << frame->getAngle(0,5,1) << " " << frame->getAngle(0,6,0) << " " << frame->getAngle(0,6,1) << " "<< XCXXCXmotifs[i][0] << " " << XCXXCXmotifs[i][1] << " " << XCXXCXmotifs[i][2] << " " << XCXXCXmotifs[i][3] << " " << XCXXCXmotifs[i][4] << " " << XCXXCXmotifs[i][5] << " "<< XCXXCXmotifs[i][6] << " " << XCXXCXmotifs[i][7] << " " << XCXXCXmotifs[i][8] << " " << XCXXCXmotifs[i][9] << " " << XCXXCXmotifs[i][10] << " " << XCXXCXmotifs[i][11] << endl;
										stringstream convert;
										string countstr;
										convert << count, countstr = convert.str();
										outFile = countstr + ".XCXXCXXCXXCX.pdb";
										pdbWriter(frame, outFile);
									}
								}
							}
						}
					}
				}
				delete thePDB;
			}
		}
	}
	closedir(pdir2);
	DIR *pdir3;
	struct dirent *pent3;
	pdir3=opendir(".");
	while ((pent3=readdir(pdir3)))
	{
		string inFrame = pent3->d_name;
		if (inFrame.find(".XCXXCXXCX.pdb") != std::string::npos)
		{
			for (UInt i = 0; i < XCXmotifs.size(); i++)
			{
				count++;
				PDBInterface* thePDB = new PDBInterface(inFrame);
				ensemble* theEnsemble = thePDB->getEnsemblePointer();
				molecule* pMol = theEnsemble->getMoleculePointer(0);
				protein* frame = static_cast<protein*>(pMol);

				frame->makeResidueSilent(0,0);
				frame->makeResidueSilent(0,13);

				frame->setDihedral(0,10,XCXmotifs[i][0],0,0);
				frame->setDihedral(0,10,XCXmotifs[i][1],1,0);
				if (XCXmotifs[i][0] > 0){ frame->mutateWBC(0,10,dA);}
				else { frame->mutateWBC(0,10,A);}

				frame->setDihedral(0,11,XCXmotifs[i][2],0,0);
				frame->setDihedral(0,11,XCXmotifs[i][3],1,0);
				if (XCXmotifs[i][2] > 0){ frame->mutateWBC(0,11,dA);}
				else { frame->mutateWBC(0,11,A);}

				frame->setDihedral(0,12,XCXmotifs[i][4],0,0);
				frame->setDihedral(0,12,XCXmotifs[i][5],1,0);
				if (XCXmotifs[i][4] > 0){ frame->mutateWBC(0,12,dA);}
				else { frame->mutateWBC(0,12,A);}

				dblVec C1coords = frame->getCoords(0, 2, "CB");
				dblVec C2coords = frame->getCoords(0, 5, "CB");
				dblVec C3coords = frame->getCoords(0, 8, "CB");
				dblVec C4coords = frame->getCoords(0, 11, "CB");
				double dist1 = CMath::distance(C1coords, C4coords);
				if (dist1 < 10.6)
				{
					double dist2 = CMath::distance(C2coords, C4coords);
					if (dist2 < 10.6)
					{
						double dist3 = CMath::distance(C3coords, C4coords);
						if (dist3 < 10.6)
						{
							double Energy = frame->protEnergy();
							if (Energy < 48)
							{
								cout << count << " " << Energy << " " << frame->getAngle(0,1,0) << " " << frame->getAngle(0,1,1) << " " << frame->getAngle(0,2,0) << " " << frame->getAngle(0,2,1) << " " << frame->getAngle(0,3,0) << " " << frame->getAngle(0,3,1) << " " << frame->getAngle(0,4,0) << " " << frame->getAngle(0,4,1) << " " << frame->getAngle(0,5,0) << " " << frame->getAngle(0,5,1) << " " << frame->getAngle(0,6,0) << " " << frame->getAngle(0,6,1) << " " << frame->getAngle(0,7,0) << " " << frame->getAngle(0,7,1) << " " << frame->getAngle(0,8,0) << " " << frame->getAngle(0,8,1) << " " << frame->getAngle(0,9,0) << " " << frame->getAngle(0,9,1) << " " << XCXmotifs[i][0] << " " << XCXmotifs[i][1] << " " << XCXmotifs[i][2] << " " << XCXmotifs[i][3] << " " << XCXmotifs[i][4] << " " << XCXmotifs[i][5] << endl;
								stringstream convert;
								string countstr;
								convert << count, countstr = convert.str();
								outFile = countstr + ".XCXXCXXCXXCX.pdb";
								pdbWriter(frame, outFile);
							}
						}
					}
				}
				delete thePDB;
			}
		}
	}
	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* frame = static_cast<protein*>(pMol);
	double Buffer2 = 0;
	double Buffer1 = 0;
	UInt res;

	res=0;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,psisL[3]+Buffer2,1,0);

	res=1;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[3]+Buffer1,0,0);
	frame->setDihedral(0,res,psisD[3]-Buffer1,1,0);

	res=2;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[6]-Buffer1,0,0);
	frame->setDihedral(0,res,psisL[3]+Buffer1,1,0);

	res=3;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[3]+Buffer2,0,0);
	frame->setDihedral(0,res,psisD[3]-Buffer2,1,0);

	res=4;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[6]-Buffer1,0,0);
	frame->setDihedral(0,res,psisL[3]+Buffer1,1,0);

	res=5;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[3]+Buffer1,0,0);
	frame->setDihedral(0,res,psisD[3]-Buffer1,1,0);

	res=6;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[6]-Buffer2,0,0);
	frame->setDihedral(0,res,psisL[3]+Buffer2,1,0);

	res=7;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[3]+Buffer1,0,0);
	frame->setDihedral(0,res,psisD[3]-Buffer1,1,0);

	res=8;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[6]-Buffer1,0,0);
	frame->setDihedral(0,res,psisL[3]+Buffer1,1,0);

	res=9;
	frame->mutateWBC(0,res, G);
	frame->setDihedral(0,res,hetPsis[3]+Buffer2,0,0);

	pdbWriter(frame, "test.pdb");

	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
    for (int phiR = -180; phiR < 180; phiR++)
	{
		for (int psiR = -180; psiR < 180; psiR++)
		{
            for (int phiL = -180; phiL < 180; phiL++)
			{
				for (int psiL = -180; psiL < 180; psiL++)
				{
					protein* frame = new protein(*start);
					frame->setDihedral(0,0,psiR,1,0);
					frame->setDihedral(0,1,phiL,0,0);
					frame->setDihedral(0,1,psiL,1,0);
					frame->setDihedral(0,2,phiR,0,0);
					frame->setDihedral(0,2,psiR,1,0);
					frame->setDihedral(0,3,phiL,0,0);
					double Energy = frame->protEnergy();
                    if (Energy < 10)
					{
						cout << Energy << " " << phiR << " " << psiR << " " << phiL << " " << psiL << endl;
					}
					delete frame;
                    psiL=psiL+4;
				}
                phiL=phiL+4;
			}
            psiR=psiR+4;
		}
        phiR=phiR+4;
    }*/
    PDBInterface* thePDB = new PDBInterface(inFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* start = static_cast<protein*>(pMol);
    double phi = -80, psi = -20;
    for (UInt h = 0; h < 10; h++)
    {
        protein* bundle = new protein(*start);
        bool even = true;
        for (UInt i = 0; i < bundle->getNumResidues(0); i++)
        {
            if (even)
            {
                bundle->setDihedral(0,i,phi+h,0,0);
                bundle->setDihedral(0,i,psi+h,1,0);
                even = false;
            }
            else
            {
                 bundle->setDihedral(0,i,(phi*-1)+h,0,0);
                 bundle->setDihedral(0,i,(psi*-1)+h,1,0);
                 even = true;
            }
        }
        cout << bundle->protEnergy() << endl;
        stringstream convert;
        string countstr;
        convert << h, countstr = convert.str();
        outFile = countstr + ".SF4trip.pdb";
        pdbWriter(bundle, outFile);
        delete bundle;
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

// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
