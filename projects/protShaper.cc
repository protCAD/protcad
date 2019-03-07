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

void buildSymmetricOligamer (protein* _prot, bool antiParallel, double _radius, double _coil, double _phaseoffset1, double _phaseoffset2, double _offset);
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
	
	/*PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* frame = static_cast<protein*>(pMol);
	double phi = -60;
	double psi = 140;
	bool odd = false;
	
	for (UInt i = 0; i < frame->getNumResidues(0); i++)
	{
		if (odd){
			frame->setDihedral(0,i,phi+15,0,0);
			frame->setDihedral(0,i,psi+15,1,0);
			frame->mutateWBC(0,i,P);
			odd  = false;
		}
		else{
			frame->setDihedral(0,i,phi*-1+15,0,0);
			frame->setDihedral(0,i,psi*-1+15,1,0);
			frame->mutateWBC(0,i,dP);
			odd = true;
		}
	}
	outFile = inFile;
	pdbWriter(frame, outFile);
	
	
	double delta = 20;
	UInt resIndex = 4;
	
	double phi = frame->getAngle(0,resIndex,0);
	double psi = frame->getAngle(0,resIndex,1);

	
	frame->setDihedral(0,resIndex,psi+delta,1,1);
	vector <double> rpts(4);
	for (UInt i = resIndex+1; i < resIndex+4; i++)
	{
		rpts[i-(resIndex+1)]= frame->getResiduesPerTurn(0,i)/2;
	}
	double aveRPT;
	for (UInt i = 0; i < rpts.size(); i++)
	{
		double rptsum = 0;
		for (UInt j = 0; j <= i; j++)
		{rptsum += rpts[j];}
		double rptsAve = rptsum/(i+1);
		cout << i+1 << " " << rpts[i] << endl;
		if (i+1 > rptsAve){aveRPT = rptsAve; break;}
	}
	cout << aveRPT << endl;
	UInt rptI = aveRPT;
	rptI = rptI+1;
	cout << rptI << endl;
	double fractionRPT = rptI-aveRPT;
	cout << fractionRPT << endl;
	double newDelta = delta-((fractionRPT/aveRPT)*delta);
	cout << newDelta << endl;
	double psi2 = frame->getAngle(0,resIndex+rptI,1);
	frame->setDihedral(0,resIndex+rptI,psi2+newDelta,1,0);
	//frame->setDihedral(0,5,psi-delta,1,1);
	
    pdbWriter(frame, "test_local.pdb");

	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
	//#pragma omp parallel for
	UInt count = 0;
	int phi, psi;
	for (phi = -180; phi < 180; phi++)
	{
		for (psi = -180; psi < 180; psi++)
		{
			count++;
			protein* frame = new protein(*start);
			frame->setDihedral(0,23,phi,0,0);
			frame->setDihedral(0,23,psi,1,0);
			//double Energy = frame->protEnergy();
			//cout << Energy << " " << count << endl;
			stringstream convert;
			string countstr;
			convert << count, countstr = convert.str();
			outFile = countstr + "_helixangle.pdb";
			pdbWriter(frame, outFile);
			psi = psi+29;
			delete frame;
		}
		phi = phi+29;
	}

    /*PDBInterface* thePDB = new PDBInterface(inFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(pMol);
    UInt count = 0;
    UInt res = 7,res2;
    UInt chain = 0;
    UInt foldD = 0;
    double sPhi, sPsi, dihedralD = 20;
    for (UInt i = 0; i < 14; i++)
    {
		sPhi = prot->getPhi(chain,res+i);
		sPsi = prot->getPsi(chain,res+i);
		if (i < 7){
			prot->setDihedral(chain,res+i,sPhi+dihedralD,0,0);
			prot->setDihedral(chain,res+i,sPsi-dihedralD,1,0);
		}
		else{
			prot->setDihedral(chain,res+i,sPhi+dihedralD*-1,0,0);
			prot->setDihedral(chain,res+i,sPsi-dihedralD*-1,1,0);
		}
		stringstream convert;
        string countstr;
        convert << i, countstr = convert.str();
        outFile = countstr + "_helixangle.pdb";
        pdbWriter(prot, outFile);
    }*/
   
    
    
    
    
    double radius;
    int count=0;
	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
    for (int d = 280; d < 290; d++)
    {
        for (int p = 340; p < 350; p++)
        {
			for (int r = 0; r < 5; r++)
			{
				for (int c = 170; c < 190; c++)
				{
					count++;
					radius = 7.5 + (r*0.1);
					protein* bundle = new protein(*start);
					buildSymmetricOligamer (bundle, true, radius, c, p, d, 5);
					cout << count << " " << radius << " " << d << " " << p << " " << bundle->getNumHardClashes() << endl;
					stringstream convert;
					string countstr;
					convert << count, countstr = convert.str();
					outFile = countstr + "hethelix.pdb";
					pdbWriter(bundle, outFile);
					delete bundle;
				}
			}
            //p=p+9;
        }
        //d=d+9;
    }
    return 0;
}
void buildSymmetricOligamer (protein* _prot, bool antiParallel, double _radius, double _coil, double _phaseoffset1, double _phaseoffset2, double _offset)
{
	UInt numChains = _prot->getNumChains();
	double radial = 0.0;
	double radialInterval = 360/numChains;
	
	//--sample phase offset starting point and fix as parallel or antiparallel
	bool odd = false;
	for (UInt i = 0; i < numChains; i++)
	{
		if (antiParallel)
		{
			if (odd)
			{
				_prot->rotate(i, Y_axis, 180);
                _prot->rotate(i,Z_axis, _phaseoffset2);
				odd = false;
			}
			else
			{
				odd = true;
			}
		}
        _prot->rotate(i,Z_axis, _phaseoffset1);
	}
	
	//--translate off z and rotate each chain into position
	odd = false;
	_prot->translate(0.0, _radius, 0.0);
	for (UInt i = 0; i < numChains; i++)
	{
		_prot->rotate(i, Z_axis, radial);
		if (_offset != 0.0)
		{
			if (odd)
			{
				_prot->translate(i, 0.0, 0.0, _offset);
				odd = false;
			}
			else
			{
				odd = true;
			}
		}
		radial += radialInterval;
	}
	
	//--supercoil oligamer
	if (_coil != 0.0){_prot->coilcoil(_coil);}
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
