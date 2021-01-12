//*******************************************************************************************************
//*******************************************************************************************************
//************************************                    ***********************************************
//************************************  protOligamer 1.0  ***********************************************
//************************************                    ***********************************************
//*******************************************************************************************************
//*****************************  build symetric protein oligamer ****************************************
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

void buildSymmetricOligamer (protein* _prot, bool antiParallel, double _radius, double _phase, double _antiphase, double _coil, double _offset);
enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,SF4,HEM,NI2,CLN,CO2,MG2,OH,OXY,CLD,HIS};
int main (int argc, char* argv[])
{
	//--Running p4arameters
    if (argc !=8)
	{
        cout << "protOligamer <inFile.pdb>" << endl;
		exit(1);
	}
	string outFile, inFile = argv[1];
	string antiparallel = argv[2];
	bool anti = false;
	if (antiparallel == "t"){anti = true;}
	double radius, phase, antiphase, coil, offset;

	sscanf(argv[3], "%lf", &radius);
	sscanf(argv[4], "%lf", &phase);
	sscanf(argv[5], "%lf", &antiphase);
	sscanf(argv[6], "%lf", &coil);
	sscanf(argv[7], "%lf", &offset);
	
	//--running parameters
	residue::setElectroSolvationScaleFactor(0.0);
	residue::setHydroSolvationScaleFactor(0.0);
	residue::setPolarizableElec(false);
	amberElec::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);
	srand (time(NULL));

	
	int count=0;
	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
	start->alignToAxis(Z_axis);
	//for (int c = 90; c < 100; c++)
	//{
		//for (int p = 100; p < 120; p++)
		//{
			//for (int a = 100; a < 120; a++)
			//{
				//for (int r = 0; r < 360; r++)
				//{
					//for (int o = 0; o < 10; o++)
					//{
						count++;
						//double r = 6;
						//int c = 180;
						//double p = 115;
						//double offset = -3.0;
						protein* bundle = new protein(*start);
						//bundle->mutateWBC(0,10,Cf);
						//bundle->mutateWBC(1,10,Cf);
						//bundle->mutateWBC(2,10,Cf);
						//bundle->mutateWBC(3,10,Cf);
						buildSymmetricOligamer (bundle, anti, radius, phase, antiphase, coil, offset);
						double chi1 = -170;
						double chi2 = -170;
						bundle->setChi(0,10,0,0,chi2);
						bundle->setChi(1,10,0,0,chi1);
						bundle->setChi(2,10,0,0,chi2);
						bundle->setChi(3,10,0,0,chi1);
						cout << bundle->getChi(0,10,0,0) << " " << bundle->getChi(1,10,0,0) << " " << bundle->getChi(2,10,0,0) << " " << bundle->getChi(3,10,0,0) << endl;
						UInt clashes = bundle->getNumHardClashes();
						cout << bundle->protEnergy() << endl;
						dblVec S1coords = bundle->getCoords(0, 10, "SG");
						dblVec S2coords = bundle->getCoords(1, 10, "SG");
						dblVec S3coords = bundle->getCoords(2, 10, "SG");
						dblVec S4coords = bundle->getCoords(3, 10, "SG");
						double dist1 = CMath::distance(S1coords, S2coords);
						double dist2 = CMath::distance(S2coords, S3coords);
						double dist3 = CMath::distance(S3coords, S4coords);
						double dist4 = CMath::distance(S4coords, S1coords);
						double dist5 = CMath::distance(S3coords, S1coords);
						double dist6 = CMath::distance(S4coords, S2coords);
						double cutoff1 = 4.0, cutoff2 = 3.0;
						cout << "radius phase coil offset meandist clashes" << endl;
						if ((dist1 < cutoff1 && dist1 > cutoff2) && (dist2 < cutoff1 && dist2 > cutoff2) && (dist3 < cutoff1 && dist3 > cutoff2) && (dist4 < cutoff1 && dist4 > cutoff2) && (dist5 < cutoff1 && dist5 > cutoff2) && (dist6 < cutoff1 && dist6 > cutoff2)){
							cout << radius << " " << phase << " " << antiphase << " " << coil << " " << offset << " " << clashes << " " << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << " " << dist5 << " " << dist6 << " hit!" << endl;
						}
						else{
							cout << radius << " " << phase << " " << antiphase << " " << coil << " " << offset << " " << clashes << " " << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << " " << dist5 << " " << dist6 << endl;
						}
						stringstream convert;
						string countstr;
						convert << count, countstr = convert.str();
						outFile = countstr + "_oligamer.pdb";
						pdbWriter(bundle, outFile);
						delete bundle;
					//}
				//}
				//p=p+4;
		//}
		//c=c+4;
	//}
	return 0;
}


void buildSymmetricOligamer (protein* _prot, bool antiParallel, double _radius, double _phase, double _antiphase, double _coil, double _offset)
{
	//--set phase offset and rotate to antiparallel
	UInt numChains = 4;
	bool odd = false;
	for (UInt i = 0; i < numChains; i++)
	{
		if (antiParallel){
			if (odd){
				_prot->rotate(i, Z_axis, _antiphase);
				_prot->rotate(i, Y_axis, 180); odd = false;
			}
			else{
				_prot->rotate(i, Z_axis, _phase); odd = true;
			}
		}
		else{_prot->rotate(i, Z_axis, _phase);}
	}
	
	//--translate off z to radius, rotate each chain into position and set helix offset
	odd = false;
	double radial = 0.0, radialInterval = 360/numChains;
	_prot->translate(0.0, _radius, 0.0);
	for (UInt i = 0; i < numChains; i++)
	{
		_prot->rotate(i, Z_axis, radial);
		if (_offset != 0.0){
			if (odd){ 
				_prot->translate(i, 0.0, 0.0, _offset); odd = false;
			}
			else{odd = true;}
		}
		radial += radialInterval;
	}
	
	//--supercoil oligamer
	if (_coil != 0.0){_prot->coilcoil(_coil);}
	return;
}

// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
