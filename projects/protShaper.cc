//*******************************************************************************************************
//*******************************************************************************************************
//************************************                    ***********************************************
//************************************  protShaper 1.0  ***********************************************
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
        cout << "protShaper <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};

	//--running parameters
	residue::setTemperature(300);
	residue::setElectroSolvationScaleFactor(0.0);
	residue::setHydroSolvationScaleFactor(0.0);
	amberElec::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	srand (time(NULL));

	string outFile, inFile = argv[1];

    
    double radius;
    int count=0;
	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
	start->alignToAxis(Z_axis);
    for (int d = 0; d < 360; d++)
    {
        for (int p = 0; p < 360; p++)
        {
			for (int r = 0; r < 10; r++)
			{
				//for (int c = 170; c < 190; c++)
				//{
				
					backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
					prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[0],0,0);
					prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[1],1,0);
					count++;
					radius = 15.5 + (r*0.1);
					protein* bundle = new protein(*start);
					buildSymmetricOligamer (bundle, true, radius, 320, p, d, 2.5);
					cout << count << " " << radius << " " << d << " " << p << " " << bundle->getNumHardClashes() << endl;
					stringstream convert;
					string countstr;
					convert << count, countstr = convert.str();
					outFile = countstr + "hethelix.pdb";
					pdbWriter(bundle, outFile);
					delete bundle;
				//}
			}
            p=p+9;
        }
        d=d+9;
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

// antiparallel_beta test.pdb 5 0 -38.5 180 0 0 test2.pdb
