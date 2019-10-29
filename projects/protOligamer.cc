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

void buildSymmetricOligamer (protein* _prot, bool antiParallel, double _radius, double _coil, double _phaseoffset1, double _phaseoffset2, double _offset);

enum structure {C,L,P,T,E,Y,A,I,D,Q,R,F,H,W,K,S};
int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
        cout << "protOligamer <inFile.pdb>" << endl;
		exit(1);
	}

	//--running parameters
	residue::setElectroSolvationScaleFactor(0.0);
	residue::setHydroSolvationScaleFactor(0.0);
	residue::setPolarizableElec(false);
	amberElec::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);
	srand (time(NULL));

	string outFile, inFile = argv[1];

	int count=0;
	double radius;
	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* start = static_cast<protein*>(pMol);
	start->alignToAxis(Z_axis);
	for (int d = 0; d < 360; d++)
	{
		for (int p = 0; p < 360; p++)
		{
			//for (int r = 5; r < 12; r++)
			//{
				count++;
				radius = 12;
				protein* bundle = new protein(*start);
				buildSymmetricOligamer (bundle, true, radius, 0.0, d, p, 0);
				cout << count << " " << d << " " << p << " " << endl;
				stringstream convert;
				string countstr;
				convert << count, countstr = convert.str();
				outFile = countstr + "_oligamer.pdb";
				pdbWriter(bundle, outFile);
				delete bundle;
			//}
			p=p+29;
		}
		d=d+29;
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
