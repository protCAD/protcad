#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>


void buildDimer(protein* _prot, double _r, double _phase, double _phase2, double _z, double _tilt);
void undoDimer(protein* _prot, double _r, double _phase, double _phase2, double _z, double _tilt);

int main (int argc, char* argv[])
{


    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);
    	protein* bundle = static_cast<protein*>(pMol);

	string basename = argv[2];
	string logfile = basename + ".log";
	ofstream fout;
	fout.open(logfile.c_str());
	residue::setCutoffDistance(6.0);
	pmf::setScaleFactor(0.0);
    	rotamer::setScaleFactor(0.0);
    	microEnvironment::setScaleFactor(0.0);
    	amberVDW::setScaleFactor(1.0);
    	amberVDW::setRadiusScaleFactor(1.0);
    	amberVDW::setLinearRepulsionDampeningOff();
    	amberElec::setScaleFactor(0.0);
    	solvation::setItsScaleFactor(0.0);
	bundle->silenceMessages();
	for (double tilt = -180.0; tilt <= 180; tilt = tilt + 2.5)
	{
		double lowEnergy = 1E20;
		for (double radius = 7.5; radius <= 9.0; radius = radius + 0.5)
		{
			for (double phase = -60.0; phase <= 60.0; phase = phase + 5.0)
			{
				for (double phase2 = -60.0; phase2 <= 60.0; phase2 = phase2 + 5.0)
				{
					for (double z = -2.0; z <= 2.0; z = z + 0.2)
					{
						buildDimer(bundle, radius, phase, phase2, z, tilt);
						double energy = bundle->intraEnergy(0,1);
						if (energy < lowEnergy)
						{
							lowEnergy = energy;
							char tilts[20];
							sprintf(tilts, "%.1lf", tilt);
							string outfile = basename + tilts + ".pdb";
							pdbWriter(bundle, outfile);
						}
						undoDimer(bundle, radius, phase, phase2, z, tilt);
					}
				}
			}
		}
		fout << tilt << " " << lowEnergy << endl;
		cout << tilt << " " << lowEnergy << endl;
	}	
	

	return 0; 

}


void buildDimer(protein* _prot, double _r, double _phase, double _phase2, double _z, double _tilt)
{
	_prot->rotate(1, Z_axis, _phase);
	_prot->rotate(0, Z_axis, _phase2);
	_prot->translate(1, 0.0, _r, _z);
	_prot->rotate(1, Y_axis, _tilt);
	return;
}

void undoDimer(protein* _prot, double _r, double _phase, double _phase2, double _z, double _tilt)
{
	_prot->rotate(1, Y_axis, -1.0*_tilt);
	_prot->translate(1, 0.0, -1.0*_r, -1.0*_z);
	_prot->rotate(1, Z_axis, -1.0*_phase);
	_prot->rotate(0, Z_axis, -1.0*_phase2);
	return;
}
