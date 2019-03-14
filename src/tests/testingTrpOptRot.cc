#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"
#include "annealer.h"

#define PI 3.14159

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		cout << "testingTrpOptRot pdbFileName " << endl;
		exit(1);
	}
	
	//open pdb file... (of monomeric helices interacting...)
	string pdbName = argv[1];
	ifstream pdbFile;
	pdbFile.open(pdbName.c_str());
	if (!pdbFile)
	{
		cout << "Quitting ... problem opening pdb file\n";
		exit(1);
	}
	
	cout << "reading in protein structure for pdb file: " << pdbName << endl;
	//read in prot structure
	PDBInterface* thePDB = new PDBInterface(pdbName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	cout << "force field data" << endl;
	//set forcefield data for that protein structure
	residue::setCutoffDistance(12.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.9);//soft vdw...
	amberVDW::setLinearRepulsionDampeningOff();//turn linear repulsion dampening on for now...
	amberElec::setScaleFactor(0.0);
	amberElec::setDielectricConstant(10.0);
	amberElec::distanceDependance = true;
	solvation::setItsScaleFactor(0.0);
	//double hydrogenBondScaleFactor = 0.0;
	
	//how many chains in protein...
	cout << "number of chains in protein = " << prot->getNumChains() << endl;
	
	//activate side chains 
	cout << "activate side chains" << endl;
	prot-> activateAllForRepacking(0);
	prot-> activateAllForRepacking(1);
	prot-> setCanonicalHelixRotamersOnly(0);
	prot-> setCanonicalHelixRotamersOnly(1);
	
	
	//optimize rotamers
	prot->optimizeRotamers();
	
	//get energy
	double energy = prot->intraEnergy();
	cout << "energy for dimer is " << energy << endl;
	
	//print out pdb file after rotamer optimization...
	pdbWriter (prot, "monomericHelicesAfterRotOpt.pdb");
	
	return 0;
}
