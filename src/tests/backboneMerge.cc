#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>

#define ITERMAX 10000
#define ANGLEVAR 15.0

double calculateEnergy(protein* _prot);
vector < vector < double > > getCurrentPhiPsi(UInt _start, UInt _stop, protein* _prot);

int main (int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};


	// READ IN PDB FILE

	string pdbIn = argv[1];
	
	PDBInterface* thePDB = new PDBInterface(pdbIn);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);

	residue::setCutoffDistance(6.0);
	pmf::setScaleFactor(0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.9);
	amberVDW::setLinearRepulsionDampeningOn();
	amberElec::setScaleFactor(0.0);

	// READ IN PHI PSI TABLES

	string inputFileName = argv[2];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to open " << inputFileName << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	vector <vector <vector <double> > > project;
	project.resize(0);

	
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		if (project.size() == 0)
		{
			project.resize(parsedStrings.size() / 2);
		}

		for (UInt i = 0; i < parsedStrings.size(); i = i + 2)
		{
			double phi, psi;
			sscanf(parsedStrings[i].c_str(), "%lf", &phi);
			sscanf(parsedStrings[i+1].c_str(), "%lf", &psi);

			vector < double > pair;
			pair.push_back(phi);
			pair.push_back(psi);

			project[i/2].push_back(pair);
		}
	}

	// SET RANDOMSEED FROM INPUT

	string seedInp = argv[3];
	UInt ranSeed;
	
	sscanf(seedInp.c_str(), "%u", &ranSeed);
	ran ranNumber;
	ranNumber.setSeed(ranSeed);

	// START MAIN SIMULATION


	double oldEnergy = calculateEnergy(prot);
	prot->saveCurrentState();  // save current rotamer configuration

	for (UInt iter = 0; iter < ITERMAX; iter ++)
	{
		UInt start = (int)(ranNumber.getNext()*(prot->getNumResidues(0)-2));
		UInt del  = (int)(ranNumber.getNext()*(prot->getNumResidues(0) - start));

		UInt stop = start + del;

		UInt model = (int)(ranNumber.getNext()*project.size());


		// save current phipsi values
		vector < vector <double> > currentPhiPsi = getCurrentPhiPsi(start, stop, prot);
		
		// choose new phipsi values
		for (UInt i = start; i <= stop; i ++)
		{
			double phi = (project[model][i][0] - ANGLEVAR) + ranNumber.getNext()*2.0*ANGLEVAR;
			double psi = (project[model][i][1] - ANGLEVAR) + ranNumber.getNext()*2.0*ANGLEVAR; 
	
			prot->setPhi(0,i,phi);
			prot->setPsi(0,i,psi);	
		}

		double newEnergy = calculateEnergy(prot);
		if (newEnergy <= oldEnergy)
		{
			oldEnergy = newEnergy;
			cout << "CHANGE ACCEPTED" << endl;
		}
		else
		{
			prot->undoState(); // revert to old rotamers
			prot->saveCurrentState();
			for (UInt i = start; i <= stop; i ++)
			{
				prot->setPhi(0,i,currentPhiPsi[i-start][0]);
				prot->setPsi(0,i,currentPhiPsi[i-start][1]);
			}
			cout << "CHANGE REJECTED" << endl;
		}
	}

	pdbWriter(prot, "final.pdb");


	return 0;
}

double calculateEnergy(protein* _prot)
{
	return _prot->intraEnergy();
}

vector < vector < double > > getCurrentPhiPsi(UInt _start, UInt _stop, protein* _prot)
{

	vector < vector < double > > currentPhiPsi;

	for (UInt x = _start; x <= _stop; x ++)
	{
		double phi = _prot->getPhi(0,x);
		double psi = _prot->getPsi(0,x);

		vector < double > pair;
		pair.push_back(phi);
		pair.push_back(psi);
		
		currentPhiPsi.push_back(pair);
	}

	return currentPhiPsi;
}

