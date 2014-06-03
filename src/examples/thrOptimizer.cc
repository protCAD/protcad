#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

int main (int argc, char* argv[])
{
    residue::setCutoffDistance(6.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

    string inputFile = argv[1];
    ifstream inFile;
    inFile.open(inputFile.c_str());
    if (!inFile)
    {
        cout << "Unable to find or open file" << endl;
        exit(1);
    }

	ran ranNumber;
	UInt randomSeed;

	string currentLine;
	vector < string> parsedStrings(0);

    if (argc == 3)
    {
        getline(inFile, currentLine, '\n'); // read in but skip
        sscanf(argv[2], "%u", &randomSeed); // use second command line arg instead
        ranNumber.setSeed(randomSeed);
    }
    else
    {
        getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
        sscanf(parsedStrings[0].c_str(), "%u", &randomSeed);
        ranNumber.setSeed(randomSeed);
    }

    UInt iterations, history;
    getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
    sscanf(parsedStrings[0].c_str(), "%u", &iterations); sscanf(parsedStrings[1].c_str(), "%u", &history);

    double startBeta, endBeta;
    getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
    sscanf(parsedStrings[0].c_str(), "%lf", &startBeta); sscanf(parsedStrings[1].c_str(), "%lf", &endBeta);

	UInt thrPos;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &thrpos);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);


    PDBInterface* thePDB1 = new PDBInterface(parsedStrings[0]);
    ensemble* theEnsemble1 = thePDB1->getEnsemblePointer();
    molecule* theMol1 = theEnsemble1->getMoleculePointer(0);
    protein* lthr = static_cast<protein*>(theMol1);

    PDBInterface* thePDB2 = new PDBInterface(parsedStrings[1]);
    ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
    molecule* theMol2 = theEnsemble2->getMoleculePointer(0);
    protein* dthr = static_cast<protein*>(theMol2);

    PDBInterface* thePDB3 = new PDBInterface(parsedStrings[2]);
    ensemble* theEnsemble3 = thePDB3->getEnsemblePointer();
    molecule* theMol3 = theEnsemble3->getMoleculePointer(0);
    protein* lallothr = static_cast<protein*>(theMol3);

	PDBInterface* thePDB4 = new PDBInterface(parsedStrings[3]);
    ensemble* theEnsemble4 = thePDB4->getEnsemblePointer();
    molecule* theMol4 = theEnsemble34->getMoleculePointer(0);
    protein* dallothr = static_cast<protein*>(theMol4);




		prot->setChi(0,thrpos,0,0,chi);
		double energy = prot->getPositionEnergy(0, thrpos);
		cout << "OG " << chi << " " << energy << endl;
	
	return 0;
}

