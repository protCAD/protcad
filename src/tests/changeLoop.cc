#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>
#include "ran.h"

double getRMSD (protein* _prot, vector <UIntVec> _posArray1, vector <UIntVec> _posArray2);

int main (int argc, char* argv[])
{

//  STEP 1:  generate bundle
    	if (argc != 4)
    	{
    	exit(1);
    	}
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);
    	protein* bundle = static_cast<protein*>(pMol);
   
	string dataFileName = argv[2];
	ifstream dataFile;
	dataFile.open(dataFileName.c_str());
	if (!dataFile)
	{
		cout << "Unable to find data input file " << dataFileName << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings(0);

	// random number seed
	string seedStr = argv[3];
	int seed; 
	sscanf(seedStr.c_str(), "%d", &seed);

	ran ranNumber;
	ranNumber.setSeed(seed);

	// get number of modifiable positoins
	getline (dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	UInt numPos;  sscanf(parsedStrings[0].c_str(), "%u", &numPos);
	vector < UInt > modChain(0);
	vector < UInt > modRes(0);
	vector < double > phiMin(0);
	vector < double > phiMax(0);
	vector < double > psiMin(0);
	vector < double > psiMax(0);
	UInt modChainTemp, modResTemp;
	double phiMinTemp, phiMaxTemp, psiMinTemp, psiMaxTemp;
	for (UInt i = 0; i < numPos; i ++)
	{
		getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[0].c_str(), "%u", &modChainTemp);
		sscanf(parsedStrings[1].c_str(), "%u", &modResTemp);
		sscanf(parsedStrings[2].c_str(), "%lf", &phiMinTemp);
		sscanf(parsedStrings[3].c_str(), "%lf", &phiMaxTemp);
		sscanf(parsedStrings[4].c_str(), "%lf", &psiMinTemp);
		sscanf(parsedStrings[5].c_str(), "%lf", &psiMaxTemp);

		modChain.push_back(modChainTemp);	
		modRes.push_back(modResTemp);
		phiMin.push_back(phiMinTemp);
		phiMax.push_back(phiMaxTemp);
		psiMin.push_back(psiMinTemp);
		psiMax.push_back(psiMaxTemp);
	}

	// get number of overlap positions
	getline (dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	UInt numComp; sscanf(parsedStrings[0].c_str(), "%u", &numComp);
	
	vector < vector < UInt > > posArray1, posArray2;

	UInt c1, r1, c2, r2;

	for (UInt i = 0; i < numComp; i ++)
	{
		getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[0].c_str(), "%u", &c1);
		sscanf(parsedStrings[1].c_str(), "%u", &r1);
		sscanf(parsedStrings[2].c_str(), "%u", &c2);
		sscanf(parsedStrings[3].c_str(), "%u", &r2);
		
    		vector <UInt> tmpPos;
		tmpPos.push_back(c1); tmpPos.push_back(r1);
		posArray1.push_back(tmpPos);
		tmpPos.resize(0);
		tmpPos.push_back(c2); tmpPos.push_back(r2);
		posArray2.push_back(tmpPos);
	}

	getline(dataFile,currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	UInt iterations;  
	sscanf(parsedStrings[0].c_str(), "%u", &iterations);

	dataFile.close();
	
	for (UInt i = 0; i < numPos; i ++)
	{
		bundle->setPhi(modChain[i], modRes[i], phiMin[i] + ranNumber.getNext() * (phiMax[i] - phiMin[i]));
		bundle->setPsi(modChain[i], modRes[i], psiMin[i] + ranNumber.getNext() * (psiMax[i] - psiMin[i]));
	}

	double bestRMSD = getRMSD(bundle, posArray1, posArray2);
	double currentRMSD = bestRMSD; 
	for (UInt i = 0; i < iterations; i ++)
	{
		UInt posToModify = UInt(ranNumber.getNext()*(numPos+1));
		if (posToModify >= numPos) posToModify --;

		double oldPhi = bundle->getPhi(modChain[posToModify],modRes[posToModify]);
		double oldPsi = bundle->getPsi(modChain[posToModify],modRes[posToModify]);

		double newPhi = phiMin[posToModify] + ranNumber.getNext() * (phiMax[posToModify] - phiMin[posToModify]);
		double newPsi = psiMin[posToModify] + ranNumber.getNext() * (psiMax[posToModify] - psiMin[posToModify]);

		bundle->setPhi(modChain[posToModify],modRes[posToModify], newPhi);
		bundle->setPsi(modChain[posToModify],modRes[posToModify], newPsi);

		double RMSD = getRMSD(bundle, posArray1, posArray2);

		if (RMSD <= currentRMSD)
		{
			currentRMSD = RMSD;
			pdbWriter(bundle, "current.pdb");
			if (RMSD < bestRMSD)
			{
				bestRMSD = RMSD;
				string fileName = "best" + seedStr + ".pdb";
				pdbWriter(bundle, fileName);
				cout << "iteration " << i << "  best = " << bestRMSD << endl;
			}
		}
		else
		{
			double beta = 1E-3 + (float)i*1e-1 /(float)iterations;
			double prob = pow(2.718, (bestRMSD - RMSD) * beta);
			if (prob > ranNumber.getNext())
			{
				currentRMSD = RMSD;
				pdbWriter(bundle, "current.pdb");
			}
			else
			{
				bundle->setPhi(modChain[posToModify],modRes[posToModify], oldPhi);
				bundle->setPsi(modChain[posToModify],modRes[posToModify], oldPsi);
			}
		}
	}

	return 0;
}
			
double getRMSD(protein* _prot, vector < UIntVec > _posArray1, vector < UIntVec > _posArray2)
{
	double deviation = 0.0;
	for (UInt i = 0; i < _posArray1.size(); i ++)
	{
		dblVec CA1 = _prot->getCoords(_posArray1[i][0], _posArray1[i][1], "CA");
		dblVec CA2 = _prot->getCoords(_posArray2[i][0], _posArray2[i][1], "CA");

		deviation += pow(2.0,CMath::distance(CA1, CA2));
	}

	deviation = deviation / float(_posArray1.size());

	return sqrt(deviation);
}
