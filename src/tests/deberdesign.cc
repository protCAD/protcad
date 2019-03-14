#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>

void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
vector <UIntVec> generateSequences (vector <UIntVec> _allowedResidues);


int main (int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

	string inputFileName = argv[1];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to open " << inputFileName << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);

	getline (inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	string pdbIn = parsedStrings[0];
	string pdbOut = parsedStrings[1];
	
	PDBInterface* thePDB = new PDBInterface(pdbIn);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);

	string parlogfile = pdbOut + "_par.log";
	ofstream fout;
	fout.open(parlogfile.c_str());
	
	prot->silenceMessages();
	prot->symmetryLinkChainAtoB(1,0);
	prot->activateAllForRepacking(0);
	prot->setCanonicalHelixRotamersOnly(0);

	residue::setCutoffDistance(6.0);
	pmf::setScaleFactor(0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOn();
	amberElec::setScaleFactor(0.0);

	double radstart, radend, radstep, phasestart, phaseend, phasestep, coil;
	getline (inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &radstart);
	sscanf(parsedStrings[1].c_str(), "%lf", &radend);
	sscanf(parsedStrings[2].c_str(), "%lf", &radstep);

	getline (inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &phasestart);
	sscanf(parsedStrings[1].c_str(), "%lf", &phaseend);
	sscanf(parsedStrings[2].c_str(), "%lf", &phasestep);
	
	getline (inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &coil);
	
	UInt threadstart, threadend;
	getline (inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &threadstart);
	sscanf(parsedStrings[1].c_str(), "%u", &threadend);

	vector < UIntVec > allowedResidues(0);
	while (getline(inFile, currentLine, '\n'))
	{
		UIntVec tempAllowed(0);
		parsedStrings = Parse::parse(currentLine);
		for (UInt j = 1; j < parsedStrings.size(); j++)
		{
			UInt resId = 1000;
			if (parsedStrings[j] == "A") resId = A;
			if (parsedStrings[j] == "R") resId = R;
			if (parsedStrings[j] == "N") resId = N;
			if (parsedStrings[j] == "D") resId = D;
			if (parsedStrings[j] == "C") resId = C;
			if (parsedStrings[j] == "Q") resId = Q;
			if (parsedStrings[j] == "E") resId = E;
			if (parsedStrings[j] == "G") resId = G;
			if (parsedStrings[j] == "H") resId = H;
			if (parsedStrings[j] == "I") resId = I;
			if (parsedStrings[j] == "L") resId = L;
			if (parsedStrings[j] == "K") resId = K;
			if (parsedStrings[j] == "M") resId = M;
			if (parsedStrings[j] == "F") resId = F;
			if (parsedStrings[j] == "P") resId = P;
			if (parsedStrings[j] == "S") resId = S;
			if (parsedStrings[j] == "T") resId = T;
			if (parsedStrings[j] == "W") resId = W;
			if (parsedStrings[j] == "Y") resId = Y;
			if (parsedStrings[j] == "V") resId = V;
			if (resId == 1000)
			{
				cout << "ERROR - bad amino acid designation ... quitting." << endl;
				exit(1);
			}
			tempAllowed.push_back(resId);
		}
		allowedResidues.push_back(tempAllowed);
	}
	
	UInt numSequences = 1;
	for (UInt i = 0; i < allowedResidues.size(); i ++)
	{
		numSequences *= allowedResidues[i].size();
	}
	vector <UIntVec> allowedSequences = generateSequences(allowedResidues);
	for (UInt seqindex = 0; seqindex < numSequences; seqindex ++)
	{
		UIntVec sequence = allowedSequences[seqindex];
		string bestoutname = "best";
		double bestoutenergy = 1e10;
		double bestoutrad = 0.0;
		double bestoutphase = 0.0;
		UInt bestoutstart = 0;
		for (UInt tmpi = 0; tmpi < sequence.size(); tmpi ++)
		{
			prot->mutate(0,0,sequence[tmpi]);
			fout << prot->getTypeStringFromResNum(0,0) << " ";
			bestoutname = bestoutname + "_" + prot->getTypeStringFromResNum(0,0);
		}
		bestoutname = bestoutname + ".pdb";
		prot->mutate(0, 0, A);
		for (UInt startpos = threadstart; startpos <= threadend; startpos ++)
		{
			for (UInt i = 0; i < prot->getNumResidues(0); i ++)
			{
				prot->mutate(0, i, A);
			}
			for (UInt i = 0; i < sequence.size(); i ++)
			{
				prot->mutate(0, startpos + i, sequence[i]);
			}
			//buildParallelDimer(prot, 10.0, 0.0, coil);
			//prot->optimizeRotamers();
			//undoParallelDimer(prot, 10.0, 0.0, coil);
			double bestrad = radstart;
			double bestphase = phasestart;
			double bestEnergy = 1E10;
			for (double rad = radstart; rad <= radend; rad = rad + radstep)
			{
				for (double phase = phasestart; phase <= phaseend; phase = phase + phasestep)
				{
					buildParallelDimer(prot, rad, phase, coil);
					prot->optimizeRotamers();
					double energy = prot->intraEnergy(0,1); 
					if (energy < bestEnergy)
					{
						bestEnergy = energy;
						bestrad = rad;
						bestphase = phase;
						if (bestEnergy < bestoutenergy)
						{
							bestoutenergy = bestEnergy;
							bestoutrad = bestrad;
							bestoutphase = bestphase;
							bestoutstart = startpos;
							pdbWriter(prot, bestoutname);
						}	
					}
					undoParallelDimer(prot, rad, phase, coil);
				}
			}
			fout << " |s| " << startpos << "  r " << bestrad << " p " << bestphase << " e " << bestEnergy;
		}
		fout << endl;
	}
	return 0;
}



vector <UIntVec> generateSequences (vector <UIntVec> _allowedResidues)
{
	vector < UIntVec > sequences(0);
	for (UInt a = 0; a < _allowedResidues[0].size(); a ++)
	{
		for (UInt b = 0; b < _allowedResidues[1].size(); b ++)
		{
			for (UInt c = 0; c < _allowedResidues[2].size(); c ++)
			{
				for (UInt d = 0; d < _allowedResidues[3].size(); d ++)
				{
					for (UInt e = 0; e < _allowedResidues[4].size(); e ++)
					{
						for (UInt f = 0; f < _allowedResidues[5].size(); f ++)
						{
							for (UInt g = 0; g < _allowedResidues[6].size(); g ++)
							{
								UIntVec tmpSequence(0);
								tmpSequence.push_back(_allowedResidues[0][a]);
								tmpSequence.push_back(_allowedResidues[1][b]);
								tmpSequence.push_back(_allowedResidues[2][c]);
								tmpSequence.push_back(_allowedResidues[3][d]);
								tmpSequence.push_back(_allowedResidues[4][e]);
								tmpSequence.push_back(_allowedResidues[5][f]);
								tmpSequence.push_back(_allowedResidues[6][g]);
								sequences.push_back(tmpSequence);
							}
						}
					}
				}
			}
		}
	}

	return sequences;
}

void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->rotate(Z_axis, _phase);
	_prot->translate(0.0, _radius, 0.0);
	_prot->rotate(1, Z_axis, 180);
	if (_coil != 0.0) _prot->coilcoil(_coil);
	return;
}

void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	if (_coil != 0.0)_prot->coilcoil(-1.0 * _coil);
	_prot->rotate(1, Z_axis, 180);
	_prot->translate(0.0, -1.0*_radius, 0.0);
	_prot->rotate(Z_axis, -1.0 * _phase);
	return;
}

