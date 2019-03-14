#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>

UInt countHbonds(protein* _prot);
double getEnergy(protein* _prot);
UInt countCAOclashes(protein* _prot);
string  intToString(int _num);
void saveState(protein* _prot, string _prefix, UInt _fileNum, UInt _step, double _beta, double _score);
bool isFixed(UIntVec _fixedPositions, UInt _res);
bool flipChirality(protein* _prot, UInt _res);
void fixCoords(protein* _prot);
double getProb(double _oldScore, double _score, double _beta);

int main (int argc, char* argv[])
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};


	//read command line
	if (argc < 2)
	{
		cout << "primalSoup inputfile" << endl;
		exit(1);
	}

	//read inputfile

	string inputFile = argv[1];
	ifstream inFile;
	inFile.open(inputFile.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}
	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	string inputFileName = parsedStrings[0];

	// read in prot structure
	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);


	residue::setCutoffDistance(4.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.9);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	ran ranNumber;
	UInt randomSeed;
	
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

	double probLAla; // enantiomeric excess of L-Ala
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &probLAla);

	//fixCoords(prot);
	pdbWriter(prot,"fixed.pdb");
	UInt size = prot->getNumResidues(0);
	for (UInt i = 0; i < size; i ++)
	{
		prot->activateForRepacking(0,i);
		prot->setPhi(0,i,180);
		prot->setPsi(0,i,180);
		if (ranNumber.getNext() > 0.5)  // randomize initial sequence
			prot->mutate(0, i, 20);
		else
			prot->mutate(0, i, 0);
	}

	UIntVec fixedPositions;
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings=Parse::parse(currentLine);
		UInt tempRes, tempID;
		sscanf(parsedStrings[0].c_str(), "%u", &tempRes);
		sscanf(parsedStrings[1].c_str(), "%u", &tempID);
		fixedPositions.push_back(tempRes);
		prot->mutate(0, tempRes, tempID);
	}


	//
	//	
	for (UInt i = 0; i < size; i++)
	{
		prot->setPhi(0,i,-65.0);
		prot->setPsi(0,i,-40.0);
		prot->mutate(0,i,0);
		cout << "RIGHT HELIX ENERGY" << getEnergy(prot) <<endl;
		pdbWriter(prot, "rightHelix.pdb");
	}
	for (UInt i = 0; i < size; i++)
	{
		prot->setPhi(0,i,65.0);
		prot->setPsi(0,i,40.0);
		prot->mutate(0,i,20);
		cout << "LEFT HELIX ENERGY" << getEnergy(prot) << endl;
		pdbWriter(prot, "leftHelix.pdb");
	}
	return 0;

}

double getProb(double _oldScore, double _score, double _beta)
{
	double dScore = _oldScore - _score;
	double prob = pow(2.718,  dScore * _beta);
	return prob;
}

double getEnergy(protein* _prot)
{

	UInt nhbs = countHbonds(_prot);
	UInt nclashes = countCAOclashes(_prot);
	double score = _prot->intraEnergy() + 5*nclashes - 5*nhbs;
	return score;
}

UInt countCAOclashes(protein* _prot)
{
	UInt size = _prot->getNumResidues(0);
	UInt clashes = 0;
	vector < dblVec > CBcoords;
	vector < dblVec > Ocoords;

	for (UInt i = 0; i < size; i ++) // do not count clashes on last residue as we cannot change its PSI anyways
	{
		dblVec tempcoord;
		UInt type = _prot->getTypeFromResNum(0,i);
		if (type == 0 || type == 20)
		{
			tempcoord =_prot->getCoords(0, i, 4);
			CBcoords.push_back(tempcoord);
		}
		tempcoord =_prot->getCoords(0, i, 3);
		Ocoords.push_back(tempcoord);
	}

	for (UInt i = 0; i < CBcoords.size(); i ++)
	{
		for (UInt j = 0; j < Ocoords.size(); j ++)
		{
			dblVec temp = CBcoords[i] - Ocoords[j];
			double distance = sqrt(CMath::dotProduct(temp,temp));
			int resSpace = (i - j);
			if (distance < 3.0 && abs(resSpace) < 2 )
			{
				clashes ++;
			}
		}
	}
	return clashes;
}



UInt countHbonds(protein* _prot)
{
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
				if (distance < 3.2 && angle > 90.0)
				{
					numHbonds ++;
				}
			}
		}
	}
	return numHbonds;
}

string intToString(int _num)
{
	ostringstream myStream;
	myStream << _num << flush;
	return (myStream.str());
}

void saveState(protein* _prot, string _prefix, UInt _fileNum, UInt _step, double _beta, double _score)
{

	string pad = "";
	if     (_fileNum < 10) pad = "000000";
	else if(_fileNum < 100) pad = "00000";
	else if(_fileNum < 1000) pad = "0000";
	else if(_fileNum < 10000) pad = "000";
	else if(_fileNum < 100000) pad = "00";
	else if(_fileNum < 1000000) pad = "0";


	string filename = _prefix + pad + intToString(_fileNum) + ".pdb";
	pdbWriter(_prot, filename);
	string logfile = _prefix + ".log";
	ofstream fout (logfile.c_str() ,ios::app);
	fout << "step " << _step << " beta " << _beta << " score " << _score << " nhbs " <<
		countHbonds(_prot)<< " clashes " <<  countCAOclashes(_prot) << " sequence ";
	UInt size = _prot->getNumResidues(0);
	for (UInt n = 0; n < size; n ++)
	{
		if (_prot->getTypeFromResNum(0,n) == 0)
		{
			fout << "L ";
		}
		else if (_prot->getTypeFromResNum(0, n) == 20)
		{
			fout << "D ";
		}
		else
		{
			fout << _prot->getTypeStringFromResNum(0,n) << " ";
		}
	}
	for (UInt n = 0; n < size; n ++)
	{
		fout << _prot->getPhi(0,n) << " ";
		fout << _prot->getPsi(0,n) << " ; ";
	}
	fout << endl;

	fout.close();
	return;
}

bool flipChirality(protein* _prot, UInt _res)
{
	if (_prot->getTypeFromResNum(0,_res) == 0)
	{
		_prot->mutate(0, _res, 20);
		return true;
	}
	else if (_prot->getTypeFromResNum(0, _res) == 20)
	{
		_prot->mutate(0, _res, 0);
		return true;
	}
	return false;
}

bool isFixed(UIntVec _fixedPositions, UInt _res)
{
	bool flag = false;
	for (UInt i = 0; i < _fixedPositions.size(); i ++)
	{
		if (_res == _fixedPositions[i])
		{
			flag = true;
		}
	}
	
	return flag;
}

void fixCoords(protein* _prot)
{
	UInt numChains = _prot->getNumChains();

	// N1 CA1 C1 O1  N2 CA2 C2 O2 N3
	double coords[] = {0.0,0.0,0.0,     1.36,0.89,0.0,    2.45,0.0,0.0,     2.45,-1.23,0.0, 
			   3.77,0.89,0.0,    4.88,0.0,0.0,    6.22,0.89,0.0,    6.22,2.12,0.0, 
			   7.29,0.0,0.0 };
       	
	for (UInt i = 0; i < numChains; i ++)
	{
		for (UInt j = 0; j < _prot->getNumResidues(i); j ++)
		{
			_prot->mutate(i, j, 7);
		}
		
		for (UInt j = 0; j < _prot->getNumResidues(i); j ++)
		{
			
			dblVec NCoords(3);
			dblVec CACoords(3);
			dblVec CCoords(3);
			dblVec OCoords(3);
			if (j % 2 == 0)
			{
				NCoords[0] = coords[0]+7.29*(float)j/2.0;
				NCoords[1] = coords[1];
				NCoords[2] = coords[2];

				CACoords[0] = coords[3]+7.29*(float)j/2.0;
				CACoords[1] = coords[4];
				CACoords[2] = coords[5];
				
				CCoords[0] = coords[6]+7.29*(float)j/2.0;
				CCoords[1] = coords[7];
				CCoords[2] = coords[8];

				OCoords[0] = coords[9]+7.29*(float)j/2.0;
				OCoords[1] = coords[10];
				OCoords[2] = coords[11];
			}
			else if (j % 2 == 1)
			{
                               	NCoords[0] = coords[12]+7.29*((float)j-1.0)/2.0;
			       	NCoords[1] = coords[13];
			       	NCoords[2] = coords[14];

                               	CACoords[0] = coords[15]+7.29*((float)j-1.0)/2.0;
                               	CACoords[1] = coords[16];
                      		CACoords[2] = coords[17];

				CCoords[0] = coords[18]+7.29*((float)j-1.0)/2.0;
				CCoords[1] = coords[19];
				CCoords[2] = coords[20];

				OCoords[0] = coords[21]+7.29*((float)j-1.0)/2.0;
				OCoords[1] = coords[22];
				OCoords[2] = coords[23];
			}
			_prot->setCoords(i,j,0,NCoords);
			_prot->setCoords(i,j,1,CACoords);
			_prot->setCoords(i,j,2,CCoords);
			_prot->setCoords(i,j,3,OCoords);
		}
	}
	return;
}
