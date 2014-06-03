#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>
#include<cstdlib>
#include<ctime>



void printTable(protein* _prot);
vector < double > getPhiPsi (vector < vector < double > > _table, ran &_ranNum);
UInt countHbonds(protein* _prot);
double getEndToEndDistance(protein* _prot);
double getMaxDistance(protein* _prot);

double getRadiusOfGyration(protein* _prot);
void mapCoords(protein* _prot, vector <dblVec> _itsCoords);
vector <dblVec> getCoords(protein* _prot);

// usage:  evenMC file.pdb ranseed steps 0 20 0 20 20 ...

int main (int argc, char* argv[])
{

        enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X};

	string inputFileName = argv[1];

       	// read in prot structure
        PDBInterface* thePDB = new PDBInterface(inputFileName);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* theMol = theEnsemble->getMoleculePointer(0);
        protein* prot = static_cast<protein*>(theMol);

	// set energy function parameters - key for this program is radius scale, set to 0.9
        residue::setCutoffDistance(8.0);
        pmf::setScaleFactor(0.0);
        rotamer::setScaleFactor(0.0);
        microEnvironment::setScaleFactor(0.0);
        amberVDW::setScaleFactor(1.0);
        amberVDW::setRadiusScaleFactor(1.0);
        amberVDW::setLinearRepulsionDampeningOff();
        amberElec::setScaleFactor(0.0);
        solvation::setItsScaleFactor(0.0);

	prot->silenceMessages();
	prot->activateAllForRepacking(0);

	string stub = argv[2];
	string outFile = stub + ".out";
	ofstream fout;
	fout.open(outFile.c_str());

	string outPDBFile = stub + ".pdb";

	vector < vector < double > > table;

	string phipsifile = argv[3];
	ifstream inFile;
	inFile.open(phipsifile.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open " << phipsifile << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		vector < double > entry;
		for (UInt i = 0; i < parsedStrings.size(); i ++)
		{
			double num;
			sscanf(parsedStrings[i].c_str(), "%lf", &num);
			entry.push_back(num);
		}
		table.push_back(entry);
	}

 

	UInt steps, seed;
	sscanf(argv[4], "%u", &seed);
	sscanf(argv[5], "%u", &steps);	

	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		UInt ID;
		prot->setPhi(0,i, 180.0);
		prot->setPsi(0,i, 180.0);
		sscanf (argv[6+i], "%u", &ID);
		prot->mutate(0,i,ID);
	}

	ran ranNumber;
	ranNumber.setSeed(seed);

	vector <dblVec> itsCoords;

	itsCoords = getCoords(prot);

	double bestE = 1E10;
	UInt goodCts = 0;
	for (UInt iter = 0; goodCts < steps; iter ++ )
	{
		if (iter > 0 && iter%10000 == 0) // prevent warping from repeated rotations
		{
			mapCoords(prot, itsCoords);
		}
		for (UInt i = 0; i < prot->getNumResidues(0); i ++)
		{
			vector < double > phipsi = getPhiPsi(table, ranNumber);

			double sign = 1.0;
			if (prot->getTypeFromResNum(0,i) == 20) sign = -1.0;
			if (prot->getTypeFromResNum(0,i) == 7 )
	 		{
				sign = ranNumber.getNext();
        			if (sign < 0.5) sign = -1.0; else sign = 1.0;
			}
		
			prot->setPhi(0,i,sign*phipsi[0]);
			prot->setPsi(0,i,sign*phipsi[1]);
		}

		double energy = prot->intraEnergy();
		if (energy < bestE)
		{
			bestE = energy;
			pdbWriter (prot, outPDBFile.c_str());
		}

		cout << "ENERGY " << energy  << endl;
		double rG = getRadiusOfGyration(prot);
		fout << "Eng " << energy << " rG " << rG  << endl ;
		goodCts ++;
	}	
				

	return 0;
}

double getRadiusOfGyration(protein* _prot)
{
	dblVec centroid = _prot->getCoords(0,0,"CA");
	for (UInt i = 1; i < _prot->getNumResidues(0); i ++)
	{
		centroid = centroid + _prot->getCoords(0,i, "CA");
	}

	centroid = centroid / (double)_prot->getNumResidues(0);

	double sum = 0;
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		dblVec diff = centroid - _prot->getCoords(0,i,"CA");
		sum = sum + CMath::dotProduct(diff,diff);
	}

	sum = sum / (double)_prot->getNumResidues(0);
	return sqrt(sum);
}

double getEndToEndDistance(protein* _prot)
{
	dblVec begin = _prot->getCoords(0,0, "CA");
	dblVec end = _prot->getCoords(0, _prot->getNumResidues(0) - 1, "CA");

	dblVec diff = begin - end;	

	return sqrt(CMath::dotProduct(diff,diff));
}

double getMaxDistance(protein* _prot)
{
	vector <dblVec> CAcoords;
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		CAcoords.push_back(_prot->getCoords(0,i,"CA"));
	}

	double maxDist = 0.0;

	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		for (UInt j = i + 1; j < _prot->getNumResidues(0); j ++)
		{
			dblVec diff = CAcoords[i] - CAcoords[j];
			double dist = sqrt(CMath::dotProduct(diff, diff));
			if (dist > maxDist) maxDist = dist;
		}
	}

	return maxDist;
}


vector < double > getPhiPsi (vector < vector < double > > _table, ran &_ranNum)
{
	double percentage = 0;
	int index = -1;
	UInt i = 0;

	double num = _ranNum.getNext();

	while (index < 0 && i < _table.size())
	{
		percentage += _table[i][2];
		//cout << percentage << " ";
		if ( (100.0 * num) < percentage) index = i;
		i ++;
	}



	if (index == -1) index = _table.size() - 1;

	//cout << num << " -> " << index << " " << _table[index][0] << " " << _table[index][1] <<  endl;

	double sign = _ranNum.getNext();
	if (sign < 0.5) sign = -1.0; else sign = 1.0;

	double phiHarmonic = _table[index][0] + sign * pow(_ranNum.getNext(), 2.0) * 30.0;

	sign = _ranNum.getNext();
	if (sign < 0.5) sign = -1.0; else sign = 1.0;

	double psiHarmonic = _table[index][1] + sign * pow(_ranNum.getNext(), 2.0) * 30.0; 

	vector < double > phipsi;

	//cout << "phipsi " << phiHarmonic - _table[index][0]  << " " << psiHarmonic - _table[index][1] << endl;

	phipsi.push_back(phiHarmonic);
	phipsi.push_back(psiHarmonic);

	return phipsi;

}		


void printTable(protein* _prot)
{
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		cout << i << " PHI " << _prot->getPhi(0,i) << " PSI " << _prot->getPsi(0,i) << endl;
	}
	return;
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
				if (distance < 3.2 && angle > 120.0)
				{
					numHbonds ++;
				}
			}
		}
	}
	return numHbonds;
}

void mapCoords(protein* _prot, vector <dblVec> _itsCoords)
{
	UInt n = 0;

	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		for (UInt j = 0; j < _prot->getNumAtoms(0,i); j ++)
		{
			_prot->setCoords(0,i,j,_itsCoords[n]);
			n ++;
		}
	}

	return;
}			
	

vector <dblVec> getCoords(protein* _prot)
{
	vector <dblVec> itsCoords;

	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		for (UInt j = 0; j < _prot->getNumAtoms(0,i); j ++)
		{
			itsCoords.push_back(_prot->getCoords(0,i,j));
		}
	}
	
	return itsCoords;
}

	





