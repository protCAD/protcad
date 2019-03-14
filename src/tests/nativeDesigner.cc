#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>

UInt countHbonds(protein* _prot);
double getEnergy(protein* _prot);
UInt countCAOclashes(protein* _prot);
string  intToString(int _num);

double getProb(double _oldScore, double _score, double _beta);

int main (int argc, char* argv[])
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

	string inputFileName = argv[1];
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


	UInt numChains = prot->getNumChains();
	for (UInt j=0; j < numChains; j ++)
	{
		UInt size = prot->getNumResidues(numChains);
		for (UInt i = 0; i < size; i ++)
		{
			prot->activateForRepacking(0,i);
			prot->setPhi(0,i,180);
			prot->setPsi(0,i,180);
		}
	}

	int randomSeed;
	sscanf(argv[2], "%d", &randomSeed);
	ran ranNumber;
	ranNumber.setSeed(randomSeed);

	UInt iterations;
	sscanf(argv[3], "%u", &iterations);
	
	UInt history;
	sscanf(argv[4], "%u", &history);

	prot->silenceMessages();
	double lowestScore = 1000;
	double oldScore = lowestScore;
	double startBeta = 1e-6;
	double endBeta = 1;
	UInt lowCounter = 0;

	UInt counter = 0;
	for (UInt i = 0; i < iterations; i ++)
	{

		counter ++;
		double beta = startBeta + i * (endBeta - startBeta) / iterations;


		UInt res = UInt(ranNumber.getNext()*size);
		if (res >= size) res = size -1;

		double choice = ranNumber.getNext();
		if (choice > 0.75)
		{
			if (prot->getTypeFromResNum(0,res) == 20)
			{
				prot->mutate(0, res, 0);
			}
			else
			{
				prot->mutate(0, res, 20);
			}

			double score = getEnergy(prot);
			if (score < oldScore)
			{
				oldScore = score;
				//pdbWriter(prot, "current.pdb");
				if (score < lowestScore)
				{
					lowestScore = score;
					cout << "LOW ";
					lowCounter ++;
					string filename = "low" + intToString(lowCounter) + ".pdb";
					pdbWriter(prot, filename);
					ofstream fout("low.log",ios::app);
					fout << "step " << i << " lowcount " << lowCounter << " beta " << beta << " score " << oldScore << " low " << lowestScore << " nhbs " <<
				 		countHbonds(prot)<< " clashes " <<  countCAOclashes(prot) << " sequence ";

					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getTypeStringFromResNum(0,n) << " ";
					}
					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getPhi(0,n) << " ";
						fout << prot->getPsi(0,n) << "  ";
					}
					fout.close();
				}
			}
			else
			{
				double prob = getProb(oldScore, score, beta);
			        if ( prob > ranNumber.getNext() )
				{
					oldScore = score;
					//pdbWriter(prot, "current.pdb");
				}
				else
				{
					if (prot->getTypeFromResNum(0,res) == 20)
					{
						prot->mutate(0, res, 0);
					}
					else
					{
						prot->mutate(0, res, 20);
					}
				}
			}
		}
		else if (choice > 0.25)
		{
			double oldPhi = prot->getPhi(0,res);
			double oldPsi = prot->getPsi(0,res);
			double newPhi = 180 - 360*ranNumber.getNext();
			double newPsi = 180 - 360*ranNumber.getNext();
			prot->setPhi(0,res,newPhi);
			prot->setPsi(0,res,newPsi);

			double score = getEnergy(prot);
			if (score < oldScore)
			{
				oldScore = score;
				//pdbWriter(prot, "current.pdb");
				if (score < lowestScore)
				{
					lowestScore = score;
					cout << "LOW ";
					lowCounter ++;
					string filename = "low" + intToString(lowCounter) + ".pdb";
					pdbWriter(prot, filename);
					pdbWriter(prot, filename);
					ofstream fout("low.log",ios::app);
					fout << "step " << i << " lowcount " << lowCounter << " beta " << beta << " score " << oldScore << " low " << lowestScore << " nhbs " <<
				 		countHbonds(prot)<< " clashes " <<  countCAOclashes(prot) << " sequence ";
					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getTypeStringFromResNum(0,n) << " ";
					}
					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getPhi(0,n) << " ";
						fout << prot->getPsi(0,n) << "  ";
					}
					fout << endl;

					fout.close();
				}
			}
			else
			{
				double prob = getProb (oldScore, score, beta);
				if ( prob > ranNumber.getNext() )
				{
					oldScore = score;
					//pdbWriter(prot, "current.pdb");
				}
				else
				{
					prot->setPhi(0, res, oldPhi);
					prot->setPsi(0, res, oldPsi);
				}
			}
		}
		else
		{
                        double oldPhi = prot->getPhi(0,res);
			double oldPsi = prot->getPsi(0,res);
			double newPhi = 180 - 360*ranNumber.getNext();
			double newPsi = 180 - 360*ranNumber.getNext();
			prot->setPhi(0,res,newPhi);
			prot->setPsi(0,res,newPsi);
			if (prot->getTypeFromResNum(0,res) == 20)
			{
				prot->mutate(0, res, 0);
			}
			else
			{
				prot->mutate(0, res, 20);
			}

			double score = getEnergy(prot);
			if (score < oldScore)
			{
				oldScore = score;
				//pdbWriter(prot, "current.pdb");
				if (score < lowestScore)
				{
					lowestScore = score;
					cout << "LOW ";
					lowCounter ++;
					string filename = "low" + intToString(lowCounter) + ".pdb";
					pdbWriter(prot,  filename);
					pdbWriter(prot, filename);
					ofstream fout("low.log",ios::app);
					fout << "step " << i << " lowcount " << lowCounter << " beta " << beta << " score " << oldScore << " low " << lowestScore << " nhbs " <<
						countHbonds(prot)<< " clashes " <<  countCAOclashes(prot) << " sequence ";
					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getTypeStringFromResNum(0,n) << " ";
					}
					for (UInt n = 0; n < size; n ++)
					{
						fout << prot->getPhi(0,n) << " ";
						fout << prot->getPsi(0,n) << "  ";
					}
					fout << endl;

					fout.close();
				}
			}
			else
			{
				double prob = getProb(oldScore, score, beta);
			        if ( prob > ranNumber.getNext() )
				{
					oldScore = score;
					//pdbWriter(prot, "current.pdb");
				}
				else
				{
					if (prot->getTypeFromResNum(0,res) == 20)
					{
						prot->mutate(0, res, 0);
					}
					else
					{
						prot->mutate(0, res, 20);
					}
					prot->setPhi(0, res, oldPhi);
					prot->setPsi(0, res, oldPsi);
				}
			}
		}
		cout << "Score " << i << " " << oldScore << " temp " << beta << " lowest " << lowestScore << endl;
		if (counter == history)
		{
			counter = 0;
			ofstream fout("out.log",ios::app);
			fout << "step " << i << " beta " << beta << " score " << oldScore << " low " << lowestScore << " nhbs " <<
				 countHbonds(prot)<< " clashes " <<  countCAOclashes(prot) << " sequence ";
			for (UInt n = 0; n < size; n ++)
			{
				fout << prot->getTypeStringFromResNum(0,n) << " ";
			}
			for (UInt n = 0; n < size; n ++)
			{
				fout << prot->getPhi(0,n) << " ";
				fout << prot->getPsi(0,n) << "  ";
			}
			fout << endl;
			fout.close();
			string num = intToString(i);
			string filename = "struct" + num + ".pdb";
			pdbWriter(prot, filename);
		}
	}

	cout << "Score " << oldScore << endl;

	pdbWriter(prot,"out.pdb");
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

	for (UInt i = 0; i < size; i ++)
	{
		dblVec tempcoord =_prot->getCoords(0, i, 4);
		CBcoords.push_back(tempcoord);
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

	for (UInt i = 0; i < size; i ++)
	{
		dblVec tempcoord =_prot->getCoords(0, i, 0);
		Ncoords.push_back(tempcoord);
		tempcoord =_prot->getCoords(0, i, 3);
		Ocoords.push_back(tempcoord);
	}

	for (UInt i = 0; i < Ncoords.size(); i ++)
	{
		for (UInt j = 0; j < Ocoords.size(); j ++)
		{
			dblVec temp = Ncoords[i] - Ocoords[j];
			double distance = sqrt(CMath::dotProduct(temp,temp));
			int resSpace = (i - j);
			if (distance < 3.25 && abs(resSpace) > 2 )
			{
				numHbonds ++;
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
