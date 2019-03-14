#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
double calcCaDME (protein* _prot1, protein* _prot2);
void mapGCN4DimerSequence (protein* _prot);

int main (int argc, char* argv[])
{
	cout << "starting" << endl;
//  STEP 1:  generate bundle
    	if (argc != 4)
    	{
        	cout << "dimerBuilder  <infilename.pdb> <target.pdb>  <outfilename.pdb> " << endl;
    	exit(1);
    	}
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);

    	protein* bundle = static_cast<protein*>(pMol);

	string targetfile = argv[2];
	PDBInterface* thePDB2 = new PDBInterface(targetfile);
	ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
	molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
	
	protein* target = static_cast<protein*>(pMol2);

	string bestFile = argv[3];

    	residue::setCutoffDistance(6.0);
    	pmf::setScaleFactor(0.0);
    	rotamer::setScaleFactor(1.0);
    	microEnvironment::setScaleFactor(0.0);
    	amberVDW::setScaleFactor(1.0);
    	amberVDW::setRadiusScaleFactor(1.0);
    	amberVDW::setLinearRepulsionDampeningOn();
    	amberElec::setScaleFactor(0.0);
    	solvation::setItsScaleFactor(0.0);

	cout << "DME " << calcCaDME(bundle, target);

	double bestDME = 1e20;
	double coil = 150.0;
	string logfile = bestFile + ".log";
	ofstream fout;
	fout.open(logfile.c_str(), ios::app);

	bundle->silenceMessages();
	bundle->symmetryLinkChainAtoB(1,0);
	bundle->activateAllForRepacking(0);
	bundle->setCanonicalHelixRotamersOnly(0);

	fout << "mutations:  ";
	for (UInt x = 0; x < bundle->getNumResidues(0); x ++)
	{
		fout << x << " ";
		for (UInt y = x + 1; y < bundle->getNumResidues(0); y ++)
		{
			fout << x << "," << y << " ";
		}
	}
	fout << endl;
	
	for (double radius = 3.0; radius <= 7.0; radius = radius + 1.0)
	{
		for (double phase = -180.0; phase < 170.0; phase = phase + 10.0)
		{
			cout << "r " << radius << " p " << phase << " c " << coil << " ";
			fout << "r " << radius << " p " << phase << " c " << coil << " ";
			buildParallelDimer(bundle, radius, phase, coil);
			double DME = calcCaDME(bundle, target);
			cout << "DME " << DME << " ";
			fout << "DME " << DME << " ";
			mapGCN4DimerSequence(bundle);
			bundle->optimizeRotamers();
			bundle->saveCurrentState();
			if (bestDME > DME)
			{
				string bestName = bestFile + ".pdb";
				bestDME = DME;
				pdbWriter(bundle, bestName.c_str());
			}
			double wtEnergy = bundle->intraEnergy();
			cout << wtEnergy;
			fout << wtEnergy;
			for (UInt i = 0; i < bundle->getNumResidues(0); i ++)
			{
				double ebefore = bundle->getPositionEnergy(0,i);
				bundle->mutate(0, i, A);
				double eafter = bundle->getPositionEnergy(0,i);
				fout << eafter-ebefore << " ";
				bundle->undoState();
				bundle->saveCurrentState();
				for (UInt j = i + 1; j < bundle->getNumResidues(0); j ++)
				{
					double e2before = bundle->getPositionEnergy(0,i) + bundle->getPositionEnergy(0,j);
					bundle->mutate(0,i, A);
					bundle->mutate(0,j, A);
					double e2after = bundle->getPositionEnergy(0,i) + bundle->getPositionEnergy(0,j);
					fout << e2after-e2before << " ";
					bundle->undoState();
					bundle->saveCurrentState();
				}	
			}
			fout << endl;
			cout << endl;			
			undoParallelDimer(bundle, radius, phase, coil);
		}
	}

	fout.close();
	fout.clear();
	
	ofstream phenotype;
	string phile = bestFile + "_phen.dat";
	phenotype.open(phile.c_str(), ios::app);

	target->silenceMessages();
	target->symmetryLinkChainAtoB(1,0);
	target->activateAllForRepacking(0);
	target->saveCurrentState();
	
	for (UInt i = 0; i < target->getNumResidues(0); i ++)
	{
		double ebefore = target->getPositionEnergy(0,i);
		target->mutate(0, i, A);
		double eafter = target->getPositionEnergy(0,i);
		phenotype << i << " " << eafter-ebefore << " ";
		if (eafter > ebefore) phenotype << " D ";
		else phenotype << " N ";
		target->undoState();
		target->saveCurrentState();
		for (UInt j = i + 1; j < target->getNumResidues(0); j ++)
		{
			double e2before = target->getPositionEnergy(0,i) + target->getPositionEnergy(0,j);
			target->mutate(0,i, A);
			target->mutate(0,j, A);
			double e2after = target->getPositionEnergy(0,i) + target->getPositionEnergy(0,j);
			phenotype << i << "," << j << " " << e2after-e2before << " ";
			if (e2after > e2before) phenotype << " D ";
			else phenotype << " N ";
			target->undoState();
			target->saveCurrentState();
		}	
	}
	phenotype << endl;
	phenotype.close();
	phenotype.clear();
	return 0; 

}

void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->rotate(Z_axis, _phase);
	_prot->translate(0.0, _radius, 0.0);
	_prot->rotate(1, Z_axis, 180);
	_prot->coilcoil(_coil);
	return;
}

void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->coilcoil(-1.0 * _coil);
	_prot->rotate(1, Z_axis, 180);
	_prot->translate(0.0, -1.0*_radius, 0.0);
	_prot->rotate(Z_axis, -1.0 * _phase);
	return;
}


double calcCaDME (protein* _prot1, protein* _prot2)
{

	vector < vector <dblVec> > prot1trace(0);
	vector < vector <dblVec> > prot2trace(0);
	for (UInt i = 0; i < _prot1->getNumChains(); i ++)
	{
		vector < dblVec > tmp1(0);
		vector < dblVec > tmp2(0);
		for (UInt j = 0; j < _prot1->getNumResidues(i); j ++)
		{
			dblVec itsCoords = _prot1->getCoords(i, j, "CA");
			tmp1.push_back(itsCoords);
			itsCoords = _prot2->getCoords(i, j, "CA");
			tmp2.push_back(itsCoords);
		}
		prot1trace.push_back(tmp1);
		prot2trace.push_back(tmp2);
	}

	if (prot1trace.size() != prot2trace.size())
	{
		cout << " proteins dont contain same number of residues. aborting DME calc" << endl;
		return -1.0;
	}
	for (UInt i = 0; i < prot1trace.size(); i ++)
	{
		if (prot1trace[i].size() != prot2trace[i].size())
		{
			cout << " proteins dont contain same number of residues. aborting DME calc" << endl;
			return -1.0;
		}
	}


	UInt N = 0;
	double dij = 0.0;
	for (UInt i = 0; i < prot1trace.size() - 1; i ++) 
	{
		for (UInt j = i + 1; j < prot1trace.size(); j ++)
		{
			for (UInt k = 0; k < prot1trace[i].size(); k ++)
			{
				for (UInt l = 0; l < prot1trace[j].size(); l ++)
				{
					N++;
					dblVec diff1 = prot1trace[i][k] - prot1trace[j][l];
					dblVec diff2 = prot2trace[i][k] - prot2trace[j][l];

					dij += fabs(CMath::dotProduct(diff1, diff1) - CMath::dotProduct(diff2,diff2));
				}
			}
		}
	}
	double DME = sqrt( (2.0 /(double(N)*(double(N)-1.0))) * dij);

	return DME;
} 

void mapGCN4DimerSequence (protein* _prot)
{

    enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
    int _resID[] = {L,E,D,K,V,E,E,L,L,S,K,N,Y,H,L,E,N,E,V,A,R,L,K,K,L};
    int arraySize = sizeof(_resID)/sizeof(_resID[0]);
    cout << arraySize << endl;

    for (int i = 0; i < arraySize; i ++)
    {
        _prot->mutate(0, i, 0);
        _prot->mutate(0, i, (UInt)_resID[i]);
    }

}

