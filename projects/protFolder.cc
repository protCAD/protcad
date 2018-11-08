//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protFolder     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
#include <unistd.h>

UInt protSampling(protein* _prot);
int main (int argc, char* argv[])
{
	if (argc !=2)
	{	cout << "protFolder <inFile.pdb>" << endl;
		exit(1); }

	//--Assign model nomenclature
	string infile = argv[1];
	stringstream convert;
	string startstr;
	srand (getpid());
	UInt name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string foldModel = startstr + "_fold.pdb";
	
	//--Build protein object and calculate starting energy
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	double Energy = 1E10;
	double newEnergy = Energy;
	
	//--sample folds and save energetic minima
	while (true)
	{
		UInt misses = protSampling(prot);
		newEnergy = prot->protEnergy();
		cout << startstr << " " << prot->getNumResidues(0)*5 << " " << misses << " " << newEnergy;
		if (newEnergy < Energy){
			pdbWriter(prot, foldModel);
			cout << " minima!" << endl;
			Energy = newEnergy;
		}
		else{cout << endl;}
	}
	return 0;
}

UInt protSampling(protein* _prot)
{
	//--Initialize variables for loop, calculate starting energy and build energy vectors---------------
	UInt randchain, randres, totalres = 0, chainNum = _prot->getNumChains(), foldD, numres, type, GisL;
	UInt alphaBias, startNumClashes = _prot->getNumHardClashes(), numClashes, sampled = 0, count = 0;
	double sPhi, sPsi;
	
	for (UInt i = 0; i < chainNum; i++)
	{
		totalres+= _prot->getNumResidues(i);
	}
	
	//--Run optimizaiton loop to relative minima, determined by _plateau----------------------------
	do
	{   //--choose random residue
		randchain = rand() % chainNum;
		numres = _prot->getNumResidues(randchain);
		randres = rand() % numres;
		count++;
	
		//--Backbone sampling-----------------------------------------------------------------------
		if (randres != 0 && randres != numres-1)
		{
			// Get angles and types from res
			sPhi = _prot->getPhi(randchain,randres);
			sPsi = _prot->getPsi(randchain,randres);
			type = _prot->getTypeFromResNum(randchain,randres);
			foldD = rand() % 2;
	
			// Correct D-amino acid handedness if needed
			if (type > 26 && type < 53 && sPhi < 0){
				_prot->setDihedral(randchain,randres,sPhi*-1,0,foldD);
				_prot->setDihedral(randchain,randres,sPsi*-1,1,foldD);
			}
	
			// Randomly flip glycine handedness
			GisL = rand() % 2;
			if (type == 26 && sPhi < 0 && GisL == 1){
				_prot->setDihedral(randchain,randres,sPhi*-1,0,foldD);
				_prot->setDihedral(randchain,randres,sPsi*-1,1,foldD);
			}
			if (type == 26 && sPhi > 0 && GisL == 0){
				_prot->setDihedral(randchain,randres,sPhi*-1,0,foldD);
				_prot->setDihedral(randchain,randres,sPsi*-1,1,foldD);
			}
	
			//--flip secondary structure channels
			alphaBias = rand() % 2;
			if (sPhi < 0) // right handed
			{
				if (sPsi < 60 /*&& alphaBias < 1*/){  // flip alpha to beta or polyproline
					_prot->setDihedral(randchain,randres,sPsi+180,1,foldD);
				}
				if (sPsi >= 60){
					if (alphaBias < 1){
						if (sPhi < -90) {  // flip beta to polyproline
							_prot->setDihedral(randchain,randres,sPhi+90,0,foldD);
						}
						else{ // flip polyproline to beta
							_prot->setDihedral(randchain,randres,sPhi-90,0,foldD);
						}
					}
					else{_prot->setDihedral(randchain,randres,sPsi-180,1,foldD);} // flip beta or polyproline to alpha
				}
			}
			else{  // left handed
				if (sPsi > -60 /*&& alphaBias < 1*/){  // flip alpha to beta or polyproline
					_prot->setDihedral(randchain,randres,sPsi-180,1,foldD);
				}
				if (sPsi <= -60){
					if (alphaBias < 1){
						if (sPhi > 90) {  // flip beta to polyproline
							_prot->setDihedral(randchain,randres,sPhi-90,0,foldD);
						}
						else{ // flip polyproline to beta
							_prot->setDihedral(randchain,randres,sPhi+90,0,foldD);
						}
					}
					else{_prot->setDihedral(randchain,randres,sPsi+180,1,foldD);} // flip beta or polyproline to alpha
				}
			}
			
			//--test for acceptable confirmation and if so optimize further
			numClashes = _prot->getNumHardClashes();
			if (numClashes <= startNumClashes)
			{
				sampled++;
				_prot->protOpt(true);
			}
			else{
				_prot->setDihedral(randchain,randres,sPhi,0,foldD);
				_prot->setDihedral(randchain,randres,sPsi,1,foldD);
			}
		}
	} while (sampled < totalres*5);
	return count;
}
