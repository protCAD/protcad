#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"
#include <pthread.h>


void *intraChainEnergy(chain* _chain);

int main(int argc, char* argv[])
{

	string inputFile = argv[1];

	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	prot->silenceMessages();

	chain one(*(prot->getChain(0)));
	chain two(*(prot->getChain(1)));

	pthread_t mythread;
	pthread_create(&mythread, NULL, intraChainEnergy, &one);
	
	double energy=two->intraEnergy();

	cout << "main thread energy for two is " << energy << endl;
}

void *intraChainEnergy(chain* _chain)
{
	double energy = _chain->intraEnergy();
	cout << "child thread energy for one is " << energy << endl;
	return;
}
