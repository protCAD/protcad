#include "pdbReader.h"
#include "atomIterator.h"
#include "annealer.h"
#include "ran.h"

int main (int argc, char* argv[])
{
	string fileName = argv[1];
	molecule* pTheMolecule = pdbReader(fileName);
	cout << fileName << " " <<  pTheMolecule->getVolume(0) << endl;
	return 0;
}
