#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"

int main (int argc, char* argv[])
{
	string fileName = argv[1]; 
	molecule* pTheMolecule = pdbReader(fileName);
	protein* pTheProtein = static_cast<protein*>(pTheMolecule);
	cout << "protein read in" << endl;
	if (pTheProtein == 0) return 1;


	return 0;
}
