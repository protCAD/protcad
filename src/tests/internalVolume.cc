#include <pdbReader.h>
#include <pdbWriter.h>
#include <atomIterator.h>
#include <ran.h>

int main (int argc, char* argv[])
{
	if (argc == 0) 
	{
		cout << "Please enter filename" << endl;
		return 1;
	}

	molecule* pTheMolecule = pdbReader(argv[1]);
	if (pTheMolecule == 0)
	{
		cout << "protein file empty" << endl;
		return 1;
	}

	protein* pTheProtein = static_cast<protein*>(pTheMolecule);

	cout << "number of chains: " << chain::getHowMany() << endl;
	cout << "number of residues: " << residue::getHowMany() << endl;

	if (argc < 2)
	{
		cout << "need active residues ... list residue numbers separated by a space." << endl;
	vector <UInt> residueList;
	for (UInt i = 2; i <= argc  
