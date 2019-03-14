#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argc, char* argv[])
{	
	if (argc != 3)
	{	cout << "Usage: addMissingSidechains inputpdb outputpdb" << endl;
		exit (1);
	}	
		
	string filename = argv[1];
	string outfilename = argv[2];

	cout << "Inputfile: " << filename << endl;
	cout << "Outputfile: " << outfilename << endl;
	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	
		cout << "Error reading pdbfile " << filename << endl;
		return 1;
	}

	pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
	
	delete pTheProtein1;
	return 0;
}
