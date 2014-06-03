#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main(int argv, char* args[])
{	if( argv < 2 )
	{	cout << " need to provide the pdb filename " << endl;
		exit(1);
	}

	string filename = args[1];

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}
	cout << "number of residues generated: " << residue::getHowMany() << endl;
	cout << "number of chains generated:" << chain::getHowMany() << endl;

	ensemble* pEnsemble = new ensemble;
	pEnsemble->add(pTheProtein1);
	
	cout << " energy : " << pEnsemble->energy() << endl;

	delete pTheProtein1;

	return 0;
}
