#include "enums.h"
#include "chain.h"
#include <vector>

int main()
{	chain* pTheChain = new chain;

	for (unsigned int i=0; i<20; i++)
	{	residue* temp= new residue(i);
		pTheChain->add(temp);
		cout << pTheChain->intraEnergy() << endl;
	}
	
	delete pTheChain;
	return 0;
}
