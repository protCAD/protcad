#include "enums.h"
#include "protein.h"
#include <vector>

int main()
{
		protein* pTheProtein = new protein;
        chain* pTheChain = new chain;
		pTheProtein->add(pTheChain);
        for (unsigned int i=0; i<20; i++)
        {
                residue* temp1= new residue(i);
                pTheChain->add(temp1);
        }
        cout << residue::getHowMany() << endl;
        cout << atom::getHowMany() << endl;
        pTheChain->mutate(1,5);

        delete pTheProtein;

	return 0;
}
