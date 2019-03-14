#include "enums.h"
#include "protein.h"
#include "vanDerWaals.h"
#include "atomIterator.h"
#include <vector>

int main()
{

	chain* pChain = new chain;
	
	for(unsigned int i=0; i < 20 ; i++)
	{	residue* pResidue = new residue(i);
		cout << "residue " << i << " is built " << endl;
		pChain->add(pResidue);
	}

	protein* pProtein = new protein;

	pProtein->add(pChain);

	atomIterator HuiZhang(pProtein);
	
	for (;!(HuiZhang.last());HuiZhang++) 
	{
		atomIterator HuiZhang2(HuiZhang);
		HuiZhang2++;
		for (;!(HuiZhang2.last());HuiZhang2++)
		{
			if (HuiZhang.getResiduePointer())
			cout << (HuiZhang.getResiduePointer())->getType() << "  ";
			if (HuiZhang.getAtomPointer())
			cout << (HuiZhang.getAtomPointer())->getName() << "  ";
			if (HuiZhang2.getResiduePointer())
			cout << (HuiZhang2.getResiduePointer())->getType() << " ";
			if (HuiZhang2.getAtomPointer())
			cout << (HuiZhang2.getAtomPointer())->getName() << endl;;
		}
	}
	
	delete pProtein;

	return 0;
}
