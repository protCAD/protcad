#include "residue.h"
#include "atom.h"
#include "pdbData.h"
#include <vector>

int main()
{	residue* tempRes = new residue("ALA");
	delete tempRes;
	cout << "Initialization finished... " << endl;
	cout << "# residues built: " << residue::getHowMany() << endl;
	cout << "# atoms build: " << atom::getHowMany() << endl;
	cout << endl;
	vector<residue*>* temp = new vector<residue*>(0);

	cout << " starting to build residues " << endl;
	temp->push_back(new residue("ALA"));
/*	temp->push_back(new residue("ARG"));
	temp->push_back(new residue("ASN"));
	temp->push_back(new residue("ASP"));
	temp->push_back(new residue("CYS"));
	temp->push_back(new residue("GLN"));
	temp->push_back(new residue("GLU"));
	temp->push_back(new residue("GLY"));
	temp->push_back(new residue("HIS"));
	temp->push_back(new residue("ILE"));
	temp->push_back(new residue("LEU"));
	temp->push_back(new residue("LYS"));
	temp->push_back(new residue("MET"));
	temp->push_back(new residue("PHE"));
	temp->push_back(new residue("PRO"));
	temp->push_back(new residue("SER"));
	temp->push_back(new residue("THR"));
	temp->push_back(new residue("TRP"));
	temp->push_back(new residue("TYR"));
	temp->push_back(new residue("VAL"));
	temp->push_back(new residue("BHA"));
*/
	temp->push_back(new residue( *(*temp)[0] ));

	cout << "# residues built: " << residue::getHowMany() << endl;
	cout << "# atoms build: " << atom::getHowMany() << endl;

	for (unsigned int i=0;i<temp->size();i++)
	{
		cout << (*temp)[i]->getType() << endl;
		(*temp)[i]->printMainChain();
		(*temp)[i]->printBranchPoints();
		cout << "bpt 1 : " << endl;
		cout << "Chi1 " << (*temp)[i]->getChi(0,0) << endl;
		cout << "Chi2 " << (*temp)[i]->getChi(0,1) << endl;
		cout << "Chi3 " << (*temp)[i]->getChi(0,2) << endl;
		cout << "Chi4 " << (*temp)[i]->getChi(0,3) << endl;
		cout << "bpt 2 : " << endl;
		cout << "Chi1 " << (*temp)[i]->getChi(1,0) << endl;
		cout << "Chi2 " << (*temp)[i]->getChi(1,1) << endl;
		cout << "Chi3 " << (*temp)[i]->getChi(1,2) << endl;
		cout << "Chi4 " << (*temp)[i]->getChi(1,3) << endl;
	}

/*	(*temp)[1]->rotate(1,4,100);
	cout << "after chi1 rotation of ARG" << endl;
	cout << "Chi1 " << (*temp)[1]->getChi(0,0) << endl;

	(*temp)[2]->rotate(1,2,20);
*/


	for(unsigned int i=0;i<temp->size();i++)
	{	delete (*temp)[i];
	}

	delete temp;

	return 0;
}
