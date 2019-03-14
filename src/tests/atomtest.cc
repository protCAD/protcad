#include <iostream>
#include "atom.h"
//#include "svmt.h"

int main()
{	atom* Jason = new atom("N");
     	cout << "Jason's X: " << Jason->getX() << endl;
	cout <<"atom's type " << Jason->getType() << endl;

	delete Jason;
	cout << "No of Atoms Left: " << atom::getHowMany() << endl;

	atom* Bill = new atom("C");
	atom chris;
	cout << "No of Atoms Left: " << atom::getHowMany() << endl;
	chris = *Bill;
	cout << "No of Atoms Left: " << atom::getHowMany() << endl;
	atom gary = *Bill;
	cout << "No of Atoms Left: " << atom::getHowMany() << endl;

	delete Bill;
	return 1;
}
