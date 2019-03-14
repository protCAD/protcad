#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);
	string filename = "temp.pdb";

	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}

    atomIterator it1(pTheProtein1);

	double shortestDistance = 20.0;
	int numshort = 0;
	atom* pShortest = 0;

    for (;!(it1.last());it1++)
    {
		atom* pAtom1 = it1.getAtomPointer();
        atomIterator it2(pTheProtein1);
		atom* pAtom2 = 0
		double tempDistance;
        for (;!(it2.last());it2++)
        {
            pAtom2 = it2.getAtomPointer();
			tempDistance = pAtom1->distance(pAtom2);
            if (tempDistance < shortestDistance && pAtom1 != pAtom2)
			{ 
			   shortestDistance = tempDistance; 
			   pShortest = pAtom2;
			}
        }
		cout <<  pAtom1->getName() << "  " << pAtom2->getName() << "  " << shortestDistance << endl;
    }

	delete pTheProtein1;

	time(&endTime);
	cout << "startTime = " << startTime << "\n";
	cout << "endTime = " << endTime << "\n";
	runTime = endTime - startTime;
	cout << "runTime = " << runTime << "\n";

	return 0;
}
