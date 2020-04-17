#include "protein.h"
#include "pdbWriter.h"
#include <iomanip>
#include <iostream>
#include <fstream>


unsigned int pdbWriter(protein* _pProtein, const string& _pdbFile)
{	ofstream outFile( _pdbFile.c_str());
        if(!outFile)
	{	cout << "ERROR: cannot write to file: " << _pdbFile << endl;
                return 0;
        }

		_pProtein->updateResiduesPerTurnType();
		renumberAtoms(_pProtein);
        atomIterator theIterator(_pProtein);
	for (; !(theIterator.last()); theIterator++)
	{
		atom* pCurrentAtom = theIterator.getAtomPointer();
		residue* pCurrentResidue = theIterator.getResiduePointer();
		chain* pCurrentChain = theIterator.getChainPointer();
		outFile << "ATOM";
		outFile << " ";

	// output the serial number
		outFile.width(6);
		outFile << pCurrentAtom->getSerialNumber();
		outFile << " ";

	// output the atom name
		int len = pCurrentAtom->getName().size();
		if (len == 1)
		{
			outFile << " ";
			outFile << pCurrentAtom->getName();
			outFile << "  ";
		}
		if (len == 2)
		{
			outFile << " ";
			outFile << pCurrentAtom->getName();
			outFile << " ";
		}
		if (len == 3)
		{
			outFile << " ";
			outFile << pCurrentAtom->getName();
		}
		if (len == 4)
		{
//			outFile.width(4);
			outFile << pCurrentAtom->getName();
		}
		outFile << " ";

	// output the residue name converting terminal and d-amino acid types to standard amino acid nomenclature
		UInt resTypeIndex = pCurrentResidue->getTypeIndex();
		if (pCurrentResidue->isDNterm(resTypeIndex) || pCurrentResidue->isLNterm(resTypeIndex) || pCurrentResidue->isGNterm(resTypeIndex)){resTypeIndex = resTypeIndex-Nterm;}
		if (pCurrentResidue->isDCterm(resTypeIndex) || pCurrentResidue->isLCterm(resTypeIndex) || pCurrentResidue->isGCterm(resTypeIndex)){resTypeIndex = resTypeIndex-Cterm;}
		if (pCurrentResidue->isD(resTypeIndex)){resTypeIndex = resTypeIndex-Daa;}
		outFile << pCurrentResidue->getType(resTypeIndex);
		outFile << " ";

	// output the chain ID
		outFile << pCurrentChain->getChainID();
		outFile << " ";

	// output the residue number
		outFile.width(3);
		outFile << pCurrentResidue->getResNum();
		outFile << "    ";

	// output the coordinates
		dblVec coords = pCurrentAtom->getCoords();
		for (unsigned int i=0; i<3; i++)
		{
			if (coords[i] >= 100.0)
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
				outFile << setprecision(6) << coords[i];
			}
			if (coords[i] <= -100.0)
			{
				outFile.width(7);
				outFile.setf(ios::showpoint);
				outFile << setprecision(6) << coords[i];
			}
			if ((coords[i] < 100.0 && coords[i] >= 10.0) || (coords[i] > -100.0 && coords[i] <= -10.0))
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
            	     outFile << setprecision(5) << coords[i];
			}
			if ((coords[i] < 10.0 && coords[i] >= 1.0) || (coords[i] > -10.0 && coords[i] <= -1.0))
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
                    outFile << setprecision(4) << coords[i];
			}
			if ((coords[i] < 1.0 && coords[i] >= 0.1) || (coords[i] > -1.0 && coords[i] <= -0.1))
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
                    outFile << setprecision(3) << coords[i];
			}
			if ((coords[i] < 0.1 && coords[i] >= 0.01) || (coords[i] > -0.1 && coords[i] <= -0.01))
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
                    outFile << setprecision(2) << coords[i];
			}
			if ((coords[i] < 0.01 && coords[i] > 0.0) || (coords[i] > -0.01 && coords[i] < 0.0))
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
                    outFile << setprecision(1) << coords[i];
			}
			if (coords[i] == 0.0)
			{
				outFile << " ";
				outFile.width(7);
				outFile.setf(ios::showpoint);
				outFile << setprecision(4) << coords[i];
			}
		}
		outFile << " ";

		// output the dielectric in b-factor
		double dielectric = pCurrentAtom->getDielectric();
		if (dielectric < 10.0)
		{
			outFile << " ";
			outFile.width(4);
			outFile.setf(ios::showpoint);
			outFile << setprecision(3) << dielectric;
		}
		else
		{
			outFile.width(4);
			outFile.setf(ios::showpoint);
			outFile << setprecision(4) << dielectric;
		}

		//output residues per turn
		double RPT = pCurrentAtom->getRPTType();
		outFile.width(6);
		outFile << RPT;
		outFile << endl;
	}
	outFile.close();
	return 1;
}

void renumberAtoms(protein* _pProtein)
{ 
        atomIterator theIterator(_pProtein);
	atom* pCurrentAtom = theIterator.getAtomPointer();
	unsigned int Counter = pCurrentAtom->getSerialNumber();
        for (; !(theIterator.last()); theIterator++)
        {	
		atom* pCurrentAtom = theIterator.getAtomPointer();
		pCurrentAtom->setSerialNumber(Counter);
		Counter++;
	}
	return;
}
