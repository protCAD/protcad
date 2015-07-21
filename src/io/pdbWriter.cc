#include "protein.h"
#include "ligand.h"
#include "pdbWriter.h"
//#include "svmt.h"
#include <iomanip>
#include <iostream>
#include <fstream>


unsigned int pdbWriter(vector<protein*> _protVec, vector<ligand*> _ligVec, const string& _pdbFile)
{
    ofstream outFile(_pdbFile.c_str());
    outFile.close();  // These two lines to solve an appending problem
    outFile.close(); 
    for(UInt i=0; i<_protVec.size(); i++)
    {
        pdbWriter(_protVec[i],_pdbFile);
    }

    for (UInt j=0; j<_ligVec.size(); j++)
    {
        pdbWriter(_ligVec[j],_pdbFile);
    }
    return 0;
}

unsigned int pdbWriter(protein* _pProtein, const string& _pdbFile)
{	ofstream outFile( _pdbFile.c_str());
        if(!outFile)
	{	cout << "ERROR: cannot write to file: " << _pdbFile << endl;
                return 0;
        }

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
	//	outFile << pCurrentResidue->getResNum();
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
	// output the residue name
		outFile << pCurrentAtom->getResType();
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

		// output solvation energy
		double solvationEnergy = pCurrentAtom->getSolvationEnergy();
		if ((solvationEnergy > -10.0 && solvationEnergy <= -1.0) || (solvationEnergy < 10.0 && solvationEnergy >= 1.0))
		{
			outFile.width(6);
			outFile.setf(ios::showpoint);
			   outFile << setprecision(4) << solvationEnergy;
		}
		else if ((solvationEnergy > -1.0 && solvationEnergy <= -0.1) || (solvationEnergy < 1.0 && solvationEnergy >= 0.1))
		{
			outFile.width(6);
			outFile.setf(ios::showpoint);
			   outFile << setprecision(3) << solvationEnergy;
		}
		else if ((solvationEnergy > -0.1 && solvationEnergy <= -0.01) || (solvationEnergy < 0.1 && solvationEnergy >= 0.01))
		{
			outFile.width(6);
			outFile.setf(ios::showpoint);
			   outFile << setprecision(2) << solvationEnergy;
		}
		else if ((solvationEnergy > -0.01 && solvationEnergy > -0.001) || (solvationEnergy < 0.01 && solvationEnergy > 0.001))
		{
			outFile.width(6);
			outFile.setf(ios::showpoint);
			  outFile << setprecision(1) << solvationEnergy;
		}
		outFile << endl;
	}
	outFile.close();
	return 1;
}

unsigned int pdbWriter(ligand* _pLigand, const string& _pdbFile)
{       
	ofstream outFile( _pdbFile.c_str(),ios::app);
        if(!outFile)
        {       cout << "ERROR: cannot write to file: " << _pdbFile << endl;
                return 0;
        }

	//cout <<"This is the pdbWriter(ligand*) function." << endl;
	//cout <<"Atoms in ligand= " <<_pLigand->atomCount() << endl;

	renumberAtoms(_pLigand);

        atomIterator theIterator(_pLigand);

        for (; !(theIterator.last()); theIterator++)
        {
                atom* pCurrentAtom = theIterator.getAtomPointer();
                //residue* pCurrentResidue = theIterator.getResiduePointer();
                //chain* pCurrentChain = theIterator.getChainPointer();
                outFile << "HETATM";
                
        // output the serial number
                outFile.width(5);
                outFile << pCurrentAtom->getSerialNumber();
        //	outFile << pCurrentResidue->getResNum();
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
//                      outFile.width(4);
                        outFile << pCurrentAtom->getName();
                }

                outFile << " ";
        // output the residue name
                string tempResType= pCurrentAtom->getResType();
                outFile <<tempResType;
                
                for(UInt i=tempResType.size(); i<4;i++){outFile << " ";}
                //this for loop takes care of resNames less than size 4
                
                outFile << "    ";
        // output the chain ID
		outFile << pCurrentAtom->getLigChainID();
		UInt tempSize=(pCurrentAtom->getLigChainID()).size();
		//cout << "number of spaces to add= " << 8-tempSize  << endl;
		for(UInt i=0; i<5-tempSize;i++){outFile << " ";}
		
               
        // output the residue number
                //outFile.width(3);
                //outFile << pCurrentResidue->getResNum();
                //outFile << "   ";
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
                if ((coords[i] < 100.0 && coords[i] >= 10.0) ||
                    (coords[i] > -100.0 && coords[i] <= -10.0))
 			{
                        outFile << " ";
                        outFile.width(7);
                        outFile.setf(ios::showpoint);
                        outFile << setprecision(5) << coords[i];
                        }
                if ((coords[i] < 10.0 && coords[i] >= 1.0) ||
                    (coords[i] > -10.0 && coords[i] <= -1.0))
                        {
                        outFile << " ";
                        outFile.width(7);
                        outFile.setf(ios::showpoint);
                        outFile << setprecision(4) << coords[i];
                        }
                if ((coords[i] < 1.0 && coords[i] >= 0.1) ||
                    (coords[i] > -1.0 && coords[i] <= -0.1))
                        {
                        outFile << " ";
                        outFile.width(7);
                        outFile.setf(ios::showpoint);
                        outFile << setprecision(3) << coords[i];
                        }
                if ((coords[i] < 0.1 && coords[i] >= 0.01) ||
                    (coords[i] > -0.1 && coords[i] <= -0.01))
                        {
                        outFile << " ";
                        outFile.width(7);
                        outFile.setf(ios::showpoint);
                        outFile << setprecision(2) << coords[i];
                        }
                if ((coords[i] < 0.01 && coords[i] > 0.0) ||
                    (coords[i] > -0.01 && coords[i] < 0.0))
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

		//Add the element to the end of the line (column 76 and 77).
		//This is important b/c the atom name is not recognized by pdbAtomRecord
		//and so it instead takes the atom type from the "element" value.

		for (int i=0; i<21; i++){outFile << " ";}
	
		string typeName;
		typeName=pCurrentAtom->getType();

		if(typeName.size() ==0){ outFile << "   ";}
		if(typeName.size()==1){ outFile << "  ";}
		if(typeName.size()==2){outFile << " ";}
		outFile << typeName;
		
		//Output the charge (if any) to column 78 and 79. Format is +/- in 78 and int in 79
		//Will implement after talking to Vikas.  At present the atom constructor sets the charge to 0.0

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

void renumberAtoms(ligand* _pLigand)
{
        atomIterator theIterator(_pLigand);
        atom* pCurrentAtom = theIterator.getAtomPointer();
        UInt Counter = pCurrentAtom->getSerialNumber();
        for (; !(theIterator.last()); theIterator++)
        {
                atom* pCurrentAtom = theIterator.getAtomPointer();
                pCurrentAtom->setSerialNumber(Counter);
		Counter++;
        }

        return;
}
