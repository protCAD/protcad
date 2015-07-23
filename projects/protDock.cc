//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                *************************************************
//**************************************  protDock 2.0  *************************************************
//**************************************                *************************************************
//*******************************************************************************************************
//*************************** -six-axis interchain docking simulation- **********************************
//*******************************************************************************************************

/////// Just specify a 2 chain infile and it will result in docked pdbs.

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

dblVec carbonCentroid(protein* _prot, UInt _chain);
void joinComplex(protein* _prot);
//void diffuseComplex(protein* _prot);

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protDock <inFile.pdb>" << endl;
		exit(1);
	}
	string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,HCe,HCd};
    PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
    residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--initialize variables for loop
    double Energy, pastEnergy, bestEnergy, rotx, roty, rotz, transx, transy, transz;
    UInt name, totalsize = 250, nobetter, test,buffer;//, acceptedStep;
    name = rand() % 1000000;
    stringstream convert;
    string countstr;
    convert << name, countstr = convert.str();
    string bestFile = countstr + ".tempBest.pdb";
    string outFile;

    //generate mutant
    UInt chainNum = bundle->getNumChains();
    UInt randres, randchain;
    UInt activeResidues[] = {46,61,92,95};
    UInt activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]);
    delete thePDB;

	//--Loop for multiple simulations
    for (int a = 0; a < 100; a++)
	{
        PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
        for (UInt i = 0; i < 3; i++)
        {
            randchain = rand() % chainNum;
            randres = activeResidues[rand() % activeResiduesSize];
            bundle->activateForRepacking(randchain, randres);
            bundle->mutateWBC(randchain, randres, He);
        }
        bundle->protOptSolvent(200);
        pdbWriter(bundle, bestFile);
        joinComplex(bundle);
        pastEnergy = bundle->interSoluteEnergy(true, 0, 1);
        nobetter = 0, test = 1, buffer = 0, bestEnergy = 1E10; //acceptedStep = 0,
		//--Run optimizaiton loop till grand minimum---------------------------------------------------
		do
		{  
			//--Move randomly constrained by energy and distance condition
            nobetter++, buffer++;
            rotx = (rand() % 7)-3, roty = (rand() % 7)-3, rotz = (rand() % 7)-3;
            transx = (rand() % 3)-1, transy = (rand() % 3)-1, transz = (rand() % 3)-1;

            do
            {
                //--Make random motion until Energy rises-----------

                bundle->rotateChain(0, X_axis, rotx), bundle->rotateChain(0, Y_axis, roty), bundle->rotateChain(0, Z_axis, rotz);
                bundle->translateChain(0, transx, 0, 0), bundle->translateChain(0, 0, transy, 0), bundle->translateChain(0, 0, 0, transz);
                test++;

                //--Get new distance and Energy
                Energy = bundle->interSoluteEnergy(true, 0, 1);
                if ((Energy < pastEnergy + (buffer/5)) && Energy < 0)
                {
                    //cout << Energy << endl;
                    if (Energy < bestEnergy)
                    {
                        bestEnergy = Energy;
                        pdbWriter(bundle, bestFile);
                    }
                    pastEnergy = Energy;
                    buffer = 0, nobetter--, test = 0;
                    /*acceptedStep++;
                    stringstream convert;
                    string countstr;
                    convert << acceptedStep, countstr = convert.str();
                    outFile = countstr + ".pdb";
                    pdbWriter(bundle, outFile);*/
                }

                //--Revert if no improvement
                if (test != 0)
                {
                    bundle->translateChain(0, 0, 0, (transz * -1)), bundle->translateChain(0, 0, (transy * -1), 0), bundle->translateChain(0, (transx * -1), 0, 0);
                    bundle->rotateChain(0, Z_axis, (rotz * -1)), bundle->rotateChain(0, Y_axis, (roty * -1)), bundle->rotateChain(0, X_axis, (rotx * -1));
                }
            }while (test == 0);
            if (nobetter == totalsize)
            {
                bundle->protOptSolvent(50);
            }
		}while (nobetter < (totalsize * 1.5));
        delete thePDB;
		
		//--Print final energy and write a pdb file----------------------------------------------------
        PDBInterface* theModelPDB = new PDBInterface(bestFile);
        ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
        molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
        protein* model = static_cast<protein*>(modelMol);
        name = rand() % 1000000;
        stringstream convert;
        string countstr;
        convert << name, countstr = convert.str();
        outFile = countstr + ".dock.pdb";
        pdbWriter(model, outFile);
        cout << name << " " << bundle->interSoluteEnergy(true, 0, 1) << endl;
        delete theModelPDB;
	}
	cout << "Complete" << endl << endl;
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///// get average coordinates of all carbons in chain ///////////////////////////////////////////////////
dblVec carbonCentroid(protein* _prot, UInt _chain)
{
	//--initialize and clear variables
	double number = 0;
	UInt numRes, numAtoms;
	string atomType;
	dblVec coords, coordsSum(3), coordsAve(3);
	coordsSum[0] = 0.0, coordsSum[1] = 0.0, coordsSum[2] = 0.0;
	coordsAve[0] = 0.0, coordsAve[1] = 0.0, coordsAve[2] = 0.0;
	
	//--loop through all atoms of all residues in search of carbon
	numRes = _prot->getNumResidues(_chain);
	for (UInt i = 0; i < numRes; i++)
	{
		numAtoms = _prot->getNumAtoms(_chain, i);
		for (UInt j = 0; j < numAtoms; j++)
		{
			atomType = _prot->getTypeStringFromAtomNum(_chain, i, j);
			if (atomType == "C")
			{
				number++;
				coords = _prot->getCoords(_chain, i, j);
				coordsSum[0] += coords[0];
				coordsSum[1] += coords[1];
				coordsSum[2] += coords[2];
			}
		}
	}

	//--get average of all carbon coordinates
	coordsAve[0] = ((coordsSum[0])/number);
	coordsAve[1] = ((coordsSum[1])/number);
	coordsAve[2] = ((coordsSum[2])/number);

	return coordsAve;
}

///translate along center of mass vector till and if interaction energy is not 0
void joinComplex(protein* _bundle)
{
    dblVec substrate, ligandCentroid;
    double dx, dy, dz, mag, nx, ny, nz, energy = 0, dist = 1;
    energy = _bundle->interSoluteEnergy(true, 0, 1);
    if (energy > -1.0)
    {
        do
        {
            substrate =  _bundle->getCoords(0, 47, 54);
            ligandCentroid = _bundle->getCoords(1, 47, 54);
            dx = (ligandCentroid[0]-substrate[0]);
            dy = (ligandCentroid[1]-substrate[1]);
            dz = (ligandCentroid[2]-substrate[2]);
            mag = sqrt(dx*dx+dy*dy+dz*dz);
            nx=dx/mag,ny=dy/mag,nz=dz/mag;
            dx=nx*dist;
            dy=ny*dist;
            dz=nz*dist;
            _bundle->translateChain(0, dx, dy, dz);
            energy = _bundle->interSoluteEnergy(true, 0, 1);
        }
        while (energy > -1.0);
    }
}

/*void diffuseComplex(protein* _prot)
{
	//--Initialize variables for loop
	double Energy, rotx, roty, rotz, transx, transy, transz;
	rotx = (rand() % 90)-45, roty = (rand() % 90)-45, rotz = (rand() % 90)-45;
	do 
	{
		transx = (rand() % 9)-4;
		transy = (rand() % 9)-4;
		transz = (rand() % 9)-4;
	} while (transx == 0 && transy == 0 && transz == 0);
	
	//--Make random motion and test for null interaction
    do
	{	
		_prot->rotateChain(1, X_axis, rotx);
		_prot->rotateChain(1, Y_axis, roty);
		_prot->rotateChain(1, Z_axis, rotz);
		_prot->translateChain(1, transx, 0, 0);
		_prot->translateChain(1, 0, transy, 0);
		_prot->translateChain(1, 0, 0, transz);
        Energy = _prot->bindingSoluteEnergy(0,1);
    } while (Energy != 0.0)
	return;
}*/
