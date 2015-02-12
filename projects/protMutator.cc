//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  protMutator 1.5  ************************************************
//***********************************                    ************************************************
//*******************************************************************************************************
//**************   -point mutations, then backbone and sidechain optimization-   ************************
//*******************************************************************************************************

/////// Just specify a infile and preferred outfile name.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "protMutator <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Hce};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Hce"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    residue::setCutoffDistance(9.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);
    solvation::setItsScaleFactor(0.0);
    srand (time(NULL));
    delete thePDB;
	
    //--parameters
    int chains[] = {0};
	int chainsSize = sizeof(chains)/sizeof(chains[0]);
    int residues[] = {1,3,5,7,9,13,15,17,19,21,28,30,32,34,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,166,168,170,172,174,181,183,185};//{78,116,54,110,133,112,3,1,21,30,166};
	int residuesSize = sizeof(residues)/sizeof(residues[0]);
    int resID[] = {Hce};
	int resIDsize = sizeof(resID)/sizeof(resID[0]);
    int replicates = 1;
	int count = 0;

    //--variables
    UInt mutant, restype;
    double pastEnergy, Energy;
    vector < vector <double> > currentRot;
    UIntVec allowedRots;
    cout << endl << "pdb " << "residue " << "site " << "mutant " << "energy " << endl;

	//--Mutations
	for (int i = 0; i < residuesSize; i++)
	{
		for (int j = 0; j < resIDsize; j++)
		{
			for (int k = 0; k < replicates; k++)
			{
				PDBInterface* thePDB = new PDBInterface(infile);
				ensemble* theEnsemble = thePDB->getEnsemblePointer();
				molecule* pMol = theEnsemble->getMoleculePointer(0);
				protein* bundle = static_cast<protein*>(pMol);

                //--mutate and optimize mutation
				mutant = resID[j];
				for (int l = 0; l < chainsSize; l++)
				{
                    //--make point mutation(s)
					bundle->activateForRepacking(chains[l], residues[i]);
					bundle->mutate(chains[l], residues[i], mutant);
                    restype = bundle->getTypeFromResNum(chains[l], (UInt)residues[i]);

                    //--find lowest Rotamer and optimize neighbors
                    allowedRots = bundle->getAllowedRotamers(chains[l], residues[i], restype, 0);
                    for (UInt a = 0; a < 3; a++)
                    {
                        pastEnergy = bundle->intraSoluteEnergy(true);
                        for (UInt j = 0; j < allowedRots.size(); j ++)
                        {
                            currentRot = bundle->getSidechainDihedrals(chains[l], residues[i]);
                            bundle->setRotamerWBC(chains[l], residues[i], 0, allowedRots[j]);
                            Energy = bundle->intraSoluteEnergy(false);
                            if (Energy < pastEnergy)
                            {
                                pastEnergy = Energy;
                            }
                            else
                            {
                                bundle->setSidechainDihedralAngles(chains[l], residues[i], currentRot);
                            }
                        }
                        bundle->protOptSolvent(200);
                    }
				}

                //--print pdb and data to output
                Energy = bundle->intraSoluteEnergy(false);
				count++;
				stringstream convert; 
				string countstr;
				convert << count, countstr = convert.str();
				string outFile = countstr + ".pdb";
                pdbWriter(bundle, outFile);
                cout << count << " " << aminoAcidString[restype] << " " << residues[i]+1 << " " << aminoAcidString[mutant] << " " << Energy << endl;
                delete thePDB;
			}
		}
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

