//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  protNetwork 1.0  *************************************************
//***********************************                    ************************************************
//*******************************************************************************************************
//**********   -point mutation network build, then backbone and sidechain optimization-   ***************
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
        cout << "protNetwork <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Hce};
    //string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Hce"};
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
    int set1[] = {70,72};
    int set2[] = {110,133,112};
    int set3[] = {127,129,116};
    int set4[] = {1,3,21,30,155};
    int set1Size = sizeof(set1)/sizeof(set1[0]);
    int set2Size = sizeof(set2)/sizeof(set2[0]);
    int set3Size = sizeof(set3)/sizeof(set3[0]);
    int set4Size = sizeof(set4)/sizeof(set4[0]);
    int resID[] = {Hce};

    //--variables
    UInt mutant = resID[0];
    double Energy;
    cout << endl << "pdb " << "energy " << endl;

	//--Mutations
    for (int i = 0; i < set1Size; i++)
	{
        for (int j = 0; j < set2Size; j++)
		{
            for (int k = 0; k < set3Size; k++)
			{
                for (int l = 0; l < set4Size; l++)
                {
                    PDBInterface* thePDB = new PDBInterface(infile);
                    ensemble* theEnsemble = thePDB->getEnsemblePointer();
                    molecule* pMol = theEnsemble->getMoleculePointer(0);
                    protein* bundle = static_cast<protein*>(pMol);

                    //--make point mutation(s)
                    bundle->activateForRepacking(0, set1[i]);
                    bundle->activateForRepacking(0, set2[j]);
                    bundle->activateForRepacking(0, set3[k]);
                    bundle->activateForRepacking(0, set4[l]);
                    bundle->mutate(0, set1[i], mutant);
                    bundle->mutate(0, set2[j], mutant);
                    bundle->mutate(0, set3[k], mutant);
                    bundle->mutate(0, set4[l], mutant);

                    //--optimize
                    bundle->protOptSolvent(500);

                    //--print pdb and data to output
                    Energy = bundle->intraSoluteEnergy(true);
                    stringstream convert1, convert2, convert3, convert4;
                    string set1Str, set2Str, set3Str, set4Str;
                    convert1 << set1[i]+1, convert2 << set2[j]+1, convert3 << set3[k]+1, convert4 << set4[l]+1;
                    set1Str = convert1.str(), set2Str = convert2.str(), set3Str = convert3.str(), set4Str = convert4.str();
                    string outFile = set1Str + "_" + set2Str + "_" + set3Str + "_" + set4Str + ".pdb";
                    pdbWriter(bundle, outFile);
                    cout << outFile << " " << Energy << endl;
                    delete thePDB;
                }
			}
		}
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

