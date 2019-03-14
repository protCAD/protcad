//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  protMutator 1.5  *************************************************
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
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

double getAverageDielectric(protein* _bundle, UInt _resIndex);
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "protMutator <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    srand (getpid());
    UInt numres = bundle->getNumResidues(0);
    UInt numchains = bundle->getNumChains();
    //delete thePDB;
	
    //--parameters
    //int chains[] = {0};
    //int chainsSize = sizeof(chains)/sizeof(chains[0]);
    //int residues[] = {1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};//{46,47,61,92,95};//61,30,9,129/46,47,61,92,95/,87,91,110,150,152{1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};
    //int residuesSize = sizeof(residues)/sizeof(residues[0]);
    //int resID[] = {Csf}, count = 0;
    //UInt resIDsize = sizeof(resID)/sizeof(resID[0]);
    //double Energy, bestE=1E100, chi2;

	//--Mutations
    //for (UInt h = 0; h < resIDsize; h++)
    /*{
		for (UInt l = 0 ; l < 1; l++)
        {
            for (UInt i = 0; i < numres; i++)
            {

                UInt restype = bundle->getTypeFromResNum(l,i);
				if (restype == dK)
                {
					bundle->mutate(l, i, dR);
                }
				if (restype == K)
				{
					bundle->mutate(l, i, R);
				}
               /*if (restype == A)
                {
                    bundle->activateForRepacking(l, i);
                    bundle->mutate(l, i, Y);
                    UIntVec allowedRots = bundle->getAllowedRotamers(l, i, D, 0);
                    bundle->setRotamerWBC(l, i, 0, allowedRots[0]);
                    chi2 = bundle->getChi(l,i,0,1);
                    bundle->setChi(l,i,0,0,chi2+20);
                    //bundle->setChi(l,i,0,0,95);
                   // bundle->activateForRepacking(0, 4);
                    //bundle->mutate(0, 7, dC);
                   // bundle->setChi(0,7,0,0,-60);
                }
                if (restype == Csf || restype == Cf)
                {
                    cout << i << " " << bundle->getChi(l,i,0,0) << endl;
                }
                /*UIntVec allowedRots = bundle->getAllowedRotamers(l, i, resID[0], 0);
                if (i == 12)
                {
                    bundle->setChi(l,i,0,0,-160);
                }
                if (i == 16)
                {
                    bundle->setChi(l,i,0,0,-80);
                }
                if (i == 60)
                {
                    bundle->setChi(l,i,0,0,-160);
                }
                if (i == 64)
                {
                    bundle->setChi(l,i,0,0,-80);
                }

                     //UIntVec allowedRots = bundle->getAllowedRotamers(0, 16, resID[0], 0);
                    //for (int j = -80; j < -40; j++)
                    /*{
                        count++;
                        bundle->setChi(0,12,0,0,-170);
                        bundle->setChi(0,60,0,0,-170);

                            //bestE = Energy;
                            //count++;
                            stringstream convert;
                            string countstr;
                            convert << count, countstr = convert.str();
                            string outFile = countstr + ".good.pdb";
                            pdbWriter(bundle, outFile);
                            //cout << i+l << " " << aminoAcidString[resID[0]] << " " << Energy << endl;
                        //}
                   // }


                //delete thePDB;z*/
           // }
       // }
    //}
    bundle->mutateWBC(0,0,R);
    bundle->mutateWBC(0,2,W);
    bundle->mutateWBC(0,3,dR);
    bundle->mutateWBC(0,5,dR);
    bundle->mutateWBC(0,6,R);
    bundle->mutateWBC(0,8,R);
    bundle->mutateWBC(0,9,dR);
    bundle->mutateWBC(0,11,dR);
    string outFile = infile + ".CRW.pdb";
    pdbWriter(bundle, outFile);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

double getAverageDielectric(protein* _bundle, UInt _resIndex)
{
	UInt count = 0;
	double totalDielectric = 0.0, dielectric;
	for (UInt i = 0; i < _bundle->getNumChains(); i++)
	{
		for (UInt j = 0; j < _bundle->getNumAtoms(i, _resIndex); j++)
		{
			dielectric = _bundle->getDielectric(i,_resIndex,j);
			totalDielectric += dielectric;
			count++;
		}
	}
	return totalDielectric/count;
} 

