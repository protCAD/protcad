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
#include "ensemble.h"
#include "PDBInterface.h"

double getAverageDielectric(protein* _bundle, UInt _resIndex);
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "protRotamer <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hcd,Pch,Csf};
    //string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hcd","Pch","Csf"};
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
	
    //--parameters
    int chains[] = {1};
    int chainsSize = sizeof(chains)/sizeof(chains[0]);
    int residues[] = {95};//61,30,9,129/46,47,61,92,95/,87,91,110,150,152{1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};
    int residuesSize = sizeof(residues)/sizeof(residues[0]);
    //int count = 0;
    //--variable
    UInt chain = 0;
    UInt res = 99;
    bool better;
    double bestE = 1E100, chi1, chi2, chi3;
    //string outFile = "opt." + infile;
    cout << "pdb " << "energy " << endl;
    //UInt randres, randchain, randrot;
    //UIntVec currentRot;
    //UIntVec allowedRots = bundle->getAllowedRotamers(chain,res,Hcd,0);
    //bundle->setChi(chain,res,0,0,-111.592);
    double chi1start = bundle->getChi(chain,res,0,0);
    double chi2start = bundle->getChi(chain,res,0,1);
    //double startEnergy = bundle->protEnergy();

    //for (int h = -30; h < 30; h++)
    //{
        //for (int i = -15; i < 15; i++)
        //{
            better = false;
            //chi1 = chi1start+h;
            //chi2 = chi2start+i;
            //bundle->setChi(chain,res,0,0,chi1);
            //bundle->setChi(chain,res,0,1,chi2);
            //currentRot = bundle->getCurrentRotamer(chain, res);
            //randrot = allowedRots[h];
            //bundle->setRotamerWBC(chain, res, 0, randrot);
            for (UInt l = 0; l < 360; l++)
            {
                chi3 = l;
                bundle->setChi(chain,res,1,0,chi3);
                double Energy = bundle->intraSoluteEnergy(true);
                //if (Energy < bestE)
                //{
                    bestE = Energy;
                    better = true;
                    stringstream convert;
                    string countStr, outFile;
                    convert << l+1, countStr = convert.str();
                    outFile = countStr + ".opt.pdb";
                    pdbWriter(bundle, outFile);
                    cout << l+1 << " " << chi3 << " " << Energy << " improved!" << endl;
                    //break;
                //}
            }
            if (!better)
            {
                bundle->setChi(chain,res,0,0,(chi1*-1));
                bundle->setChi(chain,res,0,1,(chi2*-1));
                bundle->setChi(chain,res,1,0,(chi3*-1));
                //bundle->setRotamerWBC(chain, res, 0, currentRot[0]);
            }
        //}
    //}
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

