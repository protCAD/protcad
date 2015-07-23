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
    rotamer::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);
    srand (time(NULL));
    delete thePDB;
	
    //--parameters
    int set1[] = {10,30,61,130};//{1,3,21,30,54,59,61,72,78,90,92,110,112,116,129,133,166,174,181};
    int set1Size = sizeof(set1)/sizeof(set1[0]);
    int resID[] = {Hce};

    //--variables
    UInt mutant = resID[0];
    double pastEnergy, Energy, bestAngle, angle;
    vector < vector <double> > currentRot, bestRot;
    UIntVec allowedRots;
    cout << endl << "pdb " << "energy " << endl;

	//--Mutations
    for (int i = 0; i < set1Size; i++)
	{
        for (int j = i+1; j < set1Size; j++)
		{
            for (int k = j+1; k < set1Size; k++)
            {
                #pragma omp parallel for
                for (int l = k+1; l < set1Size; l++)
                {
                    PDBInterface* thePDB = new PDBInterface(infile);
                    ensemble* theEnsemble = thePDB->getEnsemblePointer();
                    molecule* pMol = theEnsemble->getMoleculePointer(0);
                    protein* bundle = static_cast<protein*>(pMol);

                    //--make point mutation(s)
                    bundle->activateForRepacking(0, set1[i]);
                    bundle->activateForRepacking(0, set1[j]);
                    bundle->activateForRepacking(0, set1[k]);
                    bundle->activateForRepacking(0, set1[l]);
                    bundle->mutate(0, set1[i], mutant);
                    bundle->mutate(0, set1[j], mutant);
                    bundle->mutate(0, set1[k], mutant);
                    bundle->mutate(0, set1[l], mutant);

                    int restype = bundle->getTypeFromResNum(0, set1[i]);

                    //--find lowest Rotamer and optimize neighbors
                    for (UInt a = 0; a < 2; a++)
                    {
                        for (UInt b = 0; b < 4; b++)
                        {
                            pastEnergy = 1E10;
                            int res;
                            if (b == 0)
                            {
                                res = set1[i];
                            }
                            if (b == 1)
                            {
                                res = set1[j];
                            }
                            if (b == 2)
                            {
                                res = set1[k];
                            }
                            if (b == 3)
                            {
                                res = set1[l];
                            }
                            allowedRots = bundle->getAllowedRotamers(03, res, restype, 0);
                            pastEnergy = bundle->intraSoluteEnergy(true);
                            for (UInt m = 0; m < allowedRots.size(); m ++)
                            {
                                bundle->setRotamerWBC(0, res, 0, allowedRots[m]);
                                currentRot = bundle->getSidechainDihedrals(0, res);
                                angle = 0;
                                for (UInt n = 0; n < 8; n++)
                                {
                                    angle = angle+45;
                                    bundle->setChi(0, res, 1, 0, angle);
                                    Energy = bundle->intraSoluteEnergy(true);
                                    if (Energy < pastEnergy)
                                    {
                                        pastEnergy = Energy;
                                        bestAngle = angle;
                                        bestRot = currentRot;
                                    }
                                    bundle->setChi(0, res, 1, 0, (angle*-1));
                                }
                            }
                            bundle->setSidechainDihedralAngles(0, res, bestRot);
                            bundle->setChi(0, res, 1, 0, bestAngle);
                        }
                        //bundle->protOptSolvent(200);
                    }

                    //--print pdb and data to output
                    Energy = bundle->intraSoluteEnergy(true);
                    stringstream convert1, convert2, convert3, convert4;
                    string set1Str, set2Str, set3Str, set4Str;
                    convert1 << set1[i]+1, convert2 << set1[j]+1, convert3 << set1[k]+1, convert4 << set1[l]+1;
                    set1Str = convert1.str(), set2Str = convert2.str(), set3Str = convert3.str(), set4Str = convert4.str();
                    string outFile = set1Str + "_" + set2Str + "_" + set3Str + "_" + set4Str +".pdb";
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

