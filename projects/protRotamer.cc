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
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
    //string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hcd","Pch","Csf"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    residue::setCutoffDistance(8.0);
    residue::setElectroSolvationScaleFactor(0.0);
    residue::setHydroSolvationScaleFactor(0.0);
	amberElec::setScaleFactor(1.0);
    amberVDW::setScaleFactor(1.0);
    srand (time(NULL));
	
	UInt _chainIndex = 0;
	UInt randres = 7;
    UInt bestrot;
	UInt count = 0;

    double pastEnergy = 1E100, Energy;
    //--Get current rotamer and allowed

	//double chi = bundle->getChi(_chainIndex,randres,0,0);
	//double chi2 = bundle->getChi(_chainIndex,randres,0,2);
	//bundle->setChi(_chainIndex,randres,0,0,-142.9);
	//bundle->setChi(_chainIndex,randres,0,0,chi+20);
	//bundle->setChi(_chainIndex,randres,0,2,chi2+10);
	//bundle->mutateWBC(_chainIndex,randres,dDh);
	UInt randrestype = bundle->getTypeFromResNum(_chainIndex,randres);
	//bundle->mutateWBC(2,randres,dD);
	vector <UIntVec> allowedRots = bundle->getAllowedRotamers(_chainIndex, randres, randrestype);
    UInt b = 0;
    for (UInt j = 0; j < allowedRots[b].size(); j++)
    {
        count++;

		bundle->setRotamerWBC(_chainIndex, randres, b, allowedRots[b][j]);
		//bundle->setRotamerWBC(2, randres, b, allowedRots[b][j]);
        Energy = bundle->protEnergy();
        stringstream convert;
        string countstr;
        convert << count, countstr = convert.str();
		string outFile = countstr + ".rot.pdb";
        pdbWriter(bundle, outFile);
        cout << count << " " << Energy << endl;
        if (Energy < pastEnergy)
        {
            bestrot = j;
            pastEnergy = Energy;
        }
	}
	/*double betachi = bundle->getBetaChi(_chainIndex, randres);
	for (int j = -10; j < 0; j++)
	{
		count++;
		bundle->setBetaChi(_chainIndex, randres, betachi+j);
		Energy = bundle->protEnergy();
		stringstream convert;
		string countstr;
		convert << count, countstr = convert.str();
		string outFile = countstr + ".rot.pdb";
		pdbWriter(bundle, outFile);
		cout << count << " " << Energy << endl;
		if (Energy < pastEnergy)
		{
			bestrot = j;
			pastEnergy = Energy;
		}
	}*/
    //cout << infile << " " << pastEnergy << endl;
	//bundle->setRotamerWBC(_chainIndex, randres, 0, allowedRots[bestrot]);
	//pdbWriter(bundle, infile);
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

