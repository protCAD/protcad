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

//double getAverageDielectric(protein* _bundle, UInt _resIndex);
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
	
	/*UInt _chainIndex = 0;
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
	//UInt randrestype = bundle->getTypeFromResNum(_chainIndex,randres);
	//bundle->mutateWBC(2,randres,dD);
	/*vector <UIntVec> allowedRots = bundle->getAllowedRotamers(_chainIndex, randres, randrestype);
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
	}*/
	UInt count = 0;
	double Ybetachi, Hbetachi, betachiY, betachiH;
	double betachiYS = bundle->getBetaChi(0, 90);
	double betachiHS = bundle->getBetaChi(0, 14);
	 dblVec Ocoords = bundle->getCoords(0, 90, "OH");
     dblVec Ncoords = bundle->getCoords(0, 14, "NE2");
	 dblVec Hcoords = bundle->getCoords(0,90,"HO");
	 dblVec CYcoords = bundle->getCoords(0,90,"CZ");
	 dblVec CHcoords = bundle->getCoords(0,14,"CE1");
    double dist = CMath::distance(Ocoords, Ncoords);
	double Yangle = CMath::angle(CYcoords, Ocoords, Ncoords);
	double Hangle = CMath::angle(CHcoords, Ncoords, Ocoords);
	
	//cout << betachiY << " Y " << betachiH << " H " << dist << " dist " << Yangle << " Yangle " << Hangle << " Hangle " << endl;
	for (int j = -6; j < 6; j++)
	{
		for (int k = -6; k < 6; k++)
		{
			count++;
			Ybetachi = betachiYS+j;
			bundle->setBetaChi(0, 90, Ybetachi);
			Hbetachi = betachiHS+k;
			bundle->setBetaChi(0, 14, Hbetachi);
			betachiY = bundle->getBetaChi(0, 90);
			betachiH = bundle->getBetaChi(0, 14);
			Ocoords = bundle->getCoords(0, 90, "OH");
			Ncoords = bundle->getCoords(0, 14, "NE2");
			Hcoords = bundle->getCoords(0,90,"HO");
			CYcoords = bundle->getCoords(0,90,"CZ");
			CHcoords = bundle->getCoords(0,14,"CE1");
			dist = CMath::distance(Ocoords, Ncoords);
			Yangle = CMath::angle(CYcoords, Ocoords, Ncoords);
			Hangle = CMath::angle(CHcoords, Ncoords, Ocoords);
			cout << count << " " << betachiY << " Y " << betachiH << " H " << dist << " dist " << Yangle << " Yangle " << Hangle << " Hangle " << endl;
			stringstream convert;
			string countstr;
			convert << count, countstr = convert.str();
			string outFile = countstr + ".BC.pdb";
			pdbWriter(bundle, outFile);
		}
	}
    //cout << infile << " " << pastEnergy << endl;
	//bundle->setRotamerWBC(_chainIndex, randres, 0, allowedRots[bestrot]);
	//pdbWriter(bundle, infile);*/
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

