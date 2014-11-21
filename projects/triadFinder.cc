//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************   amberAnalyzer 1.0   ************************************************
//********************************                       ************************************************
//*******************************************************************************************************
//***************   -Stability Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure and it will evolve in hetero-oligameric stability

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"


//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=1)
	{
		cout << "amberAnalyzer" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));

	//--initialize variables for loop
	UInt triad, chainNum, resNum, restype, name;
	dblVec SOcoords, EO1coords, EO2coords, HDcoords, HEcoords;
	double SOHD, SOHE, EO1HD, EO1HE, EO2HD, EO2HE;
	stringstream convert;
	string inFrame, outFile;
    ofstream proteinFile("triad6");
    double dist = 6;
	//ofstream residueFile("residue");
	//delete bundle;

	//--Run multiple independent evolutions
	proteinFile << "frame time triad interaction" << endl;
    for (UInt a = 1; a < 26000; a++)
	{ 
		name = a;
		stringstream convert; 
		string countstr;
		convert << name, countstr = convert.str();
		inFrame = "pC_Frame." + countstr + ".pdb";
		PDBInterface* theFramePDB = new PDBInterface(inFrame);
		ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
		molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
		protein* frame = static_cast<protein*>(frameMol);
		frame->silenceMessages();
		triad = 0, chainNum = frame->getNumChains();
		
		proteinFile << a << " " << (a*.02) << " ";

		//serine
		for (UInt i = 0; i < chainNum; i++)
		{
			resNum = frame->getNumResidues(i);
			for (UInt j = 0; j < resNum; j++)
			{
				restype = frame->getTypeFromResNum(i, j);
				if (restype == S)
				{
					SOcoords = frame->getCoords(i, j, "OG");

					//histidine
					for (UInt k = 0; k < chainNum; k++)
					{
						for (UInt l = 0; l < resNum; l++)
						{
							restype = frame->getTypeFromResNum(k, l);
							if (restype == He)
							{
								HDcoords = frame->getCoords(k, l, "ND1");
					 			HEcoords = frame->getCoords(k, l, "NE2");
								SOHD = CMath::distance(SOcoords, HDcoords);
								SOHE = CMath::distance(SOcoords, HEcoords);
								if (SOHD < dist || SOHE < dist)
								{
									//glutamate
									for (UInt m = 0; m < chainNum; m++)
									{
										for (UInt n = 0; n < resNum; n++)
										{
											restype = frame->getTypeFromResNum(m, n);
											if (restype == E)
											{
												EO1coords = frame->getCoords(m, n, "OE1");
												EO2coords = frame->getCoords(m, n, "OE2");
												EO1HD = CMath::distance(EO1coords, HDcoords);
												EO2HD = CMath::distance(EO2coords, HDcoords);
												EO1HE = CMath::distance(EO1coords, HEcoords);
												EO2HE = CMath::distance(EO2coords, HEcoords);
												if (EO1HD < dist || EO2HD < dist || EO1HE < dist || EO2HE < dist)
												{
													triad++;
													proteinFile << triad << " " << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << " " << m+1 << " " << n+1 << " ";
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}	
		if (triad == 0)
		{
			proteinFile << triad;
		}
		if (triad > 0)
		{
            outFile = "triad6_" + countstr + ".pdb";
			pdbWriter(frame, outFile);
		}
		proteinFile << endl;
		delete theFramePDB;
	}
	cout << "Complete" << endl << endl;
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
