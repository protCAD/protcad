//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************     triadFiner 1.0    ************************************************
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
        cout << "triadFinder" << endl;
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
    dblVec OHcoords, Hcoords, Ocoords;
    double OHHe, HeO;
	string inFrame, outFile;
    ofstream donors("donors");
    ofstream histidines("histidines");
    ofstream acceptors("acceptors");
    ofstream proteinFile("triads5all");
    double dist = 5;
	//ofstream residueFile("residue");
	//delete bundle;

	//--Run multiple independent evolutions
	proteinFile << "frame time triad interaction" << endl;
    for (UInt a = 1; a < 40000; a++)
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

        //donor
		for (UInt i = 0; i < chainNum; i++)
		{
			resNum = frame->getNumResidues(i);
			for (UInt j = 0; j < resNum; j++)
			{
				restype = frame->getTypeFromResNum(i, j);
                if (restype == S || restype == T)
				{
                    if (restype == S)
                    {
                        OHcoords = frame->getCoords(i, j, "OG");
                    }
                    if (restype == T)
                    {
                        OHcoords = frame->getCoords(i, j, "OG1");
                    }

					//histidine
					for (UInt k = 0; k < chainNum; k++)
					{
						for (UInt l = 0; l < resNum; l++)
						{
							restype = frame->getTypeFromResNum(k, l);
							if (restype == He)
							{
                                Hcoords = frame->getCoords(k, l, "CE1");
                                OHHe = CMath::distance(OHcoords, Hcoords);
                                if (OHHe < dist)
								{
                                    //acceptor
									for (UInt m = 0; m < chainNum; m++)
									{
										for (UInt n = 0; n < resNum; n++)
										{
											restype = frame->getTypeFromResNum(m, n);
                                            if (restype == E || restype == D || n == resNum-1)
											{
                                                if (restype == E)
                                                {
                                                    Ocoords = frame->getCoords(m, n, "CD");
                                                }
                                                if (restype == D)
                                                {
                                                    Ocoords = frame->getCoords(m, n, "CG");
                                                }
                                                if (n == resNum-1)
                                                {
                                                    Ocoords = frame->getCoords(m, n, "C");
                                                }
                                                HeO = CMath::distance(Ocoords, Hcoords);
                                                if (HeO < dist)
												{
													triad++;
                                                    proteinFile << triad;
                                                    donors << j+1 << endl;
                                                    histidines << l+1 << endl;
                                                    acceptors << n+1 << endl;
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
            donors << triad << endl;
            histidines << triad << endl;
            acceptors << triad << endl;
		}
		if (triad > 0)
		{
            outFile = "triad5more_" + countstr + ".pdb";
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

	
