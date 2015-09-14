//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************   database_phipsi 1.0   ************************************************
//********************************                       ************************************************
//*******************************************************************************************************
//***************   -Stability Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure and it will evolve in hetero-oligameric stability

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <dirent.h>
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
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--load data
	vector < vector <double> > data;
	ifstream infile( "sugarCoords" );
	while (infile)
	{
	    string s;
	    if (!getline( infile, s )) break;

	    istringstream ss( s );
	    vector <double> record;

	    while (ss)
	    {
		 string s;
		 if (!getline( ss, s, ',' )) break;
		 double f;
  		 ss >> f;
		 record.push_back( f );
	    }

	    data.push_back( record );
	    //cout << record[2] << endl;
	}

	//--initialize variables for loop
	//UInt chainNum, resNum;
	//UInt count = 0;
	double dist;
	dblVec coords, coords1(3);
	string inFrame;
	int skip = 0, name;
	string startstr;
	DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");

	while ((pent=readdir(pdir)))
	{ 
		inFrame = pent->d_name;
		if (inFrame.find(".pdb") != std::string::npos)
		{
			//count++;
			//cout << inFrame << " " << count << endl;
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			frame->silenceMessages();
			skip = 0;
			//cout << "test1" << endl;
			for (UInt i = 541; i < 757; i++)
			{
				//cout << "test2" << endl;
				coords = frame->getCoords(0, i, "CA");
				//cout << "test3" << endl;
				for (UInt j = 0; j < data.size(); j++)
				{
					//cout << "test3.5" << endl;
					coords1[0] = data[j][0];
					coords1[1] = data[j][1];
					coords1[2] = data[j][2];
					//cout << coords1[0] << " " << coords1[1] << " " << coords1[2] << endl;
					dist = CMath::distance(coords1, coords);
					if (dist < 8)
					{
						skip = 1;
						break;
					}
				}
				//cout << i << endl;
			}
			if (skip == 0)
			{
				cout << "good" << endl;
				stringstream convert;
				name = rand() % 1000000;
				convert << name, startstr = convert.str();
				string bestModel = startstr + "_good";
				pdbWriter(frame, bestModel); 
			}
			delete theFramePDB;
		}	
	}
	closedir(pdir);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
