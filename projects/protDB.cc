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
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Hca,Oec};
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf","dCf"};
	residue::setCutoffDistance(10.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--load data
    /*vector < vector <double> > data;
	ifstream infile( "sugarCoords" );
	while (infile)
	{
	    string s;
	    if (!getline( infile, s )) break;

	    istringstream ss( s );
        vector <string> record;

	    while (ss)
	    {
		 string s;
         if (!getline( ss, s, '.' )) break;
		 double f;
  		 ss >> f;
		 record.push_back( f );
	    }

	    data.push_back( record );
	    //cout << record[2] << endl;
    }*/

	//--initialize variables for loop
	//UInt chainNum, resNum;
    UInt counter = 0;
    //double dist;
	string inFrame;
	UInt totalres = 0;
    //int skip = 0, name;
    //string startstr;
	DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");
    vector <int> count(28);
    fill(count.begin(), count.end(),0);
    cout << "A R N D Dh C Cx Cf Q E Eh Hd He Hn Hp I L K M F P O S T W Y V G" << endl;
	while ((pent=readdir(pdir)))
	{ 
		inFrame = pent->d_name;
		if (inFrame.find(".his") != std::string::npos)
		{
            counter++;       
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
            for (UInt i = 0; i < frame->getNumChains(); i++)
            {
                for (UInt j = 0; j < frame->getNumResidues(i); j++)
                {
                    totalres++;
                    UInt restype = frame->getTypeFromResNum(i,j);
                    count[restype] = count[restype]+1;
                }
            }
            cout << inFrame << " ";
            for (UInt i = 0; i < count.size(); i++)
            {
                cout << count[i] << " ";
            }
            cout << endl;
			fill(count.begin(), count.end(),0);
			delete theFramePDB;
		}	
	}
	closedir(pdir);
	/*cout << "A R N D Dh C Cx Cf Q E Eh Hd He Hn Hp I L K M F P O S T W Y V G" << endl;
    for (UInt i = 0; i < count.size(); i++)
    {
        cout << count[i] << " ";
    }
	cout << endl << totalres;*/
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
