//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************      protSorter 1.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//***************   -                                                         -   ***********************
//*******************************************************************************************************


//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

vector < vector < UInt > > buildSequencePool();
//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Hem};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};
    if (argc !=1)
	{
                cout << "protSorter" << endl;
		exit(1);
	}
    /*vector < vector < UInt > > sequencePool;
    sequencePool = buildSequencePool();

    for (UInt i = 0; i < sequencePool.size(); i++)
    {
        for (UInt j = 0; j < sequencePool[i].size(); j++)
        {
            cout << sequencePool[i][j] << " ";
        }
        cout << endl;
    }*/

    residue::setCutoffDistance(8.0);
	residue::setTemperature(300);
    residue::setElectroSolvationScaleFactor(1.0);
    residue::setHydroSolvationScaleFactor(1.0);
    amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);

    /*/create evo data file
    string inFrame;
    DIR *pdir;
    struct dirent *pent;
    pdir=opendir(".");

    //--get sequence evolution results for position
	while ((pent=readdir(pdir)))
	{
		inFrame = pent->d_name;
		if (inFrame.find("XCXXCX.pdb") != std::string::npos)
		{
			PDBInterface* theModelPDB = new PDBInterface(inFrame);
			ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
			molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
			protein* model = static_cast<protein*>(modelMol);
			double phi;
			phi = model->getAngle(0,1,0);
			if (phi < 0)
			{
				phi = model->getAngle(0,2,0);
				if(phi > 0)
				{
					phi = model->getAngle(0,3,0);
					if (phi < 0)
					{
						phi = model->getAngle(0,4,0);
						if (phi > 0)
						{
							phi = model->getAngle(0,5,0);
							if (phi < 0)
							{
								phi = model->getAngle(0,6,0);
								if (phi > 0)
								{
									replace (inFrame.begin(), inFrame.end(),'b', 'h');
									string outfile = inFrame;
									pdbWriter(model, outfile);
								}
							}
						}
					}
				}
			}
			delete theModelPDB;
		}
	}
    closedir(pdir);*/

    //create evo final file
    string inFrame;
    DIR *pdir;
    struct dirent *pent;
    pdir=opendir(".");
    UInt count = 0;

    //--get sequence evolution results for position
    while ((pent=readdir(pdir)))
    {
        inFrame = pent->d_name;
        if (inFrame.find("evo.pdb") != std::string::npos)
        {
            count++;
            PDBInterface* theModelPDB = new PDBInterface(inFrame);
            ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
            molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
            protein* model = static_cast<protein*>(modelMol);
            //vector < double> Energy = model->chainBindingEnergy();
            fstream fs;
            fs.open ("seqpool.out", fstream::in | fstream::out | fstream::app);
            for (UInt i = 0; i < model->getNumResidues(0); i++)
            {
                 fs << model->getTypeFromResNum(0,i) << ",";
            }
            fs << endl;
            //Energy.clear();
            fs.close();
            delete theModelPDB;
        }
    }
    closedir(pdir);

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
vector < vector < UInt > > buildSequencePool()
{
    ifstream file("finalsequences.out");
    string item, line;
    vector < UInt > sequence;
    vector < vector < UInt > > sequencePool;
    while(getline(file,line))
    {
        stringstream stream(line);
        while(getline(stream,item,','))
        {
            stringstream aaString(item);
            int aaIndex;
            aaString >> aaIndex;
            sequence.push_back(aaIndex);
        }
        sequencePool.push_back(sequence);
        sequence.clear();
    }
    file.close();
    return sequencePool;
}
