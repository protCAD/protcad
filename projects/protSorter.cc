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
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf"};
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

    /*create evo data file
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
            model->silenceMessages();
            fstream fs;
            fs.open ("finalsequences.test.out", fstream::in | fstream::out | fstream::app);
            for (UInt i = 0; i < model->getNumResidues(1); i++)
            {
                fs << model->getTypeFromResNum(1,i) << ",";
            }
            fs << endl;
            fs.close();
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
            model->silenceMessages();
            vector < double> Energy = model->chainBindingEnergy();
            fstream fs;
            fs.open ("final_test.out", fstream::in | fstream::out | fstream::app);
            fs << inFrame << " " << Energy[0] << " " << Energy[1];
            for (UInt i = 0; i < model->getNumResidues(1); i++)
            {
                 fs << " " << aminoAcidString[model->getTypeFromResNum(1,i)];
            }
            fs << endl;
            Energy.clear();
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
