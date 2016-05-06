//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************     protEvolver 3.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//*******   -Folding Selective Protein Genetic Algorithm Based Evolution in Implicit Solvent -   ********
//*******************************************************************************************************

/////// Just specify infile structure, active chains and residues indexes, and it will evolve a sequence favorable for folding

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

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UInt* _activeChains, UInt _activeChainsSize, UInt* _activeResidues, UInt _activeResiduesSize);
UInt getProbabilisticMutation(vector < vector < UInt > > _sequencePool, vector <UInt> _mutantPosition, UInt *_aminoacids, UInt aaSize);
vector < vector < UInt > > buildSequencePool();

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    //--Running parameters
    if (argc !=3)
    {
        cout << "protEvolver <inFile.pdb>" << endl;
        exit(1);
    }
    string infile = argv[1];
    string localout = argv[2];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf"};
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
    srand (getpid());

    //--inputs for mutation
    UInt activeChains[] = {1};
    UInt allowedLResidues[] = {A,R,N,D,Q,E,I,L,K,M,F,S,T,W,Y,V};
    UInt activeResidues[] = {1,2,3,4,5,6,7,8,9,10};
    UInt randomResidues[] = {1,2,3,4,5,6,7,8,9,10};
    UInt allowedDResidues[] = {G};

    double phi, bestEnergy, pastEnergy, Energy;
    UInt nobetter = 0, activeChainsSize = sizeof(activeChains)/sizeof(activeChains[0]), randomResiduesSize = sizeof(randomResidues)/sizeof(randomResidues[0]), activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]);
    UInt lResidues = sizeof(allowedLResidues)/sizeof(allowedLResidues[0]), dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]);
    UInt name, timeid, sec, mutant = 0, numResidues, plateau = activeResiduesSize;
    vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
    vector < vector < UInt > > proteinSequence, finalSequence, sequencePool;
    stringstream convert;
    string startstr, outFile;
    name = rand() % 100000000;
    convert << name, startstr = convert.str();
    string tempModel = startstr + "_temp.pdb";
    delete thePDB;

    //--Run multiple independent evolution cycles-----------------------------------------------------
    for (UInt a = 1; a < 500; a++)
    {
        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* bundle = static_cast<protein*>(pMol);
        sequencePool = buildSequencePool();

        //--load in initial pdb and mutate in random starting sequence on active chains and random residues
        proteinSequence.clear(), chainSequence.clear(), nobetter = 0;
        for (UInt i = 0; i < activeChainsSize; i++)
        {
            for (UInt j = 0; j < randomResiduesSize; j++)
            {
                bundle->activateForRepacking(activeChains[i], randomResidues[j]);
                randomPosition.push_back(activeChains[i]);
                randomPosition.push_back(randomResidues[j]);
                phi = bundle->getPhi(activeChains[i], randomResidues[j]);
                if (phi > 0 && phi < 180)
                {
                    mutant = getProbabilisticMutation(sequencePool, randomPosition, allowedDResidues, dResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
                }
                if (phi < 0 && phi > -180)
                {
                    mutant = getProbabilisticMutation(sequencePool, randomPosition, allowedLResidues, lResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
                }
                randomPosition.clear();
            }
            chainSequence = getChainSequence(bundle, activeChains[i]);
            proteinSequence.push_back(chainSequence);
        }
        bundle->protOpt(false);

        //--Determine next mutation position
        mutantPosition.clear();
        mutantPosition = getMutationPosition(bundle, activeChains, activeChainsSize, activeResidues, activeResiduesSize);
        pdbWriter(bundle, tempModel);

        //--set Energy startpoint
        Energy = bundle->protEnergy();
        pastEnergy = Energy;
        bestEnergy = Energy;
        delete thePDB;

        //--Run through a single evolutionary path (ancestral line) till hitting plateau
        do
        {
            PDBInterface* thePDB = new PDBInterface(infile);
            ensemble* theEnsemble = thePDB->getEnsemblePointer();
            molecule* pMol = theEnsemble->getMoleculePointer(0);
            protein* bundle = static_cast<protein*>(pMol);

            //--Mutate current sequence, new mutant and optimize system
            nobetter++;
            for (UInt i = 0; i < activeChainsSize; i++)
            {
                numResidues = bundle->getNumResidues(activeChains[i]);
                for (UInt j = 0; j < numResidues; j++)
                {
                    bundle->setMoved(activeChains[i], j,1);
                    bundle->activateForRepacking(activeChains[i],j);
                    if (activeChains[i] == mutantPosition[0] && j == mutantPosition[1])
                    {
                        //--new mutant
                        sequencePosition.push_back(i);
                        sequencePosition.push_back(j);
                        phi = bundle->getPhi(mutantPosition[0], mutantPosition[1]);
                        if (phi > 0 && phi < 180)
                        {
                            mutant = getProbabilisticMutation(sequencePool, mutantPosition, allowedDResidues, dResidues);
                            bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
                        }
                        if (phi < 0 && phi > -180)
                        {
                            mutant = getProbabilisticMutation(sequencePool, mutantPosition, allowedLResidues, lResidues);
                            bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
                        }
                    }
                    else
                    {
                        bundle->mutateWBC(activeChains[i],j, proteinSequence[i][j]);
                    }
                }
                //randomizeSideChains(bundle, activeChains[i]);
            }
            bundle->protOpt(false);
            protein* tempBundle = new protein(*bundle);

            //--Determine next mutation position
            mutantPosition.clear();
            mutantPosition = getMutationPosition(bundle, activeChains, activeChainsSize, activeResidues, activeResiduesSize);

            //--Energy test
            Energy = bundle->protEnergy();
            if (Energy < pastEnergy)
            {
                //cout << Energy << " " << nobetter << endl;
                if (Energy < bestEnergy)
                {
                    bestEnergy = Energy;
                    pdbWriter(tempBundle, tempModel);
                }
                proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, pastEnergy = Energy;
                if (nobetter > 0) { nobetter--;
                }
                else{ nobetter = 0;
                }
            }
            sequencePosition.clear();
            delete thePDB;
            delete tempBundle;
        }while (nobetter < plateau);

        //--Print final energy and write a pdb file----------------------------------------------------
        PDBInterface* theModelPDB = new PDBInterface(tempModel);
        ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
        molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
        protein* model = static_cast<protein*>(modelMol);
        vector <double> bindingEnergy = model->chainBindingEnergy();
        //UInt numchains = model->getNumChains();
        if (bindingEnergy[0] < 0)
        {
            name = rand() % 100;
            sec = time(NULL);
            timeid = name + sec;
            stringstream convert;
            string countstr;
            convert << timeid, countstr = convert.str();
            outFile = countstr + ".evo.pdb";
            pdbWriter(model, outFile);
            finalSequence.clear(), chainSequence.clear();
            for (UInt i = 0; i < activeChainsSize; i++)
            {
                chainSequence = getChainSequence(model, activeChains[i]);
                finalSequence.push_back(chainSequence);
            }
            fstream finalline;
            finalline.open ("final.out", fstream::in | fstream::out | fstream::app);
            finalline << timeid << " " << bindingEnergy[0] << " " << bindingEnergy[1] << " ";

        fstream finallocal;
        finallocal.open (localout, fstream::in | fstream::out | fstream::app);
            finallocal << timeid << " " << bindingEnergy[0] << " " << bindingEnergy[1] << " ";

            fstream fs;
            fs.open ("finalsequences.out", fstream::in | fstream::out | fstream::app);
            for (UInt i = 0; i < activeChainsSize; i++)
            {
                for (UInt j = 0; j < finalSequence[i].size(); j++)
                {
                    finalline << aminoAcidString[finalSequence[i][j]] << " ";
                    fs << finalSequence[i][j] << ",";
                }
            }
            fs << endl;
            finalline << endl;
            finalline.close();
            finallocal << endl;
            finallocal.close();
            fs.close();
        }
        delete theModelPDB;
        sequencePool.clear();
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void randomizeSideChains(protein* _prot, UInt _chainIndex)
{
    UInt allowedRotsSize, randrot, restype, numResidues;
    UIntVec allowedRots;
    numResidues = _prot->getNumResidues(_chainIndex);
    for (UInt j = 0; j < numResidues; j++)
    {
        restype = _prot->getTypeFromResNum(_chainIndex, j);
        allowedRots = _prot->getAllowedRotamers(_chainIndex, j, restype, 0);
        allowedRotsSize = allowedRots.size();
        if (allowedRotsSize > 2)
        {
            randrot = rand() % allowedRotsSize;
            _prot->setRotamerWBC(_chainIndex, j, 0, allowedRots[randrot]);
        }
    }
    return;
}

vector < UInt > getChainSequence(protein* _prot, UInt _chainIndex)
{
    UInt restype, numResidues;
    vector < UInt > sequence;
    numResidues = _prot->getNumResidues(_chainIndex);
    for (UInt j = 0; j < numResidues; j++)
    {
        restype = _prot->getTypeFromResNum(_chainIndex, j);
        sequence.push_back(restype);
    }
    return sequence;
}

vector <UInt> getMutationPosition(protein* _prot, UInt *_activeChains, UInt _activeChainsSize, UInt *_activeResidues, UInt _activeResiduesSize)
{
    //--get median residue energy
    UInt randres, randchain;
    double posE, medE;
    vector <UInt> _mutantPosition;
    medE = _prot->getMedianResEnergy();

    //--find random position with worse than median energy
    do
    {
        randchain = _activeChains[rand() % _activeChainsSize];
        randres = _activeResidues[rand() % _activeResiduesSize];
        posE = _prot->resEnergy(randchain,randres);
    }while (posE < medE);
    _mutantPosition.push_back(randchain);
    _mutantPosition.push_back(randres);
    return _mutantPosition;
}

UInt getProbabilisticMutation(vector < vector < UInt > > _sequencePool, vector <UInt> _mutantPosition, UInt *_aminoacids, UInt aaSize)
{
    UInt mutant, chance, entropy;
    double acceptance;
    vector <UInt> resFreqs(55,1);
    UInt count = _sequencePool.size();

    //--get sequence evolution results for position
    for (UInt i = 0; i < _sequencePool.size(); i++)
    {
        UInt restype = _sequencePool[i][_mutantPosition[1]];
        resFreqs[restype] = resFreqs[restype] + 1;
    }

    //--determine population based chance of mutation acceptance or a random mutation, via linear regression of sequence entropy
    do
    {
        chance = rand() % 100;
        entropy = rand() % 100; //sequence entropy determined by pooling linear decline to resolve minima after suitable diversity
        mutant = _aminoacids[rand() % aaSize];
        UInt pooling = -0.09 * count + 140; //500 sequences equals 5% chance of pooling sequences, 95% chance of it being random, 1000=50% ...
        if (entropy > pooling)
        {
            acceptance = (resFreqs[mutant]/count)*100; //chance of accepting given amino acid at position is proportional to population
        }
        else
        {
            acceptance = 100;  //random mutation
        }
    }while (chance > acceptance);
    return mutant;
}

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
