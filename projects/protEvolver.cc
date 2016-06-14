//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************     protEvolver 3.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//********   -Fold Selective Protein Genetic Algorithm Based Evolution in Implicit Solvent -   **********
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

vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues);
UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UInt *_activeResidues);
void createPossibleMutantsDatabase(protein* bundle, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, bool _homosymmetric);
vector < vector < UInt > > buildSequencePool();
vector < vector < UInt > > buildPossibleMutants();
enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf"};

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    //--Running parameters
    if (argc !=2)
    {
        cout << "protEvolver <inFile.pdb>" << endl;
        exit(1);
    }

    //-- user inputs for evolution run
    UInt _activeChains[] = {0};                                                                 // chains active for mutation
    UInt _allowedLResidues[] = {A,R,Q,E,I,L,K,M,F,W,Y,V};                                       // amino acids allowed with phi < 0
    UInt _allowedDResidues[] = {G};                                                             // amino acids allowed with phi > 0
    UInt _activeResidues[] = {1,2,4,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27};  // positions active for mutation
    UInt _randomResidues[] = {1,2,4,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27};  // positions active for a random start sequence initially
    UInt _frozenResidues[] = {3,5,6,13};                                                        // positions that cannot move at all
    bool homoSymmetric = true;                                                                  // if true all chains are structurally symmetrical to the one listed active chain above
    bool backboneRelaxation = false;                                                            // if true allow minor backbone relaxation in structural optimization

    //--running parameters
    residue::setCutoffDistance(9.0);
    rotamer::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);

    //convert input arrays to vectors
    UInt activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]), randomResiduesSize = sizeof(_randomResidues)/sizeof(_randomResidues[0]), activeResiduesSize = sizeof(_activeResidues)/sizeof(_activeResidues[0]);
    UInt allowedLResiduesSize = sizeof(_allowedLResidues)/sizeof(_allowedLResidues[0]), allowedDResiduesSize = sizeof(_allowedDResidues)/sizeof(_allowedDResidues[0]), frozenResiduesSize = sizeof(_frozenResidues)/sizeof(_frozenResidues[0]);
    UIntVec activeChains, allowedLResidues, allowedDResidues, activeResidues, randomResidues, frozenResidues;
    for (UInt i = 0; i < activeChainsSize; i++)
    {
        activeChains.push_back(_activeChains[i]);
    }
    for (UInt i = 0; i < allowedLResiduesSize; i++)
    {
        allowedLResidues.push_back(_allowedLResidues[i]);
    }
    for (UInt i = 0; i < allowedDResiduesSize; i++)
    {
        allowedDResidues.push_back(_allowedDResidues[i]);
    }
    for (UInt i = 0; i < activeResiduesSize; i++)
    {
        activeResidues.push_back(_activeResidues[i]);
    }
    for (UInt i = 0; i < randomResiduesSize; i++)
    {
        randomResidues.push_back(_randomResidues[i]);
    }
    for (UInt i = 0; i < frozenResiduesSize; i++)
    {
        frozenResidues.push_back(_frozenResidues[i]);
    }

    //--set initial variables
    srand (getpid());
    double phi, bestEnergy, pastEnergy, Energy, randStartE;
    UInt timeid, sec, mutant = 0, numResidues, plateau = activeResiduesSize, nobetter = 0;
    vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
    vector < vector < UInt > > sequencePool, proteinSequence, finalSequence, possibleMutants;
    vector < double > bindingEnergy;
    stringstream convert;
    string infile = argv[1];
    string startstr, outFile;
    UInt name = rand() % 100000000;
    convert << name, startstr = convert.str();
    string tempModel = startstr + "_temp.pdb";


    //--determine which allowed amino acids are possible for each position from activeResidues
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    possibleMutants = buildPossibleMutants();
    if(possibleMutants.size() < activeResidues.size())
    {
        createPossibleMutantsDatabase(bundle, activeChains, activeResidues, allowedLResidues, homoSymmetric);
        possibleMutants = buildPossibleMutants();
    }
    delete thePDB;

    //--Run multiple independent evolution cycles-----------------------------------------------------
    for (UInt a = 1; a < 500; a++)
    {
        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* bundle = static_cast<protein*>(pMol);
        sequencePool = buildSequencePool();
        if (homoSymmetric)
        {
            for (UInt i = 1; i < bundle->getNumChains(); i++)
            {
                bundle->symmetryLinkChainAtoB(i, activeChains[0]);
            }
        }

        //--load in initial pdb and mutate in random starting sequence on active chains and random residues
        nobetter = 0;
        for (UInt i = 0; i < activeChains.size(); i++)
        {
            for (UInt j = 0; j < randomResidues.size(); j++)
            {
                bundle->setMoved(activeChains[i],randomResidues[j],1);
                bundle->activateForRepacking(activeChains[i], randomResidues[j]);
                randomPosition.push_back(activeChains[i]);
                randomPosition.push_back(randomResidues[j]);
                phi = bundle->getPhi(activeChains[i], randomResidues[j]);
                if (phi > 0 && phi < 180)
                {
                    mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition, _activeResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
                }
                if (phi < 0 && phi > -180)
                {
                    mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition, _activeResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
                }
                randomPosition.clear();
            }
            chainSequence = getChainSequence(bundle, activeChains[i]);
            proteinSequence.push_back(chainSequence);
        }
        if (homoSymmetric)
        {
            bundle->protOpt(backboneRelaxation, frozenResidues, activeChains[0]);
        }
        else
        {
            bundle->protOpt(backboneRelaxation);
        }
        randStartE = bundle->protEnergy();

        //--Determine next mutation position
        mutantPosition.clear();
        mutantPosition = getMutationPosition(bundle, activeChains, activeResidues);
        pdbWriter(bundle, tempModel);

        //--set Energy startpoint
        if (homoSymmetric)
        {
            Energy = bundle->intraSoluteEnergy(true, activeChains[0]);
        }
        else
        {
            Energy = bundle->protEnergy();
        }
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
            if (homoSymmetric)
            {
                for (UInt i = 1; i < bundle->getNumChains(); i++)
                {
                    bundle->symmetryLinkChainAtoB(i, activeChains[0]);
                }
            }

            //--Mutate current sequence, new mutant and optimize system
            nobetter++;
            for (UInt i = 0; i < activeChains.size(); i++)
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
                            mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition, _activeResidues);
                            bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
                        }
                        if (phi < 0 && phi > -180)
                        {
                            mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition, _activeResidues);
                            bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
                        }
                    }
                    else if (j != frozenResidues[0] && j != frozenResidues[1])
                    {
                        bundle->mutateWBC(activeChains[i],j, proteinSequence[i][j]);
                    }
                }
            }
            if (homoSymmetric)
            {
                bundle->protOpt(backboneRelaxation, frozenResidues, activeChains[0]);
            }
            else
            {
                bundle->protOpt(backboneRelaxation);
            }
            protein* tempBundle = new protein(*bundle);

            //--Determine next mutation position
            mutantPosition.clear();
            mutantPosition = getMutationPosition(bundle, activeChains, activeResidues);

            //--Energy test
            if (homoSymmetric)
            {
                Energy = bundle->intraSoluteEnergy(true, activeChains[0]);
            }
            else
            {
                Energy = bundle->protEnergy();
            }
            if (Energy < pastEnergy)
            {
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
        bindingEnergy.clear();
        bindingEnergy = model->chainBindingEnergy();
        if (bindingEnergy[0] < randStartE && bindingEnergy[0] < 0)
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
            for (UInt i = 0; i < activeChains.size(); i++)
            {
                chainSequence = getChainSequence(model, activeChains[i]);
                finalSequence.push_back(chainSequence);
            }
            fstream finalline;
            finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
            finalline << timeid << " " << bindingEnergy[0] << " " << bindingEnergy[1] << " ";

            fstream fs;
            fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
            for (UInt i = 0; i < activeChains.size(); i++)
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
            fs.close();
        }
        delete theModelPDB;
        bindingEnergy.clear(),sequencePool.clear(),proteinSequence.clear(), chainSequence.clear(), mutantPosition.clear(), chainSequence.clear(), sequencePosition.clear(), randomPosition.clear();
        bindingEnergy.resize(0),sequencePool.resize(0),proteinSequence.resize(0), chainSequence.resize(0), mutantPosition.resize(0), chainSequence.resize(0), sequencePosition.resize(0), randomPosition.resize(0);
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

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

vector <UInt> getMutationPosition(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues)
{
    //--get median residue energy
    UInt randres, randchain;
    double posE, medE;
    vector <UInt> _mutantPosition;
    medE = _prot->getMedianResEnergy(_activeChains, _activeResidues);

    //--find random position with worse than median energy
    do
    {
        randchain = _activeChains[rand() % _activeChains.size()];
        randres = _activeResidues[rand() % _activeResidues.size()];
        posE = _prot->resEnergy(randchain,randres);
    }while (posE < medE);
    _mutantPosition.push_back(randchain);
    _mutantPosition.push_back(randres);
    return _mutantPosition;
}

UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UInt *_activeResidues)
{
    float mutant, chance, entropy, acceptance, pooling, resFreqAccept;
    vector <UInt> resFreqs(58,1);
    UInt position;
    float count = _sequencePool.size();

    //--get sequence evolution results for position
    for (UInt i = 0; i < _sequencePool.size(); i++)
    {
        UInt restype = _sequencePool[i][_mutantPosition[1]];
        resFreqs[restype] = resFreqs[restype] + 1;
    }

    //--find mutant in vector
    for (UInt i = 0; i < _possibleMutants.size(); i++)
    {
        if (_activeResidues[i] == _mutantPosition[1])
        {
            position = i;
        }
    }
    UInt positionPossibles = _possibleMutants[position].size();

    //--determine population based chance of mutation acceptance or a random mutation, via linear regression of sequence entropy
    do
    {
        chance = (rand() % 100) + 1;
        entropy = (rand() % 100) + 1; //sequence entropy determined by pooling linear decline to resolve minima after suitable diversity
        mutant = _possibleMutants[position][rand() % positionPossibles];
        pooling = (-0.316 * count) + 285; //600 sequences equals 5% chance of pooling sequences, 100% at >= 900
        if (entropy > pooling)
        {
            resFreqAccept = resFreqs[mutant];
            acceptance = (resFreqAccept/(count-1))*100; //chance of accepting given amino acid at position is proportional to population
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
    ifstream file("sequencepool.out");
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

vector < vector < UInt > > buildPossibleMutants()
{
    ifstream file("possiblemutants.out");
    string item, line;
    vector < UInt > _position;
    vector < vector < UInt > > _possibleMutants;
    while(getline(file,line))
    {
        stringstream stream(line);
        while(getline(stream,item,','))
        {
            stringstream aaString(item);
            int aaIndex;
            aaString >> aaIndex;
            _position.push_back(aaIndex);
        }
        _possibleMutants.push_back(_position);
        _position.clear();
    }
    file.close();
    return _possibleMutants;
}

void createPossibleMutantsDatabase(protein* _bundle, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, bool _homoSymmetric)
{
    double Energy;
    fstream pm;
    pm.open ("possiblemutants.out", fstream::in | fstream::out | fstream::app);
    if (_homoSymmetric)
    {
        for (UInt i = 1; i < _bundle->getNumChains(); i++)
        {
            _bundle->symmetryLinkChainAtoB(i, _activeChains[0]);
        }
    }
    for (UInt i = 0; i < _activeChains.size(); i++)
    {
        for (UInt j = 0; j <_activeResidues.size(); j++)
        {
            for (UInt k = 0; k <_allowedLResidues.size(); k++)
            {
                _bundle->activateForRepacking(_activeChains[i], _activeResidues[j]);
                _bundle->mutateWBC(_activeChains[i], _activeResidues[j], _allowedLResidues[k]);
                UIntVec allowedRots = _bundle->getAllowedRotamers(_activeChains[i], _activeResidues[j], _allowedLResidues[k], 0);
                if (allowedRots.size() > 0)
                {
                    for (UInt l = 0; l < allowedRots.size(); l++)
                    {
                        _bundle->setRotamerWBC(_activeChains[i], _activeResidues[j], 0, allowedRots[l]);
                        _bundle->setMoved(_activeChains[i], _activeResidues[j], 1);
                        if (_homoSymmetric)
                        {
                            Energy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
                        }
                        else
                        {
                            Energy = _bundle->protEnergy();
                        }
                        if (Energy < 0)
                        {
                            pm << _allowedLResidues[k] << ",";
                            break;
                        }
                    }
                }
                else
                {
                    if (_homoSymmetric)
                    {
                        Energy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
                    }
                    else
                    {
                        Energy = _bundle->protEnergy();
                    }
                    if (Energy < 0)
                    {
                        pm << _allowedLResidues[k] << ",";
                    }
                }
            }
            _bundle->mutateWBC(_activeChains[i], _activeResidues[j], A);
            pm << endl;
        }
    }
}
