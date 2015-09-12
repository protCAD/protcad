//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  protNetwork 1.0  *************************************************
//***********************************                    ************************************************
//*******************************************************************************************************
//**********   -point mutation network build, then backbone and sidechain optimization-   ***************
//*******************************************************************************************************

/////// Just specify a infile and preferred outfile name.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

void protOptNetwork(protein* _prot, UInt _plateau, bool _backbone, int *positions);

int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=3)
	{
        cout << "protNetwork <inFile.pdb> <outfile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    string outfile = argv[2];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Hce};
    //string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Hce"};
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
    srand (time(NULL));
	
    int set1[] = {46,47,61,92,95};//{1,3,21,30,54,59,61,72,78,90,92,110,112,116,129,133,166,174,181};
    cout << bundle->intraSoluteEnergy(true) << endl;
    protOptNetwork(bundle, 500, false, set1);
    cout << bundle->intraSoluteEnergy(true) << endl;
    pdbWriter(bundle, outfile);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void protOptNetwork(protein* _prot, UInt _plateau, bool _backbone, int* positions)
{	// Sidechain and backbone optimization with a polarization based dielectric scaling of electrostatics and corresponding implicit solvation score
    //    _plateau: the number of consecutive optimization cycles without an energy decrease.
    //	    	     (250 is recommended for a full minimization without excessive calculation)
    // -pike 2013

    //--Initialize variables for loop and calculate starting energy-------------------------------------
    bool dielectrics = true;
    if (residueTemplate::itsAmberElec.getScaleFactor() == 0.0)
    {
        dielectrics = false;
    }
    double deltaTheta = 0, totalpreposE = 0, avepreposE = -1E10, randangle;
    double Energy, preposE, currentposE, pastEnergy = _prot->intraSoluteEnergy(dielectrics);
    UInt randchain, randres, randrestype, allowedRotsize, randrot, number = 0, nobetter = 0;
    UInt randtype, chainNum = _prot->getNumChains(), rotbetter = 0;
    vector < vector <double> > currentRot;
    int thisone;
    int positionSize = sizeof(positions)/sizeof(positions[0]);
    UIntVec allowedRots;
    srand (time(NULL));

    //--Run optimizaiton loop to energetic minima, determined by _plateau-------------------------------
    do
    {
        //--Generate random residue
        randchain = rand() % chainNum;
        randres = positions[rand() % positionSize];
        randrestype = _prot->getTypeFromResNum(randchain, randres);
        preposE = _prot->getPositionSoluteEnergy(randchain, randres, dielectrics);
        if (randrestype == 0 || randrestype == 19 || randrestype == 20 || randrestype == 26 || randrestype == 27 || randrestype == 46 || randrestype == 47)
        {nobetter++;
        }
        else
        {nobetter++,nobetter++;
        }

        //--backbone optimization----------------------------------------------------------------------
        if (rotbetter > _plateau && preposE > avepreposE && _backbone)
        {
            //--choose phi or psi and angle, for a local transformation
            randtype = rand() % 2;
            do
            { deltaTheta = ((rand() % 3) -1);
            } while (deltaTheta == 0);

            //--transform angles while energy improves, until energy degrades, then revert one step
            do
            {
                _prot->setDihedralLocal(randchain, randres, deltaTheta, randtype);
                currentposE = _prot->getPositionSoluteEnergy(randchain, randres, dielectrics), thisone = 0;
                //--Energy test
                if (currentposE < (preposE - .05))
                {
                    Energy = _prot->intraSoluteEnergy(dielectrics);
                    if (Energy < pastEnergy)
                    {
                        //cout << Energy << endl;
                        nobetter = 0, thisone = 1, pastEnergy = Energy, preposE = currentposE;
                    }
                }
            } while (thisone == 1);
            _prot->setDihedralLocal(randchain, randres, (deltaTheta*-1), randtype);
        }

        //--Rotamer optimization-----------------------------------------------------------------------
        if (preposE > avepreposE)
        {
            //--Get current rotamer and allowed
            currentRot = _prot->getSidechainDihedrals(randchain, randres);
            allowedRots = _prot->getAllowedRotamers(randchain, randres, randrestype, 0);
            allowedRotsize = (allowedRots.size() * 0.33), rotbetter++, rotbetter++;

            //--Try 1/3 of allowed rotamers keep first improvement or revert to previous angles
            for (UInt j = 0; j < allowedRotsize; j ++)
            {
                randrot = rand() % allowedRots.size();
                _prot->setRotamerWBC(randchain, randres, 0, allowedRots[randrot]);
                randangle = 0;
                for (UInt c = 0; c < 12; c++)
                {
                    randangle = randangle+30;
                    _prot->setChi(randchain, randres, 1, 0, randangle);
                    Energy = _prot->intraSoluteEnergy(dielectrics);
                    if (Energy < pastEnergy)
                    {
                        //cout << Energy << endl;
                        rotbetter--, rotbetter--, nobetter = 0, pastEnergy = Energy, preposE = currentposE;
                        break;
                    }
                    if (nobetter != 0)
                    {
                        _prot->setSidechainDihedralAngles(randchain, randres, currentRot);
                        _prot->setChi(randchain, randres, 1, 0, randangle*-1);
                    }
                }
            }
        }

        //--check status of optimization---------------------------------------------------------------
        if (number == _plateau)
        {
            number = 0, totalpreposE = 0;
        }
        number++,number++, totalpreposE = (totalpreposE + preposE), avepreposE = (totalpreposE/number);
    } while (nobetter < _plateau * 1.2);
    return;
}
