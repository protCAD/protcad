//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************      protSlipPlane    ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//********************************* -get sequence from pdb file- ****************************************
//*******************************************************************************************************

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
        cout << "protSlipPlane <inFile.pdb>" << endl;
		exit(1);
	}
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	_prot->silenceMessages();
	residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);

    double temp = 298.15;
    double EO = 8.854187817e-12; //Vacuum permit (F/m)
    double KbKcal = 0.0019872041; //boltzman (kcal/molK)
    double Kb = 1.3806485279e-23; //boltzman (j/K)
    double Na = 6.0221409e+23; //avogrados n
    double Ec = 1.60217662e-19; //electron charge
    double totEnergy = _prot->protEnergy();
    UInt numChains = _prot->getNumChains();
    cout << "Chain " << "Resdiue " << "Atom " << "Distance " << "DebyeLength " << endl;
    for (UInt i = 0; i < numChains; i++)
    {
        UInt numRes = _prot->getNumResidues(i);
        for (UInt j = 0; j < numRes; j++)
        {
            UInt numAtoms = _prot->getNumAtoms(i,j);
            for (UInt k = 0; k < numAtoms; k++)
            {
                for (UInt l = 3; l < 10000; l++)
                {
                    UInt atomType = k;
                    double atomDielectric = _prot->getDielectric(i,j,k);
                    double atomRadius = _prot->getRadius(i,j,k);
                    double charge = residueTemplate::itsAmberElec.getItsCharge(atomType, k);
                    double chargeSquared = charge*charge;
                    double waterDielectric = -0.3195 * (temp-274.15) + 86.115; //Malmberg and Maryott, 1956 JRNBS
                    double debyeLength=sqrt(EO*atomDielectric*Kb*temp/(Na*Ec*Ec*(0.1+0.1)))*1e10;
                    double surfacePotential=charge/(EO*atomDielectric*atomRadius*(1+atomRadius/debyeLength));
                    double potential=surfacePotential*atomRadius*exp((atomRadius-l)/debyeLength)/l;
                    double proteinSolventEnthalpy = -166 * (atomDielectric) * (chargeSquared/l);
                    if (proteinSolventEnthalpy < (-KbKcal*temp))
                    {
                        cout << proteinSolventEnthalpy <<  " " << -KbKcal*temp << " " << i << " " << j << " " << k << " " << l <<  endl;
                    }
                    /*if (abs(potential) > (8.6173324e-5*temp/Ec))
                    {
                        cout << "debye " << potential << " " << 8.6173324e-5*temp << i << " " << j << " " << k << " " << l << " " << debyeLength;
                        break;
                    }*/
                }
            }
        }
    }
	return 0;
}

