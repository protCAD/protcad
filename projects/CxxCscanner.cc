#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

using namespace std;

bool checkAlaSF4(UInt _AlaID, UInt _CLUID, protein* _prot);
void checkSulfurGeometry (protein* _prot, UInt _CLUID, UInt _C1, UInt _C2);

int main (int argc, char* argv[])
{
	enum aminoAcid {A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,ALD,AMD,CLU};

	string inFile = argv[1];
	double step; sscanf(argv[2], "%lf", &step);
	UInt CLUrot; sscanf(argv[3], "%u", &CLUrot);

	// 
	PDBInterface* thePDB = new PDBInterface(inFile);
 	ensemble* theEnsemble = thePDB->getEnsemblePointer();
 	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	prot->silenceMessages();

	residue::setCutoffDistance(5.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();


	// assume CXX-CLU-XXC
	UInt C1 = 0; UInt C2 = 6;
	UInt CLUID = 3;
	int sequence[] = {A,G,G,CLU,G,G,A}; 
	UInt seqLength = (UInt)sizeof(sequence)/sizeof(sequence[0]);
	for (UInt i = 0; i < seqLength; i ++)
	{
		prot->activateForRepacking(0,i);
		prot->mutate(0,i,sequence[i]);
	}

	// check if given rotamer matches allowed list  
	UIntVec allowedCLUrots = prot->getAllowedRotamers(0,CLUID,CLU,0);
	bool rotAllowed = false;
	for (UInt i = 0; i < allowedCLUrots.size(); i ++)
	{
		if (allowedCLUrots[i] == CLUrot) rotAllowed = true;
	}
	if (!rotAllowed)
	{
		cout << "rotamer " << CLUrot << " not allowed." << endl;
		exit(1);
	}
	prot->setRotamer(0,CLUID,0, CLUrot);

	// get energetically allowed conformers
	vector < UIntVec > Nconformer;  // format psi0, phi1, psi1 ... phi5, psi5, phi6
	for (UInt i = 0; i < prot->getNumResidues(0); i ++) // make initial extended conformation
	{
		if (i != 0) prot->setPhi(0,i,180.0);
		if (i != prot->getNumResidues(0)-1) prot->setPsi(0,i,180.0);
	}
	double extendedE = prot->intraEnergy();

	double thresholdE = 100.0;

	double startTorsion = -180.0;
	for (double psi0 = startTorsion; psi0 < 180.0; psi0 += step)
	{ prot->setPsi(0,0,psi0);
	  for (double phi1 = startTorsion; phi1 < 180.0; phi1 += step)
	  { while (phi1 >= -30.0 and phi1 <= 30.0) phi1 += step; prot->setPhi(0,1,phi1);
	    for (double psi1 = startTorsion; psi1 < 180.0; psi1 += step)
	    { prot->setPsi(0,1,psi1);
	      for (double phi2 = startTorsion; phi2 < 180.0; phi2 += step)
	      { while (phi2 >= -30.0 and phi2 <= 30.0) phi2 += step; prot->setPhi(0,2,phi2);
	        for (double psi2 = startTorsion; psi2 < 180.0; psi2 += step)
	        { prot->setPsi(0,2,psi2);
	          for (double phi3 = startTorsion; phi3 < 180.0; phi3 += step)
		  { 	while (phi3 >= -30.0 and phi3 <= 30.0) phi3 += step; prot->setPhi(0,3,phi3);
			if (checkAlaSF4(C1,CLUID,prot))			
			{
				double E = prot->intraEnergy();
				if ((E-extendedE) <= thresholdE)
				{
					UIntVec temp;
					temp.push_back(psi0);
					temp.push_back(phi1);
					temp.push_back(psi1);
					temp.push_back(phi2);
					temp.push_back(psi2);
					temp.push_back(phi3);
					Nconformer.push_back(temp);
					temp.resize(0);
				}
			}
	} } } } } }
	cout << "N " << Nconformer.size()  << endl;


	vector < UIntVec > Cconformer;
	for (UInt i = 0; i < prot->getNumResidues(0); i ++) // make initial extended conformation
	{
		if (i != 0) prot->setPhi(0,i,180.0);
		if (i != prot->getNumResidues(0)-1) prot->setPsi(0,i,180.0);
	}
	for (double psi3 = startTorsion; psi3 < 180.0; psi3 += step)
	{ prot->setPsi(0,3,psi3);
	  for (double phi4 = startTorsion; phi4 < 180.0; phi4 += step)
	  { while (phi4 >= -30.0 and phi4 <= 30.0) phi4 += step;  prot->setPhi(0,4,phi4);
	    for (double psi4 = startTorsion; psi4 < 180.0; psi4 += step)
	    { prot->setPsi(0,4,psi4);
	      for (double phi5 = startTorsion; phi5 < 180.0; phi5 += step)
	      { while (phi5 >= -30.0 and phi5 <= 30.0) phi5 += step;  prot->setPhi(0,5,phi5);
	        for (double psi5 = startTorsion; psi5 < 180.0; psi5 += step)
	        { prot->setPsi(0,5,psi5);
	          for (double phi6 = startTorsion; phi6 < 180.0; phi6 += step)
	          { 	while (phi6 >= -30.0 and phi6 <= 30.0) phi6 += step;  prot->setPhi(0,6,phi6);
			if (checkAlaSF4(C2,CLUID,prot))
			{
				double E = prot->intraEnergy();
				if ((E-extendedE) <= thresholdE)
				{
					UIntVec temp;
					temp.push_back(psi3);
					temp.push_back(phi4);
					temp.push_back(psi4);
					temp.push_back(phi5);
					temp.push_back(psi5);
					temp.push_back(phi6);
					Cconformer.push_back(temp);
					temp.resize(0);
		//			cout << "EC " << E << endl;
				}
			}
	} } } } } }	
	cout << "C " << Cconformer.size() << endl;

	vector < UIntVec > totalConformer;



	int totalCount = Nconformer.size() * Cconformer.size();
	int sum = 0;
	for (UInt i = 0; i < Nconformer.size(); i ++)
	{

		prot->setPsi(0,0,Nconformer[i][0]);
		prot->setPhi(0,1,Nconformer[i][1]);
		prot->setPsi(0,1,Nconformer[i][2]);
		prot->setPhi(0,2,Nconformer[i][3]);
		prot->setPsi(0,2,Nconformer[i][4]);
		prot->setPhi(0,3,Nconformer[i][5]);
		for (UInt j = 0; j < Cconformer.size(); j ++)
		{		   	
			prot->setPsi(0,3,Cconformer[j][0]);
			prot->setPhi(0,4,Cconformer[j][1]);
			prot->setPsi(0,4,Cconformer[j][2]);
			prot->setPhi(0,5,Cconformer[j][3]);
			prot->setPsi(0,5,Cconformer[j][4]);
			prot->setPhi(0,6,Cconformer[j][5]);

			double E = prot->intraEnergy();
			if ((E-extendedE) <=thresholdE)
			{
				UIntVec temp;
				for (UInt k = 0; k < Nconformer[i].size(); k ++) temp.push_back(Nconformer[i][k]);
				for (UInt k = 0; k < Cconformer[j].size(); k ++) temp.push_back(Cconformer[j][k]);
				totalConformer.push_back(temp);
				temp.resize(0);
			}
			sum ++;
		}
		float freq = (float)sum/(float)totalCount;
		cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << freq * 100  << "%" ;
	}
	cout << endl;

	prot->mutate(0,C1,C); 
	prot->mutate(0,C2,C);

	pdbWriter(prot, "out.pdb");
	
	cout << "N " << Nconformer.size() << " C " << Cconformer.size() << " total " << totalConformer.size() << endl;
	cout << "checking sulfur geometry for conformer" << endl;

	for (UInt i = 0; i < totalConformer.size(); i ++)
	{
		/*-----------------------------------*/ prot->setPsi(0,0,totalConformer[i][0]);
		prot->setPhi(0,1,totalConformer[i][1]); prot->setPsi(0,1,totalConformer[i][2]);
		prot->setPhi(0,2,totalConformer[i][3]); prot->setPsi(0,2,totalConformer[i][4]);
		prot->setPhi(0,3,totalConformer[i][5]); prot->setPsi(0,3,totalConformer[i][6]);
		prot->setPhi(0,4,totalConformer[i][7]); prot->setPsi(0,4,totalConformer[i][8]);
		prot->setPhi(0,5,totalConformer[i][9]); prot->setPsi(0,5,totalConformer[i][10]);
		prot->setPhi(0,6,totalConformer[i][11]); /*----------------------------------*/


		checkSulfurGeometry(prot,CLUID,C1,C2);
	}

	return 0;
}

	
bool checkAlaSF4(UInt _AlaID, UInt _CLUID, protein* _prot)
{
	dblVec alaCB = _prot->getCoords(0,_AlaID, "CB");
	dblVec Fe1 = _prot->getCoords(0, _CLUID, "FE1");
	dblVec Fe2 = _prot->getCoords(0, _CLUID, "FE2");
	dblVec Fe3 = _prot->getCoords(0, _CLUID, "FE3");
	dblVec Fe4 = _prot->getCoords(0, _CLUID, "FE4");
	
	double cutOff = 4.0;
	if (CMath::distance(alaCB,Fe1) < cutOff) return true;
	if (CMath::distance(alaCB,Fe2) < cutOff) return true;
	if (CMath::distance(alaCB,Fe3) < cutOff) return true;
	if (CMath::distance(alaCB,Fe4) < cutOff) return true;

	return false;
}



void checkSulfurGeometry (protein* _prot, UInt _CLUID, UInt _C1, UInt _C2)
{		
	double PI = 3.141592;	

	dblVec CB1 = _prot->getCoords(0,_C1,"CB");
        dblVec CB3 = _prot->getCoords(0,_C2,"CB");
	
	dblVec FE1 = _prot->getCoords(0,_CLUID,"FE1");
	dblVec FE3 = _prot->getCoords(0,_CLUID,"FE3");
	dblVec FE4 = _prot->getCoords(0,_CLUID,"FE4");
	dblVec CLUS2 = _prot->getCoords(0,_CLUID,"S2");

	for (UInt j = 0; j<=2; j++)
	{
		for (UInt k = 0; k <=2; k++)	
		{	

		_prot->setRotamer(0,_C1,0,j);
		_prot->setRotamer(0,_C2,0,k);

		//Geometrical Constraints for the binding site

		dblVec SG1 = _prot->getCoords(0,_C1,"SG");
		dblVec SG2 = _prot->getCoords(0,_CLUID,"SG");
		dblVec SG3 = _prot->getCoords(0,_C2,"SG");
	
		double SS12 = CMath::distance(SG1,SG2);
		double SS13 = CMath::distance(SG1,SG3);
		double SS23 = CMath::distance(SG2,SG3);
	
		//Determine pairwise interaction with iron-SG distance
		
		double FES1 = 0;
		double FES2 = 0;
		double CSF1 = 0;
		double CSF2 = 0;
		double SFS1 = 0;
		double SFS2 = 0;

		double FS11 = CMath::distance(FE1,SG1);
		double FS31 = CMath::distance(FE3,SG1);	
		double FS41 = CMath::distance(FE4,SG1);

		double FS13 = CMath::distance(FE1,SG3);
		double FS33 = CMath::distance(FE3,SG3);
		double FS43 = CMath::distance(FE4,SG3);

//First CYS		
		if (FS11 < FS31 && FS11 < FS41)
		{
		FES1 = FS11;
		CSF1 = CMath::angle(CB1,SG1,FE1);
		SFS1 = CMath::angle(SG1,FE1,CLUS2);
		}
		else if (FS31 < FS11 && FS31 < FS41)
		{
		FES1 = FS31;
		CSF1 = CMath::angle(CB1,SG1,FE3);
		SFS1 = CMath::angle(SG1,FE3,CLUS2);
		}
		else if (FS41 < FS11 && FS41 < FS31)
		{
		FES1 = FS41;
		CSF1 = CMath::angle(CB1,SG1,FE4);
		SFS1 = CMath::angle(SG1,FE4,CLUS2);
		}
//Second CYS

		if (FS13 < FS33 && FS13 < FS43)
		{
		FES2 = FS13;	
		CSF2 = CMath::angle(CB3,SG3,FE1);
		SFS2 = CMath::angle(SG3,FE1,CLUS2);
		}
		else if (FS33 < FS13 && FS33 < FS43)
		{
		FES2 = FS33;
		CSF2 = CMath::angle(CB3,SG3,FE3);
		SFS2 = CMath::angle(SG3,FE3,CLUS2);
		}
		else if (FS43 < FS13 && FS43 < 33)
		{
		FES2 = FS43;
		CSF2 = CMath::angle(CB3,SG3,FE4);
		SFS2 = CMath::angle(SG3,FE4,CLUS2);
		}

		//Distance and energy condition starts here

		UInt FESmax = 3.0;
		UInt FESmin = 2.0;
		UInt SSmax = 7.0;
		UInt SSmin = 5.8;
		if (SS12 < SSmax && SS13 < SSmax && SS23 < SSmax && SS12 > SSmin && SS13 > SSmin && SS23 > SSmin)
		{
		if (FES1 < FESmax && FES1 > FESmin && FES2 < FESmax && FES2 > FESmin) // When Fe-S bond is FESmax to FESmin
		{
		double intraenergy = _prot->intraEnergy();		
		double FES1energy = 0.17*(pow((4.0/FES1),12)-2*pow((4.0/FES1),6));
		double FES2energy = 0.17*(pow((4.0/FES2),12)-2*pow((4.0/FES2),6));
		double bondenergy1 = 50*(5*pow((2.3/FES1),12)-6*pow((2.3/FES1),10))*cos((PI/180)*(CSF1-109.5))*cos((PI/180)*(SFS1-109.5));
		double bondenergy2 = 50*(5*pow((2.3/FES2),12)-6*pow((2.3/FES2),10))*cos((PI/180)*(CSF2-109.5))*cos((PI/180)*(SFS2-109.5));
		double newenergy = intraenergy-FES1energy-FES2energy;
		double bondenergy = bondenergy1+bondenergy2;
		double totalenergy = newenergy+bondenergy;

		if (newenergy < 0 && bondenergy < 0)
			{	
			UIntVec  myRotamer0 = _prot->getCurrentRotamer(0,_C1);	
			UIntVec  myRotamer3 = _prot->getCurrentRotamer(0,_CLUID);	
			UIntVec  myRotamer6 = _prot->getCurrentRotamer(0,_C2);	
			
			int psi0 = _prot->getPsi(0,0);
			int phi1 = _prot->getPhi(0,1); int psi1 = _prot->getPsi(0,1);
			int phi2 = _prot->getPhi(0,2); int psi2 = _prot->getPsi(0,2);
			int phi3 = _prot->getPhi(0,3); int psi3 = _prot->getPsi(0,3);
			int phi4 = _prot->getPhi(0,4); int psi4 = _prot->getPsi(0,4);
			int phi5 = _prot->getPhi(0,5); int psi5 = _prot->getPsi(0,5);
						       int phi6 = _prot->getPhi(0,6);
					
			stringstream out; out << totalenergy ;;
			string pdbid = "energy"+ out.str() + ".pdb";
			pdbWriter(_prot, pdbid);

			cout << myRotamer3[0] << "," << myRotamer0[0] << "," << myRotamer6[0] << "," 
			<< psi0 << "," << phi1 << "," << psi1 << "," << phi2 << "," << psi2 << "," << phi3 << "," 
			<< psi3 << "," << phi4 << "," << psi4 << "," << phi5 << "," << psi5 << "," << phi6 << "," 
			<< newenergy << "," << bondenergy << "," << totalenergy << "," 
			<< FES1 << "," << FES2 << "," << CSF1 << "," << CSF2 << "," << SFS1 << "," << SFS2 << endl;		
			}
		}}
		}//for loop
	} //for loop

return;		
}
