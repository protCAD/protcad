#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>

//TAKE ONE ARGUMENT: NAME OF ORIGINAL PDB FILE
//OUTPUT= OUTPUT.PDB

int main (int argc, char* argv[])
{
        enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
        //enum Atom_Index{N=1,CA,C,O,CB};
        string inputFileName = argv[1];

        PDBInterface* thePDB = new PDBInterface(inputFileName);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* theMol = theEnsemble->getMoleculePointer(0);
        protein* prot = static_cast<protein*>(theMol);

        residue::setCutoffDistance(8.0);
        pmf::setScaleFactor(0.0);
        rotamer::setScaleFactor(0.0);
        microEnvironment::setScaleFactor(0.0);
        amberVDW::setScaleFactor(1.0);
        amberVDW::setRadiusScaleFactor(0.9);
        amberVDW::setLinearRepulsionDampeningOff();
        amberElec::setScaleFactor(0.0);
        solvation::setItsScaleFactor(0.0);

        prot->activateAllForRepacking(0);
	
	//Bond lenghts and bond angles
	const double PI= 3.1415;

	vector<double>bond_length(4);
	bond_length[0]=1.449;bond_length[1]=1.522;
	bond_length[2]=1.229;bond_length[3]=1.383;

	vector<double>bond_angle(4);
	bond_angle[0]=121.9;bond_angle[1]=109.5;
	bond_angle[2]=120.5;bond_angle[3]=122.9;

	vector <double>bond_x(4);
	bond_x[3]=(180+bond_angle[1]-116.6-bond_angle[0])/2;
	bond_x[0]=(180-(bond_x[3]+bond_angle[0]));
	bond_x[1]=(180-(bond_x[0]+bond_angle[1]));
	bond_x[2]=(bond_angle[2]-bond_x[1])+180;
	
	cout<<"RESIDUE NUMBER  "<<prot->getNumResidues(0)<<"  ATOM NUMBER  "<<prot->getNumAtoms(0,0)<<prot->getNumAtoms(0,1)<<prot->getNumAtoms(0,2);

	//Calculate coordinates for the first residue, assumes ---INDEX starts from 0---
	dblVec vec(3);dblVec temp(3);
	vec[0]=0.0; vec[1]=0.0;vec[2]=0.0;
	prot->setCoords(0,0,0,vec);
	cout<<vec[0]<<"   "<<vec[1]<<endl;
	vec[0]=bond_length[0]*cos(36.06*PI/180);
	vec[1]=bond_length[0]*sin(36.06*PI/180);
	prot->setCoords(0,0,1,vec);
	cout<<vec[0]<<"   "<<vec[1]<<endl;
	vec[0]=vec[0]+(bond_length[1]*cos(-(180-(36.06+bond_angle[1]))*PI/180));
	vec[1]=0;//vec[1]+(bond_length[1]*sin(-(180-(36.06+bond_angle[1]))*PI/180));
	temp[0]=vec[0];temp[1]=vec[1];temp[2]=vec[2];
	prot->setCoords(0,0,2,vec);
	cout<<vec[0]<<"   "<<vec[1]<<endl;
	//TEMP FOR ATOM 3= O
	temp[0]=temp[0]+(bond_length[2]*cos((bond_angle[2]+(36.06+bond_angle[1]))*PI/180));
	temp[1]=temp[1]+(bond_length[2]*sin((bond_angle[2]+(36.06+bond_angle[1]))*PI/180));
	prot->setCoords(0,0,3,temp);
	cout<<temp[0]<<"   "<<temp[1]<<endl;
	//First Atom of the second residue
	cout<<endl<<endl<<endl;
	vec[0]=vec[0]+(bond_length[3]*cos(bond_x[3]*PI/180));
	vec[1]=vec[1]+(bond_length[3]*sin(bond_x[3]*PI/180));
	prot->setCoords(0,1,0,vec);
	cout<<vec[0]<<"   "<<vec[1]<<endl;

	//Coordinates for the second residue on
	for(UInt i=1;i<prot->getNumResidues(0);i++)
	{
		for(UInt j=1;(j<prot->getNumAtoms(0,i));j++)
		{
			if(j<2)
			{
			    if(i%2==0)
			    {	
				vec[0]=vec[0]+(bond_length[j-1]*cos(bond_x[j-1]*PI/180));
        			vec[1]=vec[1]+(bond_length[j-1]*sin(bond_x[j-1]*PI/180));
				prot->setCoords(0,i,j,vec);
				cout<<"j= "<<j<<"  "<<vec[0]<<"   "<<vec[1]<<"  "<<endl;
			    }
			    else
			    {
				vec[0]=vec[0]+(bond_length[j-1]*cos(-bond_x[j-1]*PI/180));
                                vec[1]=vec[1]+(bond_length[j-1]*sin(-bond_x[j-1]*PI/180));
				prot->setCoords(0,i,j,vec);
				cout<<"j= "<<j<<"  "<<vec[0]<<"   "<<vec[1]<<endl;
			    }
			}//if (j<2)
			else
			{
			    if(j==3)
			    {
			    	cout<<"Last ATOM"<<endl;
				if(i%2==0)
				{
					temp[0]=vec[0]+(bond_length[2]*cos((bond_angle[2]+bond_x[0]+bond_angle[1])*PI/180));
                                        temp[1]=vec[1]+(bond_length[2]*sin((bond_angle[2]+bond_x[0]+bond_angle[1])*PI/180));
                                        prot->setCoords(0,i,3,temp);
					cout<<temp[0]<<"   "<<temp[1]<<endl;
				}
				else
				{
					temp[0]=vec[0]+(bond_length[2]*cos((180-(bond_angle[2]-(180-(bond_x[0]+bond_angle[1]))))*PI/180));
					temp[1]=vec[1]+(bond_length[2]*sin((180-(bond_angle[2]-(180-(bond_x[0]+bond_angle[1]))))*PI/180));;
					prot->setCoords(0,i,3,temp);
					cout<<"j= "<<j<<"  "<<temp[0]<<"   "<<temp[1]<<"   "<<endl;
				}
			    }
			    else
			    {
				if(i%2==0)
                         	{
                        		 vec[0]=vec[0]+(bond_length[j-1]*cos(-bond_x[j-1]*PI/180));
                                	 vec[1]=vec[1]+(bond_length[j-1]*sin(-bond_x[j-1]*PI/180));
                                	prot->setCoords(0,i,j,vec);
				 	 cout<<"j= "<<j<<"  "<<vec[0]<<"   "<<vec[1]<<endl;
                            	}
                            	else
                            	{
                                	vec[0]=vec[0]+(bond_length[j-1]*cos(bond_x[j-1]*PI/180));
                                	vec[1]=vec[1]+(bond_length[j-1]*sin(bond_x[j-1]*PI/180));
                                	prot->setCoords(0,i,j,vec);
					cout<<"j= "<<j<<"  "<<vec[0]<<"   "<<vec[1]<<endl;
                            	}
			   }
			}
		}//j LOOP
			//First atom of next residue
			cout<<endl<<endl<<endl;
			if(i<(prot->getNumResidues(0)-1))

			{
			    if(i%2==0)
			    {
			    	vec[0]=vec[0]+(bond_length[3]*cos(bond_x[3]*PI/180));
                                vec[1]=vec[1]+(bond_length[3]*sin(bond_x[3]*PI/180));
                                prot->setCoords(0,i+1,0,vec);
				cout<<vec[0]<<"   "<<vec[1]<<endl;
			    }		    
			    else
			    {
				vec[0]=vec[0]+(bond_length[3]*cos(-bond_x[3]*PI/180));
                                vec[1]=vec[1]+(bond_length[3]*sin(-bond_x[3]*PI/180));
                                prot->setCoords(0,i+1,0,vec);
                                cout<<vec[0]<<"   "<<vec[1]<<endl;
			    }
			}
	}//i LOOP
pdbWriter(prot, "output.pdb");

return 0;
}
