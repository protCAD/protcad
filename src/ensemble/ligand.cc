// 	_    _ ____ ____ _  _ ___   ____ ___  ___  
//	|    | | __ |__| |\ | |  \  |    |__] |__] 
//	|___ | |__] |  | | \| |__/ .|___ |    |    

// filename: ligand.cpp
// contents: class ligand implementation

#include "ligand.h"

UInt ligand::howMany=0;
bool ligand::dataBaseBuilt = false;
ligandTemplate ligand::itsLigTemplate;

pmf ligand::itsPMF("PMF_hires_symmetric.dat");
amberVDW ligand::itsAmberVDW(0);
//amberElec ligand::itsAmberElec(0);
solvation ligand::itsSolvation;
double ligand::cutoffDistance=6.0;
double ligand::cutoffDistanceSquared=cutoffDistance*cutoffDistance;

ligand::ligand() : molecule()
{	cout<< "default ligand constructor called" << endl;
        if(!dataBaseBuilt){setupDataBase();}
	itsNameString = "NUL";
	itsAtoms.resize(0);
        itsChainID="";
        itsLigTemplateType=-9999;
        hydrogensOn=false;
        symmetryLinked.resize(0);
	setMoleculeType(2);
	howMany++;
        isConnected=false;
        itsHeadNodes.resize(0);
        headNodeConnectVec.resize(0);
        ///cout <<"Exiting Ligand Constructor #1\n";
}

ligand::ligand(const string& _name) : molecule(_name)
{	cout<< "ligand constructor for " << _name << " called" << endl;
        if(!dataBaseBuilt){setupDataBase();}
        itsNameString = _name;
	itsAtoms.resize(0);
        itsChainID="";
        itsLigTemplateType=-9999;
        hydrogensOn=false;
        symmetryLinked.resize(0);
	setMoleculeType(2);
	howMany++;
        isConnected=false;
        itsHeadNodes.resize(0);
        headNodeConnectVec.resize(0);
	//cout <<"Exiting Ligand Constructor #2\n";
}

ligand::ligand(const ligand& _rhs)
{	
        if(!dataBaseBuilt){setupDataBase();}
        
	for (UInt i=0; i<_rhs.itsAtoms.size(); i++)
	{	atom* ligatom = new atom( *(_rhs.itsAtoms[i]));
		add(ligatom);
	}
	setMoleculeType(2);
	itsNameString = _rhs.itsNameString;
        itsChainID=_rhs.itsChainID;
        itsLigTemplateType=_rhs.itsLigTemplateType;
        symmetryLinked=_rhs.symmetryLinked;
        hydrogensOn=_rhs.hydrogensOn;
        isConnected=_rhs.isConnected;
        itsHeadNodes=_rhs.itsHeadNodes;
        headNodeConnectVec=_rhs.headNodeConnectVec;
        //cout <<"Exiting Ligand Constructor #3\n";
}

ligand::~ligand()
{	cout<< "ligand destructor called " << endl;
	for (UInt i=0; i < itsAtoms.size(); i++)
	{	delete itsAtoms[i];
	}
        howMany--;
}

void ligand::setupDataBase()
{
        cout << "In ligand:setupDataBase()"<<endl;
        //ligandDataBase.resize(0);
        
        string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	string ligDataFolder="/data/lig/";
        string ligPath=path+ligDataFolder;
        
	string ligandList = "ligandList.lig";
    
        //Get names of .lig files from ligandList.lig
        vector<string> ligNames=readLigFile(ligPath,ligandList);
        
        //Open each .lig file and read the lines in
        for(UInt i=0; i<ligNames.size(); i++)
        {
            //strips off all lines starting with "#" and makes a vector<string> of all the remaining lines
            vector<string> passVector=readLigFile(ligPath, ligNames[i]);
            
            //Pass .lig vector to the LigandTemplate to build template of it
            itsLigTemplate.readLigandLibFile(passVector);      
        }
        
        dataBaseBuilt=true;
        //printLigDataBaseTypes(); // info for debugging
        //printLigDataBaseConnectivity(); // info for debugging
        cout << "leaving ligand::SetupDataBase()" << endl;
    
}

void ligand::MatchToTemplate()
{
    cout << "In ligand::MatchToTemplate()." << endl;
    
    if(itsLigTemplateType==-9999)
    {
        cout << "No ligand by this name in the database.  No reorder necessary." << endl;
        return;
    }

    
    //make sure number of atoms in lig== number of atoms in that ligands type
    UInt NumAtoms=itsAtoms.size();
    UInt NumTemplateAtoms=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getNumAtomsInTemplate(hydrogensOn);
    
    if(NumAtoms != NumTemplateAtoms)
    {
        cout << "Num Atoms in Ligand doesn't match num atoms in corresponding LigandTemplate" << endl;
        exit(1);
    }


    //Create a temp vector<atoms*> that is correctly ordered.
    vector<atom*> tempArray;
    string TempName, TempName2;
    bool matchComplete;
    
    for(UInt i=0; i < NumAtoms; i++)
    {
        TempName=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getAtomName(i);
        matchComplete=false;
        
        for(UInt j=0; j<NumAtoms; j++)
        {
            TempName2=itsAtoms[j]->getName();
            if(TempName==TempName2)
            {
                tempArray.push_back(itsAtoms[j]);
                j=NumAtoms;
                matchComplete=true;
                
            }
        }
        
        if(matchComplete==false)
        {
            cout << "Error in ligand:MatchToTemplate().  Atom name not in template." << endl;
            exit(1);
        }
    }
    
    itsAtoms=tempArray;
    
    //Read the template, and if connectivity information exists, build the connectivity.
    isConnected=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getIsConnected();
    cout << "IsConnected=" << isConnected<<endl;
    
    if(isConnected){buildConnectivity();}

}

vector<string> ligand::readLigFile(string _ligpath, string _fileName)
{
    	string iFile;
	ifstream inFile;
        
        iFile = _ligpath + _fileName;
	inFile.open(iFile.c_str());

	if(!inFile)
	{	cout << "Error: unable to open input file: "
			 << iFile << endl;
		exit (1);
	}

        else{cout << "Opened " << _fileName <<" OK" << endl;}

        vector<string> theLines;
        string linebuffer;
        string firstCharacter;
        UInt i=0;
        
	while (getline(inFile,linebuffer,'\n'))
	{
		firstCharacter = linebuffer.substr(0,1);
                if(firstCharacter != "#"){theLines.push_back(linebuffer);}
                i++;
	}
        
        inFile.close();
	inFile.clear();
        
        if(theLines.size()<=0)
        {
            cout << "Error, " << _fileName << " had no readable lines. Check formatting.  Quitting out..."<<endl;
            exit(1);
        }
        
        return theLines;
}


void ligand::buildConnectivity()
{

    if(!isConnected)
    {
        cout << "ligand::buildConnectivity(). No connectivity information exists for this ligand.  Quitting out..."<<endl;
        exit(1);
    }
    
    cout << "Building Connectivity for " << itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getItsType()<<"."<<endl;
    
    //Build Main Connections
    vector<vector<string> > mainConnectVec=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getMainConnect();
        
        for(UInt i=0; i<mainConnectVec.size(); i++)
        {
            childConnect(mainConnectVec[i]);
        }
    
    mainConnectVec.resize(0);
    
    //Build Hydrogen Connections
        if(hydrogensOn)
        {
            vector<vector<string> > hydroConnectVec=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getHydroConnect();
            for(UInt i=0; i<hydroConnectVec.size(); i++)
            {
                childConnect(hydroConnectVec[i]);
            }
            hydroConnectVec.resize(0);
        }
        
    //Build List of Independent Atoms and then link them
    itsHeadNodes=nameToIndexVec(itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getIndependentAtoms());
    
        //check to make sure that treeNode classifies these as headnodes
        for(UInt i=0; i<itsHeadNodes.size(); i++)
        {
            if(!itsAtoms[itsHeadNodes[i]]->isHeadNode())
            {
                cout << "An atom in the independent atom list is NOT a headnode. Quitting out..." << endl;
                exit(1);
            }
        }
        
        //Push Index numbers into headNodeConnectVec
        vector<vector<string> > IndeConnectVec=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getLinkedIndependentAtoms();
        
        for(UInt i=0; i<IndeConnectVec.size(); i++)
        {
            headNodeConnectVec.push_back(nameToIndexVec(IndeConnectVec[i]));
        
        }
    
    printLigConnectivity(); // info for debugging 
}

vector<UInt> ligand::nameToIndexVec(vector<string> _NameVec)
{
    
    vector<UInt> indexVec;
    indexVec.resize(0);
    
    for(UInt i=0; i<_NameVec.size(); i++)
    {
        indexVec.push_back(getAtomIndexFromName(_NameVec[i]));
    }
    
    return indexVec;
}

void ligand::childConnect(vector<string> _connectVec)
{
    //match atom names in _connectVec to atom indeces and push onto vector
    vector<UInt> IndexVec= nameToIndexVec(_connectVec);
    
    cout << "vec1 size= " << _connectVec.size() << " / vec2 size= " << IndexVec.size() << endl;
    
    //IndexVec[0] is the parent.  The others are children of this parent.
    for(UInt i=1; i<IndexVec.size(); i++)
    {
    
        if(i==1)
        {
            itsAtoms[IndexVec[0]]->setChild(itsAtoms[IndexVec[i]]);
        }
            
        else
        {
            itsAtoms[IndexVec[i-1]]->setNextSib(itsAtoms[IndexVec[i]]);
        }
        
    }
    
}
void ligand::printLigDataBaseTypes()
{
    itsLigTemplate.printLigTemplates();
}

void ligand::printLigDataBaseConnectivity()
{
    itsLigTemplate.printConnectivity();
}

void ligand::printLigConnectivity()
{
    cout << "******************************************"<< endl;
    cout << "Printing Actual Connectivity for Ligand " << itsNameString << "." << endl;
    
    //HeadNodes (aka: independent Nodes)
    cout << "Number of HeadNodes= " << itsHeadNodes.size() << endl;
    cout << "They are...";
    for(UInt i=0; i< itsHeadNodes.size(); i++)
    {
        cout << itsAtoms[itsHeadNodes[i]]->getName();
        if(i!=(itsHeadNodes.size()-1)){cout <<",";}
    }

    cout << endl;
    
    //Print Children of each headnode
    for(UInt i=0; i<itsHeadNodes.size(); i++)
    {
        cout << "Connectivity for the " << i << " headNode group ..." << endl;
        
        UInt j=itsHeadNodes[i];
        cout << itsAtoms[j]->getName() << "..." << itsAtoms[j]->getNumChildren()<<endl;
    }
    cout << endl;
    
    cout << "total connect info...." << endl;
    for (UInt i=0; i<itsAtoms.size(); i++)
    {
        cout << "Atom: " << itsAtoms[i]->getName() << " ...";
        itsAtoms[i]->printImmediateChildren();
        cout << endl;
    }
    
    cout << "******************************************"<< endl;
    
}

void ligand::add(atom* _ligAtom)
{	
        int numAtoms=atomCount();
    
        if(numAtoms==0)
        {
            //cout << "first atom in this ligand, setting type." << endl;
            itsNameString=_ligAtom->getResType();
            cout << "itsLigName=*" << itsNameString<<"*" <<endl;
            itsLigTemplateType=itsLigTemplate.scanTemplates(itsNameString);
        }
        
        // just checking for hydrogens...
        string tempName=_ligAtom->getName();
        if((tempName.substr(0,1))=="H"){hydrogensOn=true;}
        
        int atomLocation=itsLigTemplate.scanTemplatesForAtomLocation(itsLigTemplateType,tempName);
        
        //Set forcefield properties from ligandTemplate
        string amberAll= itsLigTemplate.getAmberAllTypeName(itsLigTemplateType,atomLocation);
        int tempInt= itsAmberVDW.getIndexFromNameString(amberAll);
        string amberUnited=itsLigTemplate.getAmberUnitedTypeName(itsLigTemplateType,atomLocation);
        int tempInt2= itsAmberVDW.getIndexFromNameString(amberUnited);
        
        _ligAtom->setAmberAllType(tempInt);
        _ligAtom->setAmberUnitedType(tempInt2);
        _ligAtom->setAmberAllCharge(itsLigTemplate.getAmberAllCharge(itsLigTemplateType,atomLocation));    
        _ligAtom->setAmberUnitedCharge(itsLigTemplate.getAmberUnitedCharge(itsLigTemplateType,atomLocation));
        
        itsAtoms.push_back(_ligAtom);
}


atom* ligand::getAtom(UInt _num)
{
	if(_num<=atomCount()){
		//cout << "Atom does exist, getting name..." << endl;
		//cout << itsAtoms[_num]->getName()<<endl;
		return itsAtoms[_num];
	}
	else{
		cout <<"********* Atom does not exist..."<<endl;
		return itsAtoms[0];
        }
}

UInt ligand::getAtomIndexFromName(string _AtomName)
{
    for(UInt i=0; i< itsAtoms.size(); i++)
    {
        string tempName=itsAtoms[i]->getName();
        
        if(_AtomName==tempName)
        {
            return i;
        }
    }
    
    cout << "Error in ligand::getAtomIndexFromName()... atom name not in ligand. Quitting..." << endl;
    exit(1);
    
    return 1;

}

vector<double> ligand::getCoords(UInt _atomIndex)
{
	if(_atomIndex > itsAtoms.size())
        {
            cout << "ligand::getCoords()... _atomIndex > itsAtoms.size().  Quitting out..." << endl;
            exit(1);
        }
        
        vector<double> coords(3);
    
	coords[0] = itsAtoms[_atomIndex]->getX();
	coords[1] = itsAtoms[_atomIndex]->getY();
	coords[2] = itsAtoms[_atomIndex]->getZ();
	
	if (coords.size() != 3) cout << "Coords for LigAtom #" << _atomIndex<< " not found ..." << endl;
	return coords;
}
//
// Modifiers
//

void ligand::translate(const dblVec& _dblVec)
{	for (UInt i = 0; i < itsAtoms.size(); i++)
	{	itsAtoms[i]->translate(_dblVec);
	}
	return;
}

void ligand::translate(const double _x,const double _y,const double _z)
{
	dblVec vec;
#ifdef USE_SVMT
	vec.resize(3);
#else
	vec.newsize(3);
#endif
	vec[0] = _x;
	vec[1] = _y;
	vec[2] = _z;
	translate(vec);
	return;
}

void ligand::transform(const dblMat& _dblMat)
{	for (UInt i = 0; i < itsAtoms.size(); i++)
	{	itsAtoms[i]->transform(_dblMat);
	}
	return;
}

void ligand::rotate(const axis _axis, const double _theta)
{       point origin;
        // The default is to set this point to the origin
        origin.setCoords(0.0,0.0,0.0);
        dblVec vec = dblVec(3);
        for (UInt i = 0; i<3; i++)
                vec[i] = 0.0;
        if (_axis == X_axis)
                vec[0] = 1.0;
        if (_axis == Y_axis)
                vec[1] = 1.0;
        if (_axis == Z_axis)
                vec[2]  = 1.0;
        rotate(origin,vec,_theta);
        return;
}

void ligand::rotate(const point& _point, const axis _axis, const double _theta)
{
        dblVec vec = dblVec(3);
        for (UInt i = 0; i<3; i++)
                vec[i] = 0.0;
        if (_axis == X_axis)
                vec[0] = 1.0;
        if (_axis == Y_axis)
                vec[1] = 1.0;
        if (_axis == Z_axis)
                vec[2]  = 1.0;
        rotate(_point,vec,_theta);
    return;
}

void ligand::rotate(const point& _point, const dblVec& _R_axis, const double _theta)
{
    dblVec toOrigin=_point.getCoords() * -1.0;
    dblVec backHome = _point.getCoords();
    dblMat R(3,3,0.0);
    R=CMath::rotationMatrix(_R_axis, _theta);
    for (UInt i = 0; i < itsAtoms.size(); i ++)
    {
        itsAtoms[i]->translate(toOrigin);
        itsAtoms[i]->transform(R);
        itsAtoms[i]->translate(backHome);
    }
    
    return;
}

//These rotation functions use connectivity provided via treeNode
void ligand::rotate(atom* _pAtom1, atom* _pAtom2, double _theta)
{    
    if(!isConnected)
    {
        cout << "Attempting to use a rotate command that requires connectivity information that has not been provided.  Use a non-connectivity command or add the proper information to the .lig file. Quitting out..."<<endl;
        exit(1);
    }
    
    dblVec toOrigin = _pAtom1->getCoords() * (-1.0);
    dblVec backHome = _pAtom1->getCoords();
    
    //Translate atom1, atom2, and atom2's children toOrigin
    _pAtom1->translate(toOrigin);
    _pAtom2->translate(toOrigin);
    _pAtom2->translateChildren(toOrigin);
    
    //calcualate rotation matrix based on vector defined by
    //point q1 (second point)
    dblMat R(3,3,0.0);
    R = CMath::rotationMatrix(_pAtom2->getCoords(), _theta);

    _pAtom2->transformChildren(R);
    
    //Finally, translate all the atoms back via the backHome vector
    _pAtom1->translate(backHome);
    _pAtom2->translate(backHome);
    _pAtom2->translateChildren(backHome);

}

void ligand::rotate(string _atom1, string _atom2, double _theta)
{
    //Note, atom1 stays in one place while atom2 and its children rotate
    
    //Find index of atoms in itsAtoms from names
    UInt index1=getAtomIndexFromName(_atom1);
    UInt index2=getAtomIndexFromName(_atom2);
    
    //Make sure the pair are headNodes
    if(!(isHeadNode(index1)) || !(isHeadNode(index2)))
    {
        cout << "Attempting to rotate around two atoms that are not defined as headNodes (independent) in the .lig file. Quitting out... " <<endl;
    }
    
   //Make sure the pair is in headNodeConnectVec
    if(!(isHeadNodePair(index1,index2)))
    {
        cout << "Attempting to rotate around two atoms which are not in the headNodeConnectVec.  Check .lig file or your atom names.  Quitting out..."<<endl;
        exit(1);
    }
    
    //pass atom pointers to primary rotate command
    rotate(itsAtoms[index1], itsAtoms[index2],_theta);
}

void ligand::rotate(UInt _atom1, UInt _atom2, double _theta)
{
    //Make sure the pair are headNodes
    if(!(isHeadNode(_atom1)) || !(isHeadNode(_atom2)))
    {
        cout << "Attempting to rotate around two atoms that are not defined as headNodes (independent) in the .lig file. Quitting out... " <<endl;
    }
    
    //Make sure the pair is in headNodeConnectVec
    if(!(isHeadNodePair(_atom1,_atom2)))
    {
        cout << "Attempting to rotate around two atoms which are not in the headNodeConnectVec.  Check .lig file or your atom names.  Quitting out..."<<endl;
        exit(1);
    }
    
    //pass atom pointers to primary rotate command
    rotate(itsAtoms[_atom1], itsAtoms[_atom2],_theta);
    
}

void ligand::rotate(UInt _rotatePair, double _theta)
{
    //for convenience... if you know which pair you are rotating from the .lig file
    
    UInt index1= headNodeConnectVec[_rotatePair][0];
    UInt index2= headNodeConnectVec[_rotatePair][1];
    
    //pass atom pointers to primary rotate command
    rotate(itsAtoms[index1], itsAtoms[index2],_theta);

}

//
//   Energy Functions
//

int ligand::getAmberAllType(UInt _index)
{
  return itsAtoms[_index]->getAmberAllType();
}

int ligand::getAmberUnitedType(UInt _index)
{
  return itsAtoms[_index]->getAmberUnitedType();
}

double ligand::getAmberElec(UInt _index)
{
    if(hydrogensOn)
        return itsAtoms[_index]->getAmberAllCharge();
    
    else
        return itsAtoms[_index]->getAmberUnitedCharge();
}

void ligand::printAmberTypes()
{
    cout << "Starting printAmberAllTypes()" << endl;
    
    cout << "LigandName= " << itsNameString << endl;
    
    if(itsLigTemplateType==-9999)
    {
        cout << "Error in ligandTemplate::scanTemplatesForAtomLocation." << endl;
        cout << "No ligand by this name in the database.  Quitting out." << endl;
        exit(1);
    }

    string tempName;
    int ambertypenum=0;
    int amberunitednum=0;
    
    for(UInt i=0; i < itsAtoms.size(); i++)
    {
        tempName=itsAtoms[i]->getName();
        cout << "AtomName= " << tempName;
        
        ambertypenum=getAmberAllType(i);
        amberunitednum=getAmberUnitedType(i);
        
        cout << "     AA= " << ambertypenum << "    UA= " << amberunitednum << endl;
    }
        

}

double ligand::intraEnergy() 
{
    double TotalEnergy=0;
    UInt atomSize=itsAtoms.size();
    UInt Atom1, Atom2;
    double distanceSQ;
    UInt i, j;
    
    if(hydrogensOn)
    {
        for(i=0; i < atomSize; i++)
        {
            Atom1=(UInt)getAmberAllType(i);
                
            for(j=i+1; j<atomSize; j++)
            {
                Atom2=(UInt)getAmberAllType(j);
            
                distanceSQ=itsAtoms[i]->distanceSquared(itsAtoms[j]);
            
               if(distanceSQ >4){TotalEnergy+=itsAmberVDW.getEnergySQ(Atom1,Atom2,distanceSQ);}
            }
        }
    }
    
    else
    {
        for(i=0; i < atomSize; i++)
        {
            Atom1=getAmberUnitedType(i);
                
            for(j=i+1; j<atomSize; j++)
            {
                Atom2=getAmberUnitedType(j);
            
                distanceSQ=itsAtoms[i]->distanceSquared(itsAtoms[j]);
                            
                if(distanceSQ >4){TotalEnergy+=itsAmberVDW.getEnergySQ(Atom1,Atom2,distanceSQ);}
                
            }
        }
    }
    
    return TotalEnergy;
}

/*double ligand::getInterEnergy(ligand* _other)
{
    double TotalEnergy=0;
    UInt atomSize=itsAtoms.size();
    UInt otherAtomSize=_other->itsAtoms.size();
    UInt Atom1, Atom2;
    double distanceSQ;
    
    //testing variables
    double prevEnergy=0;
    double difference=0;
    string atom1name, atom2name;
    
    for(UInt i=0; i < atomSize; i++)
    {
        if(hydrogensOn){Atom1=(UInt)getAmberAllType(i);}
        else{Atom1=(UInt)getAmberUnitedType(i);}
                
        for(UInt j=0; j<otherAtomSize; j++)
        {
            if(_other->hydrogensOn){Atom2=(UInt)_other->getAmberAllType(j);}
            else{Atom2=(UInt)_other->getAmberUnitedType(j);}
            
            distanceSQ=itsAtoms[i]->distanceSquared(_other->itsAtoms[j]);
            
            if(distanceSQ < cutoffDistanceSquared)
            {
                TotalEnergy+=itsAmberVDW.getEnergySQ(Atom1,Atom2,distanceSQ);
            
            
                //Lines from here to end of if statement are for testing
                difference=TotalEnergy-prevEnergy;
                //atom1name=itsLigTemplate.itsLigandDataTypes[itsLigTemplateType]->getAmberAllName(i);
                //atom2name=_other->itsLigTemplate.itsLigandDataTypes[_other->itsLigTemplateType]->getAmberAllName(j);
                
                
                //cout << i <<","<< j << " type-> "<<Atom1 <<":" << Atom2<< "  DistSQ= " << distanceSQ << " EnDif= " << difference <<endl;  
            }
                  
            prevEnergy=TotalEnergy;
       
         }//inside loop
    }//outside loop
        
    return TotalEnergy;
}*/

//
//   Surface Area and Volume Calculations
//

void ligand::initializeSpherePoints()
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		itsAtoms[i]->initializeSpherePoints();
	}
	return;
}

void ligand::removeIntraLigandSpherePoints()
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = 0; j < itsAtoms.size(); j ++)
		{
			itsAtoms[i]->removeSpherePoints(itsAtoms[j]);
		}
	}

	return;
}

void ligand::removeInterLigandSpherePoints(ligand* _other)
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = 0; j < _other->itsAtoms.size(); j ++)
		{
			itsAtoms[i]->removeSpherePoints(_other->itsAtoms[j]);
		}
	}
	return;
}

double ligand::tabulateSurfaceArea()
{
	double surfaceArea = 0.0;
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		surfaceArea += itsAtoms[i]->calculateExposedSASA();
	}

	return surfaceArea;
}

//
// Ligand-Ligand Connectivity Code
//

void ligand::removeLigFromOtherLigs()
{
    for(UInt i=0; i< symmetryLinked.size(); i++)
    {
        symmetryLinked[i]->removeLinkedLigand(this);
    }
    
    symmetryLinked.resize(0);
    
    cout << "Ligand pointer removed from all symmetry info" << endl;
}

void ligand::removeLinkedLigand(ligand* _ligPointer)
{
    vector<ligand*> tempLigVec;
    bool tripFlag=false;
    
    for(UInt i=0; i<symmetryLinked.size(); i++)
    {
        if(_ligPointer==symmetryLinked[i])
        {
            cout << "Ligand Found: removing it from linkage info" << endl;
            tripFlag=true;
        }
        else{ tempLigVec.push_back(symmetryLinked[i]);}
    }
    
    if(tripFlag==false)
    {
        cout << "This ligand is not symmetry linked to the specified ligand"<<endl;
        cout << "Quitting from ligand::removeLinkedLigand() in error." << endl;
        exit(1);
    }
    
    symmetryLinked=tempLigVec;
}

//
// Intra-Ligand Connectivity Code
//

bool ligand::isHeadNode(UInt _index)
{
    for(UInt i=0; i<itsHeadNodes.size(); i++)
    {
        if(itsHeadNodes[i]==_index){return true;}
    }
    
    return false;
}

bool ligand::isHeadNodePair(UInt _index1, UInt _index2)
{
    //Make sure the pair is in headNodeConnectVec
    for(UInt i=0; i<headNodeConnectVec.size(); i++)
    {
                if((headNodeConnectVec[i][0]==_index1)&&(headNodeConnectVec[i][1]==_index2))
                {return true;}
                
                if((headNodeConnectVec[i][1]==_index1)&&(headNodeConnectVec[i][0]==_index2))
                {return true;}
    }
    
    return false;
}

