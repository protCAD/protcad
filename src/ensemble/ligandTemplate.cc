#include "ligandTemplate.h"

UInt ligandTemplate::howManyTemplates=0;
UInt ligandData::howManyLigData=0;

ligandTemplate::ligandTemplate()
{
	//cout << "ligandTemplate constructor called" << endl;
	initialize();
        howManyTemplates++;
        //cout << "Leaving ligandTemplate() constructor" << endl;
}

ligandTemplate::~ligandTemplate()
{
    howManyTemplates--;
}

void ligandTemplate::initialize()
{
	//cout << "called ligandTemplate::initialize()" << endl;
        itsLigandDataTypes.resize(0);
}

void ligandTemplate::readLigandLibFile(const vector<string> _theLines)
{
    //cout << "In ligandTemplate::readLigandLibFile()" << endl;
    
    string linebuffer;
    int totalLines=_theLines.size();
    int j=0;
    string tempSize="";
    int NumAtomSize=0;
    bool connected=false;
    UInt parseCounter=0;  
        //0=atom definitions with 8 columns, etc
        //1= independent atoms
        //2= linked independent atoms
        //3= Main Connection Map
        //4= Hydrogen Connection
    vector<string> independentAtoms;
    vector<vector<string> > linkedIndependentAtoms;
    vector<vector<string> > mainConnectivity;
    vector<vector<string> > hydrogenConnectivity;
    vector<string> tempVector;
        
    ligandData* tempLigData= new ligandData();
    
    while(j<totalLines)
    {
        if(parseCounter==0)
        {
            linebuffer=_theLines[j];
            
            //Parse Counter=0... atom definitions
            if(linebuffer.substr(0,1)== "!")
            {
                // Set type, numAtoms, and Connectivity from ! line
                tempVector=getTypeFromFile(linebuffer);
                tempLigData->setItsType(tempVector[0]);
                //cout <<"itsType=*" << tempVector[0]<< "*"<<endl;
                
                sscanf(tempVector[1].c_str(),"%d",&NumAtomSize);
                tempLigData->setNumAtoms(NumAtomSize);
                //cout << "NumAtomSize=" << tempLigData->getNumAtomsInTemplate(true) << endl;
                
                if(tempVector[2]=="Y"){connected=true;}
                tempLigData->setConnectivity(connected);
                //cout << "Connectivity=" << connected << endl;
                
                tempVector.resize(0);
                
                //Parse the following NumAtomSize Lines
                for(int i=1; i <= NumAtomSize; i++)
                {
                    linebuffer=_theLines[j+i];
                    tempVector=parseTabbedLine(linebuffer);
                    tempLigData->pushLigAtomData(tempVector);
                    tempVector.resize(0);
                }
                
                j+=(NumAtomSize+1);
                parseCounter++;
                
                if(connected==false){parseCounter=5; j=totalLines;}//skip connectivity building
            }
    
            else
            {
                cout << "Error in reading a .lig file in ligandTemplate::readLigandLibFile()!" << endl << "The program will now quit.  Please check the formatting on the file." << endl;
                exit(1);
            }
        } //end parseCounter==0
    
        
        
        if(parseCounter==1)
        {
            linebuffer=_theLines[j];
            
            tempVector=parseSpacedLine(linebuffer);
            if(tempVector[0]== "!INDEPENDENT_ATOMS")
            {
                UInt numIndeAtoms;
                
                sscanf(tempVector[1].c_str(),"%u",&numIndeAtoms);
                tempVector.resize(0);
                
                for(UInt i=0; i<numIndeAtoms; i++)
                {	
                    tempVector=parseSpacedLine(_theLines[j+i+1]);
                    
                    independentAtoms.push_back(tempVector[0]); //there should be only one element
                    
                    tempVector.resize(0);
                }
                
                j+=(numIndeAtoms+1);
                parseCounter++;
            }
            
            else
            {
                cout << "ReadLigandLibFile():: Independent_Atoms line formatted incorrectly."<<endl;
            }
        
        } //end parseCounter==1
        
        if(parseCounter==2)
        {
            linebuffer=_theLines[j];
            
            tempVector=parseSpacedLine(linebuffer);
            
            if(tempVector[0]== "!LINKED_INDEPENDENT_ATOMS")
            {
                UInt numLinkedIndeAtoms;
                
                sscanf(tempVector[1].c_str(),"%u",&numLinkedIndeAtoms);
                tempVector.resize(0);
                
                for(UInt i=0; i<numLinkedIndeAtoms; i++)
                {	
                    tempVector=parseSpacedLine(_theLines[j+i+1]);
                    
                    linkedIndependentAtoms.push_back(tempVector); //pushes whole parsed line
                    
                    tempVector.resize(0);
                }
                
                j+=(numLinkedIndeAtoms+1);
                parseCounter++;
            }
            
            else
            {
                cout << "ReadLigandLibFile():: Linked_Independent_Atoms line formatted incorrectly."<<endl;
            }
        
        } //end parseCounter==2
        
        if(parseCounter==3)
        {
            linebuffer=_theLines[j];
            
            tempVector=parseSpacedLine(linebuffer);
            
            if(tempVector[0]== "!MAIN_CONNECT")
            {
                UInt numMainConnectAtoms;
                sscanf(tempVector[1].c_str(),"%u",&numMainConnectAtoms);
                tempVector.resize(0);
                
                for(UInt i=0; i<numMainConnectAtoms; i++)
                {	
                    tempVector=parseSpacedLine(_theLines[j+i+1]);
                    
                    mainConnectivity.push_back(tempVector); //pushes whole parsed line
                    
                    tempVector.resize(0);
                }
                
                j+=(numMainConnectAtoms+1);
                parseCounter++;
            }
            
            else
            {
                cout << "ReadLigandLibFile():: Main_Connect line formatted incorrectly."<<endl;
            }
        
        } //end parseCounter==3
        
        if(parseCounter==4)
        {
            linebuffer=_theLines[j];
            
            tempVector=parseSpacedLine(linebuffer);
            
            if(tempVector[0]== "!HYDRO_CONNECT")
            {
                UInt numHydroConnectAtoms;
                sscanf(tempVector[1].c_str(),"%u",&numHydroConnectAtoms);
                tempVector.resize(0);
                
                for(UInt i=0; i<numHydroConnectAtoms; i++)
                {	
                    tempVector=parseSpacedLine(_theLines[j+i+1]);
                    
                    hydrogenConnectivity.push_back(tempVector); //pushes whole parsed line
                    
                    tempVector.resize(0);
                }
                
                j+=(numHydroConnectAtoms+1);
                parseCounter++;
            }
            
            else
            {
                cout << "ReadLigandLibFile():: Hydro_Connect line formatted incorrectly."<<endl;
            }
        
        } //end parseCounter==4 ... parse counter should now equal 5.
       
       }//end of outside while loop
     
    if((parseCounter!=5)&&(connected))
    {
        cout << "Error in readLigandLibFile().  parseCounter!=5."<<endl;
        cout << "Please check formatting in the .lig files.   Quitting out..."<<endl;
        exit(1);
    }
        
        
    if(connected)
    {
        tempLigData->setConnectivity(true);
        tempLigData->setIndependentAtoms(independentAtoms);
        tempLigData->setLinkedIndependentAtoms(linkedIndependentAtoms);
        tempLigData->setMainConnect(mainConnectivity);
        tempLigData->setHydroConnect(hydrogenConnectivity);
    }
    
    //tempLigData->printVecSizes(); //for debugging .lig file atom properties
    //tempLigData->printConnectivity();
     
    itsLigandDataTypes.push_back(tempLigData);
    
    cout << "Setup " << tempLigData->getItsType() << " OK." << endl;
          
}

vector<string> ligandTemplate::getTypeFromFile(string _linebuffer)
{
    //This will parse the ! line and push the types into an array
    
    string tempBuffer=_linebuffer;
    string itemString="";
    string tempstring;
    vector<string> tempVector;
    
    //Starts from the 3rd column in the line (so "!<space>NAME is very important)
    for(UInt i=2; i<tempBuffer.size(); i++)
    {
        tempstring=tempBuffer.substr(i,1);
        
        if(tempstring != "\t"){itemString+=tempstring;}
        
        else
        {  
            tempVector.push_back(itemString);
            itemString="";
        }
    }

    // Push last element... assumes line doesn't end in space
    tempVector.push_back(itemString);
    
    return tempVector;
    
}

vector<string> ligandTemplate::parseSpacedLine(string _linebuffer)
{
    string tempstring;
    vector<string> tempVector;
    
    for(UInt i=0; i<_linebuffer.size(); i++)
    {
        string singlet=_linebuffer.substr(i,1);
        if(singlet != " "){tempstring+=singlet;}
        
        else
        {
            tempVector.push_back(tempstring);
            tempstring="";
        }
    }
    
    //push last element... assumes line doesn't end in a space
    tempVector.push_back(tempstring);
    
    return tempVector;

}

vector<string> ligandTemplate::parseTabbedLine(string _linebuffer)
{
    //cout << "Starting parseTabbedLine()" << endl;
    string tempBuffer=_linebuffer;
    string itemString="";
    string tempstring;
    vector<string> tempVector;
    
    for(UInt i=0; i<tempBuffer.size(); i++)
    {
        tempstring=tempBuffer.substr(i,1);
        
        if(tempstring != "\t"){itemString+=tempstring;}
        
        else
        {  
            tempVector.push_back(itemString);
            itemString="";
        }
    }
    
    //push last item
    tempVector.push_back(itemString);
    
    return tempVector;
}

int ligandTemplate::scanTemplates(string _ligName)
{
    int numTemplates=itsLigandDataTypes.size();
    string tempLigName="";
    
    for(int i=0; i<numTemplates; i++)
    {
        tempLigName=itsLigandDataTypes[i]->getItsType();
        if(tempLigName==_ligName){return i;}
        tempLigName+=" "; //allows for 3 character names, instead of 4
        if(tempLigName==_ligName){return i;}
    }
    
    cout << "Error in ligandTemplate::scanTemplates()." << endl;
    cout << "No template for ligand of type " << _ligName << ".  Quitting out." << endl;
    exit(1);
    
    return -9999;
}

int ligandTemplate::scanTemplatesForAtomLocation(int _LigDataIndex, string _AtomName)
{
    
    int LigDataAtomIndex= itsLigandDataTypes[_LigDataIndex]->getAtomIndex(_AtomName);
    
    if(LigDataAtomIndex==-9999)
    {
        cout << "Error in ligandTemplate::scanTemplatesForAtomLocation." << endl;
        cout << "No atom of type " << _AtomName << ".  Quitting out." << endl;
        exit(1);
    }
    
    return LigDataAtomIndex;
}
    

void ligandTemplate::printLigTemplates()
{
    for(UInt i=0; i< itsLigandDataTypes.size(); i++)
    {
        cout << "LigandType=*" << itsLigandDataTypes[i]->getItsType() << "*... NumAtoms=*" <<itsLigandDataTypes[i]->getNumAtomsInTemplate(false)<<"*" <<endl;
    }
}
string ligandTemplate::getAmberAllTypeName(int _LigDataIndex, int _LigDataAtomIndex)
{

    string AmberType=(itsLigandDataTypes[_LigDataIndex])->getAmberAllName(_LigDataAtomIndex);
    
    return AmberType;

}

string ligandTemplate::getAmberUnitedTypeName(int _LigDataIndex, int _LigDataAtomIndex)
{

    string AmberUnitedType=(itsLigandDataTypes[_LigDataIndex])->getAmberUnitedName(_LigDataAtomIndex);
    
    return AmberUnitedType;

}

double ligandTemplate::getAmberAllCharge(int _LigDataIndex, int _LigDataAtomIndex)
{
    double AmberAllCharge=(itsLigandDataTypes[_LigDataIndex])->getAmberAllCharge(_LigDataAtomIndex);
    
    return AmberAllCharge;
}

double ligandTemplate::getAmberUnitedCharge(int _LigDataIndex, int _LigDataAtomIndex)
{
    double AmberUnitedCharge=(itsLigandDataTypes[_LigDataIndex])->getAmberUnitedCharge(_LigDataAtomIndex);
    
    return AmberUnitedCharge;
}

string ligandTemplate::getSummaAtomType(int _LigDataIndex, int _LigDataAtomIndex)
{
    string SummaAtomType=(itsLigandDataTypes[_LigDataIndex])->getSummaAtomType(_LigDataAtomIndex);
    
    return SummaAtomType;
}

string ligandTemplate::getSummaEnvType(int _LigDataIndex, int _LigDataAtomIndex)
{
    string SummaEnvType=(itsLigandDataTypes[_LigDataIndex])->getSummaEnvType(_LigDataAtomIndex);
    
    return SummaEnvType;
}

string ligandTemplate::getSixAtomSolvationType(int _LigDataIndex, int _LigDataAtomIndex)
{
    string SixAtomSolvationType=(itsLigandDataTypes[_LigDataIndex])->getSixAtomType(_LigDataAtomIndex);
    
    return SixAtomSolvationType;
}

void ligandTemplate::printConnectivity(UInt _index)
{
    if(!itsLigandDataTypes[_index]->getIsConnected())
    {
        cout << "ligandTemplate::printConnectivity()... no connectivity to print for " << itsLigandDataTypes[_index]->getItsType()<< endl;
        return;
    }
    
    cout << "*******************************" << endl;
    cout << "*******************************" << endl;
    cout << "Connectivity list for " << itsLigandDataTypes[_index]->getItsType() << endl;
    cout << "-------------------------------" << endl;
    itsLigandDataTypes[_index]->printConnectivity();
    cout << "*******************************" << endl;
    cout << "*******************************" << endl;
}

void ligandTemplate::printConnectivity()
{
    for(UInt i=0; i<itsLigandDataTypes.size(); i++)
    {
        printConnectivity(i);
    }
}

// ***************************************************************************************
// ***************************************************************************************
//
// Start ligandData Class.  Please put all ligandTemplate functions above this section.
//
// ***************************************************************************************
// ***************************************************************************************

ligandData::ligandData()
{
    //cout << "Default ligandData constructor called." << endl;
    LigDataInitialize();
    howManyLigData++;
}

ligandData::~ligandData()
{
    howManyLigData--;
}

void ligandData::LigDataInitialize()
{    
    itsTypeString="";
    itsNumAtomsInTemplate=0;
    itsNumHydrogens=0;
    itsAtomNameList.resize(0);
    itsAmberAllAtomNames.resize(0);
    itsAmberUnitedAtomNames.resize(0);
    itsSummaAtomTypes.resize(0);
    itsSummaEnvTypes.resize(0);
    itsSixAtomSolvationTypes.resize(0);
    itsAmberAllAtomCharges.resize(0);
    itsAmberUnitedAtomCharges.resize(0);

    hasConnectivity=false;
    itsIndependentAtoms.resize(0);
    itsLinkedIndependentAtoms.resize(0);
    itsMainConnect.resize(0);
    itsHydroConnect.resize(0);
}

void ligandData::pushLigAtomData(vector<string> _parsedLigAtomLine)
{
    // Vector should have 8 elements
    // 8 columns in order from left to right
    
    //Error Check
    if(_parsedLigAtomLine.size()!=8)
    {
        cout << "Error in pushLigAtomData(), Vector is not size 8." << endl;
        cout << "Quitting out..." << endl;
        exit(1);
    }

    itsAtomNameList.push_back(_parsedLigAtomLine[0]);  
    if((_parsedLigAtomLine[0].substr(0,1))=="H"){itsNumHydrogens++;}
    
    itsAmberAllAtomNames.push_back(_parsedLigAtomLine[1]); 
    itsAmberUnitedAtomNames.push_back(_parsedLigAtomLine[2]); 
    itsSummaAtomTypes.push_back(_parsedLigAtomLine[3]); 
    itsSummaEnvTypes.push_back(_parsedLigAtomLine[4]);  
    itsSixAtomSolvationTypes.push_back(_parsedLigAtomLine[5]); 
    
  
    double numZero=0.0;
    if((_parsedLigAtomLine[6])=="-"){itsAmberAllAtomCharges.push_back(numZero);}
    else
    {
        double tempDouble=0.0;
        sscanf(_parsedLigAtomLine[6].c_str(),"%lf",&tempDouble);
        itsAmberAllAtomCharges.push_back(tempDouble);
    } 
    
    if((_parsedLigAtomLine[7])=="-"){itsAmberUnitedAtomCharges.push_back(numZero);}
    else
    {
        double tempDouble2=0.0;
        sscanf(_parsedLigAtomLine[7].c_str(),"%lf",&tempDouble2);
        itsAmberUnitedAtomCharges.push_back(tempDouble2);
    } 
    
}

UInt ligandData::getNumAtomsInTemplate(bool _hydrogens)
{
    if(_hydrogens)//if hydrogens are on
    {
        return itsNumAtomsInTemplate;
    }
    
    else
    {
        return itsNumAtomsInTemplate-itsNumHydrogens;
    }

}


string ligandData::getAtomName(UInt _index)
{
    if(_index < itsAtomNameList.size()){return itsAtomNameList[_index];}
    else{
        cout << "_index too large in getAtomName(). Quitting..." << endl;
        exit(1);
    }
}

int ligandData::getAtomIndex(string _atomName)
{
    UInt i=0;
    
    for(; i<itsNumAtomsInTemplate; i++)
    {
        if(_atomName==getAtomName(i))
        {
            return (int)i;
        }
    }
    
    return -9999;
}

string ligandData::getAmberAllName(UInt _index)
{
    if(_index < itsAmberAllAtomNames.size()){return itsAmberAllAtomNames[_index];}
    else{
        cout << "_index too large in getAmberAllName(). Quitting..." << endl;
        exit(1);
    }
}

string ligandData::getAmberUnitedName(UInt _index)
{
    if(_index < itsAmberUnitedAtomNames.size()){return itsAmberUnitedAtomNames[_index];}
    else{
        cout << "_index too large in getAmberUnitedName(). Quitting..." << endl;
        exit(1);
    }
}

string ligandData::getSummaAtomType(UInt _index)
{
    if(_index < itsSummaAtomTypes.size()){return itsSummaAtomTypes[_index];}
    else{
        cout << "_index too large in getSummaAtomType(). Quitting..." << endl;
        exit(1);
    }
}

string ligandData::getSummaEnvType(UInt _index)
{
    if(_index < itsSummaEnvTypes.size()){return itsSummaEnvTypes[_index];}
    else{
        cout << "_index too large in getSummaEnvType(). Quitting..." << endl;
        exit(1);
    }
}

string ligandData::getSixAtomType(UInt _index)
{
    if(_index < itsSixAtomSolvationTypes.size()){return itsSixAtomSolvationTypes[_index];}
    else{
        cout << "_index too large in getSixAtomType(). Quitting..." << endl;
        exit(1);
    }
}

double ligandData::getAmberAllCharge(UInt _index)
{
    if(_index < itsAmberAllAtomCharges.size()){return itsAmberAllAtomCharges[_index];}
    else{
        cout << "_index too large in getAmberAllCharge(). Quitting..." << endl;
        exit(1);
    }
}

double ligandData::getAmberUnitedCharge(UInt _index)
{
    if(_index < itsAmberUnitedAtomCharges.size()){return itsAmberUnitedAtomCharges[_index];}
    else{
        cout << "_index=" << _index << " too large in getAmberUnitedCharge(). Quitting..." << endl;
        exit(1);
    }
}

void ligandData::printVecSizes()
{
    
    cout << "size of unitedatomcharge=" << itsAmberUnitedAtomCharges.size() << endl;
    cout << "size of unitedatomtype=" << itsAmberUnitedAtomNames.size()<<endl;
    cout << "size of all atom charge=" << itsAmberAllAtomCharges.size()<< endl;
    cout << "size of all atom names=" << itsAmberAllAtomNames.size()<< endl;
    cout << "size of atomnamelist=" << itsAtomNameList.size()<<endl;
    
    cout << "AmberAllCharges..." <<endl;
    for(UInt i=0; i< itsAmberAllAtomCharges.size(); i++)
    {
        cout << itsAmberAllAtomCharges[i] << ", ";
    }
    
    cout << endl<< "AmberUnitedCharges"<<endl;
    for(UInt i=0; i< itsAmberUnitedAtomCharges.size(); i++)
    {
        cout << itsAmberUnitedAtomCharges[i] << ", ";
    }
    cout << endl;
    
    

}

void ligandData::printConnectivity()
{
    cout << "Independent Atoms.  NumAtoms=" <<itsIndependentAtoms.size()<< endl;
    for(UInt i=0; i<itsIndependentAtoms.size();i++)
    {
        cout << itsIndependentAtoms[i];
        if(i!=(itsIndependentAtoms.size()-1)){cout<<",";}
    }
    
    cout << endl << "Independent Atom Linkage. pairs="<< itsLinkedIndependentAtoms.size() <<endl;
    for(UInt i=0; i<itsLinkedIndependentAtoms.size();i++)
    {
        //cout << "j=" << itsLinkedIndependentAtoms[i].size() << ".  ";
        for(UInt j=0; j<itsLinkedIndependentAtoms[i].size();j++)
        {
            cout <<itsLinkedIndependentAtoms[i][j];
            if(j!=(itsLinkedIndependentAtoms[i].size()-1)){cout<<",";}
        }
        cout << endl;
    }
    
    cout << endl << "Main Connectivity.  TotalConnectLines=" << itsMainConnect.size()<<endl;
    for(UInt i=0; i<itsMainConnect.size();i++)
    {
        for(UInt j=0; j<itsMainConnect[i].size();j++)
        {
            cout <<itsMainConnect[i][j];
            if(j!=(itsMainConnect[i].size()-1)){cout<<",";}
        }
        cout << endl;
    }
    
    cout << endl << "Hydrogen Connectivity.  TotalHConnect="<< itsHydroConnect.size()<<endl;
    for(UInt i=0; i<itsHydroConnect.size();i++)
    {
        for(UInt j=0; j<itsHydroConnect[i].size();j++)
        {
            cout <<itsHydroConnect[i][j];
            if(j!=(itsHydroConnect[i].size()-1)){cout<<",";}
        }
        cout << endl;
    }
    cout <<endl;
}

    
