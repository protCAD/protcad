#include "assert.h"
#include "atom.h"
#include "typedef.h"

#ifndef LIGAND_H
class ligand;
#endif



#ifndef LIGANDDATA_H
#define LIGANDDATA_H

class ligandData
{
    public:
        // Constructors and Destructors
        ligandData();
        ~ligandData();
        
        // Setup Functions
        void LigDataInitialize();
        void pushLigAtomData(vector<string> _parsedLigAtomLine);
        
        //Accessor Functions
        string getItsType(){return itsTypeString;};
        void setItsType(string _type){itsTypeString=_type;};
        UInt getNumAtomsInTemplate(bool _hydrogens);
        void setNumAtoms(UInt _numatoms){itsNumAtomsInTemplate=_numatoms;};
        int getAtomIndex(string _atomName);
        
        string getAtomName(UInt _index);
        string getAmberAllName(UInt _index);
        string getAmberUnitedName(UInt _index);
        string getSummaAtomType(UInt _index);
        string getSummaEnvType(UInt _index);
        string getSixAtomType(UInt _index);
        double getAmberAllCharge(UInt _index);
        double getAmberUnitedCharge(UInt _index);
        
        void printVecSizes(); //for debugging .lig files... prints atom properties
       
        //Connectivity Functions
        bool getIsConnected(){return hasConnectivity;};
        void setConnectivity(bool _status){hasConnectivity=_status;};
        void printConnectivity(); //for debugging purposes
        
        void setIndependentAtoms(vector<string> _indeAtoms){itsIndependentAtoms=_indeAtoms;};
        void setLinkedIndependentAtoms(vector<vector<string> > _linkedIndeAtoms){itsLinkedIndependentAtoms=_linkedIndeAtoms;};
        void setMainConnect(vector<vector<string> > _mainConnect){itsMainConnect=_mainConnect;};
        void setHydroConnect(vector<vector<string> > _hydroConnect){itsHydroConnect=_hydroConnect;};
        
        vector<string> getIndependentAtoms(){return itsIndependentAtoms;};
        vector<vector<string> > getLinkedIndependentAtoms(){return itsLinkedIndependentAtoms;};
        vector<vector<string> > getMainConnect(){return itsMainConnect;};
        vector<vector<string> > getHydroConnect(){return itsHydroConnect;};
        
    private:
    	
        // Variable declarations
	string itsTypeString;
        UInt itsNumAtomsInTemplate;
        UInt itsNumHydrogens;
	vector<string> itsAtomNameList;
        vector<string> itsAmberAllAtomNames;
        vector<string> itsAmberUnitedAtomNames;
        vector<string> itsSummaAtomTypes;
        vector<string> itsSummaEnvTypes;
        vector<string> itsSixAtomSolvationTypes;
        vector<double> itsAmberAllAtomCharges;
        vector<double> itsAmberUnitedAtomCharges;
        
        //Connectivity Variables
        bool hasConnectivity;
        vector<string> itsIndependentAtoms;
        vector<vector<string> > itsLinkedIndependentAtoms;
        vector<vector<string> > itsMainConnect;
        vector<vector<string> > itsHydroConnect;
        
        static UInt howManyLigData;
};
#endif


#ifndef LIGANDTEMPLATE_H
#define LIGANDTEMPLATE_H

class ligandTemplate
{
public:
	// Constructors and Destructors
        ligandTemplate();
	~ligandTemplate(); 
	
        // Setup Functions
        void initialize();
        void readLigandLibFile(const vector<string> _theLines);
        vector<string> getTypeFromFile(string _linebuffer);
        vector<string> parseTabbedLine(string _linebuffer);
        vector<string> parseSpacedLine(string _linebuffer);

        // Accessor Functions
        int scanTemplates(string _ligName);
        int scanTemplatesForAtomLocation(int _LigDataIndex, string _AtomName);
        void printLigTemplates();
        
        string getAmberAllTypeName(int _LigDataIndex, int _LigDataAtomIndex);
        string getAmberUnitedTypeName(int _LigDataIndex, int _LigDataAtomIndex);
        double getAmberAllCharge(int _LigDataIndex, int _LigDataAtomIndex);
        double getAmberUnitedCharge(int _LigDataIndex, int _LigDataAtomIndex);
        string getSummaAtomType(int _LigDataIndex, int _LigDataAtomIndex);
        string getSummaEnvType(int _LigDataIndex, int _LigDataAtomIndex);
        string getSixAtomSolvationType(int _LigDataIndex, int _LigDataAtomIndex);
        
        void printConnectivity(UInt _index);
        void printConnectivity();


        // Variables
        vector<ligandData*> itsLigandDataTypes;
        static UInt howManyTemplates;
};
#endif
