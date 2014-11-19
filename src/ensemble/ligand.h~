// _    _ ____ ____ _  _ ___   _  _ 
// |    | | __ |__| |\ | |  \  |__| 
// |___ | |__] |  | | \| |__/ .|  | 

// filename: ligand.h
// contents: class ligand is defined


#ifndef ASSERT_H
#include "assert.h"
#endif

#include <string.h>
#include <vector>

#ifndef TYPEDEF_H
#include "typedef.h"
#endif
//#include "CMath.h"

#ifndef MOLECULE_H
#include "molecule.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIGANDTEMPLATE_H
#include "ligandTemplate.h"
#endif

#ifndef GENERALIO_H
#include "generalio.h"
#endif

#ifndef AMBER_VDW_H
#include "amberVDW.h"
#endif

//#ifndef AMBER_ELEC_H
//#include "amberElec.h"
//#endif

#ifndef SOLVATION_H
#include "solvation.h"
#endif

#ifndef PMF_H
#include "pmf.h"
#endif

#ifndef ATOMITERATOR_H
class atomIterator;
#endif

#ifndef LIGAND_H
#define LIGAND_H

class ligand : public molecule
{
public:
	friend class atomIterator;
	
	// Constructor and Destructor declaration
	ligand();
	ligand(const string& _name);
	ligand(const ligand& _rhs);
	~ligand();

        // Ligand DataBase Setup and Functions
        void setupDataBase();
        vector<string> readLigFile(string _ligpath, string _fileName);
        void MatchToTemplate();
        void printLigDataBaseTypes();
        void printLigDataBaseConnectivity();
        
        // Intra-Ligand Connectivity
        void buildConnectivity();
        void printLigConnectivity();
        vector<UInt> nameToIndexVec(vector<string> _NameVec);
        void childConnect(vector<string> _connectVec);
        bool isHeadNode(UInt _index);
        bool isHeadNodePair(UInt _index1, UInt _index2);

        
        // Accessors
        void add(atom* _ligAtom);
	static UInt getHowMany() {return howMany; }
	UInt atomCount() {return itsAtoms.size();}
	atom* getAtom(UInt _num);
        UInt getAtomIndexFromName(string _AtomName);
        vector<double> getCoords(UInt _atomIndex);	
        bool getHydrogensOn() const {return hydrogensOn;}
	void setHydrogensOn(const bool _hydrogensOn) ;
        void setItsNameString(const string _nameString) { itsNameString = _nameString; }
        string getItsNameString(){return itsNameString;}

        // Modifiers... these work on all atoms in a molecules
	void translate(const dblVec& _dblVec);
	void translate(const double _x,const double _y,const double _z);
	void transform(const dblMat& _dblMat);

        void rotate(const axis _axis,const double _theta);
	void rotate(const point& _pPoint, const dblVec& _R_axis, const double _theta);
	void rotate(const point& _point, const axis _axis, const double _theta);
        
        // Modifiers that use connectivity provided via treeNode
        void rotate(atom* _pAtom1, atom* _pAtom2, double _theta);
        void rotate(string _atom1, string _atom2, double _theta);
        void rotate(UInt _atom1, UInt _atom2, double _theta);
        void rotate(UInt _rotatePair, double _theta);
        
        // Energy Calculations
        void printAmberTypes();
        int getAmberAllType(UInt _index);
        int getAmberUnitedType(UInt _index);
        double getAmberElec(UInt _index);
        double intraEnergy();
        double getInterEnergy(ligand* _other);
        static void setCutoffDistance( const double _cutoff ) { cutoffDistance = _cutoff; cutoffDistanceSquared = _cutoff*_cutoff; }
        
        // Volume and Surface Area
        double tabulateSurfaceArea();
        void initializeSpherePoints();
        void removeIntraLigandSpherePoints();
        void removeInterLigandSpherePoints(ligand* _other);
        
        // Inter-Ligand (Ligand-Ligand) Connectivity Code
        vector<ligand*> getLinkedLigands(){return symmetryLinked;}
        void setLinkedLigand(ligand* _ligPointer){symmetryLinked.push_back(_ligPointer);}
        void removeAllLinkedInfo(){symmetryLinked.resize(0);}
        void removeLigFromOtherLigs(); //removes lig pointer from ligs it is linked to
        void removeLinkedLigand(ligand* _ligPointer);
        
        //Unused functions from molecule.  All set to blank inline until needed.
	void translate(const UInt _index, const dblVec& _dblVec) {}
        void translate(const UInt _index, const double _x,const double _y,const double _z) {}
        void transform(const UInt _index, const dblMat& _dblMat) {}
        void rotate(unsigned int, axis, double){}
        void rotate(unsigned int, const point &, const dblVec &, double){}
	
	double getPositionEnergy(vector <int> _position) {double x=1; cout <<"Error, call to illegal getPositionEnergy() in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

        vector<int> getLastModification() {vector<int> x; int TEMP=0; cout <<"Error, call to illegal getLastModification() in ligand.h. Press enter to continue"<<endl; cin>> TEMP; return x;}

        vector<int> chooseNextTargetPosition(ran& _ran) {vector<int> x; int TEMP=0; cout <<"Error, call to illegal chooseNextTargetPosition() in ligand.h. Press enter to continue"<<endl; cin>> TEMP; return x;}

        UInt chooseNextMutationIdentity(ran& _ran, vector <int> _position) {UInt x=1; cout <<"Error, call to illegal chooseNextMutationIdentity() in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

	void acceptModification() {}
        void rejectModification() {}
        int modify(ran& _ran) {int x=1; cout <<"Error, call to illegal modify(ran& _ran) in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

        int modify(ran& _ran, vector <int> _position) {int x=1; cout <<"Error, call to illegal modify(ran,vector) in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

        int mutate(vector <int> _position, UInt _resType) {int x=1; cout <<"Error, call to illegal mutate() in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

        void setupSystem(ran& _ran) {}
        void saveState(string& _filename) {}
        void resetAllBuffers() {}
        double getVolume(UInt _method) {double x=1; cout <<"Error, call to illegal getVolume() in ligand.h. Press enter to continue"<<endl; cin>> x; return x;}

private:
	
        // Internal Variable Declarations
	static UInt howMany;
	vector<atom*> itsAtoms;
	string itsNameString;
        string itsChainID;
        int itsLigTemplateType;
        bool hydrogensOn;
        static ligandTemplate itsLigTemplate;
        static bool dataBaseBuilt;
        
        // Energy modeling variables
	vector<UInt> itsAtomEnergyTypeIndex;
	static int itsCurrentEnergyType;  //IS THIS USED????
	static pmf itsPMF;
	//static amberElec itsAmberElec;
	static amberVDW itsAmberVDW;
	static solvation itsSolvation;
        static double cutoffDistance;
        static double cutoffDistanceSquared;
        
        // Intra-Ligand Connectivity Variables
        bool isConnected;
        vector<UInt> itsHeadNodes;
        vector<vector<UInt> > headNodeConnectVec;
        
        // Inter-Ligand (Ligand-Ligand) Connectivity Variables
        vector<ligand*> symmetryLinked;
};
#endif

