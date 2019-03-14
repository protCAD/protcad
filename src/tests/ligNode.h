#include "assert.h"
#include "typedef.h"

#ifndef LIGNODE_H
#define LIGNODE_H

class ligNode
{
public:
	// Constructors and Destructors
        ligNode(UInt _ligdataIndex);
        ligNode(UInt _ligdataIndex,bool _independence);
	~ligNode(); 
		
        // Accessors        
        void setChild(ligNode* _child);
        vector<ligNode*> getChildren(){return itsChildren;};
        UInt getNumChildren(){return itsChildren.size();}

        void setParent(ligNode* _parent);
        ligNode* getParent(){return pitsParent;};
                
        UInt getLigDataIndex(){return itsLigandDataIndex;};
        void setLigDataIndex(UInt _index){itsLigandDataIndex=_index;};
            //setLigDataIndex() should not be used unless ligandData has been
            //given the ability to dynamically rearrange the atom order.
            //ie: You have implemented a command that can change atoms
            //in the database (add, delete).
        
        void setIndependent(bool _status){isIndependent=_status;}
        bool getIndependence(){return isIndependent;};
    
private:
        vector<ligNode*> itsChildren;
        ligNode* pitsParent;
        bool isIndependent; //default is false unless set to true
        UInt itsLigandDataIndex; //what atom number is it in the ligand data class?
        
        bool parentSet;
        bool childSet;
                
public:
        // Global Variables
        static UInt howManyNodes;
    
};
#endif

#ifndef LIGNODEROT_H
#define LIGNODEROT_H

class ligNodeRot
{
public:
        // Constructors and Destructors
        ligNodeRot(ligNode* _lig1, ligNode* _lig2);
        ligNodeRot(ligNode* _lig1, ligNode* _lig2, double _rotationLowBound, double _rotationHighBound);
        ~ligNodeRot();
        
        //Accessors
        void setMinRot(double _rotation){itsRotationMin=_rotation;};
        void setMaxRot(double _rotation){itsRotationMax=_rotation;};
        void setRot(double _rotationMin, double _rotationMax){itsRotationMin=_rotationMin;itsRotationMax=_rotationMax;};
        double getMinRotBound(){return itsRotationMin;};
        double getMaxRotBound(){return itsRotationMax;};
        
private:
        vector<ligNode*> itsLinkedNodes;
        double itsRotationMin;
        double itsRotationMax;
        static UInt howManyLigNodeRot;

};
#endif



