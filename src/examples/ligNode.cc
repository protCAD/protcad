#include "ligNode.h"

UInt ligNode::howManyNodes=0;

ligNode::ligNode(UInt _ligdataIndex)
{
    cout << "ligNode constructor called" << endl;
    howManyNodes++;
    itsChildren.resize(0);
    pitsParent=0;
    parentSet=false;
    childSet=false;
    isIndependent=false;
    itsLigandDataIndex=_ligdataIndex;
}

ligNode::ligNode(UInt _ligdataIndex,bool _independence)
{
    cout << "ligNode constructor called for a headNode" << endl;
    howManyNodes++;
    itsChildren.resize(0);
    pitsParent=0;
    parentSet=false;
    childSet=false;
    isIndependent=_independence;
    itsLigandDataIndex=_ligdataIndex;
}

ligNode::~ligNode()
{
    if(itsChildren.size()>0)
    {
        cout << "ligNode Destructor: Cannot delete an atom with children." <<endl;
        exit(1);
    }
    howManyNodes--;
}


// Accessor Functions...
void ligNode::setChild(ligNode* _child)
{
    itsChildren.push_back(_child);
    childSet=true;
    _child->setParent(this);
    childSet=false;
}

void ligNode::setParent(ligNode* _parent)
{
    if(!parentSet)
    {
        pitsParent=_parent;
        parentSet=true;
    }
    
    else
    {
        cout<< "Error in ligNode::setParent(). This node already has a parent." << endl;
    }
    
    if(!_parent->childSet)
    {	
        _parent->setChild(this);
    }
}


// ***************************************************************************************
// ***************************************************************************************
//
// Start ligNodeRot Class
//
// ***************************************************************************************
// ***************************************************************************************

ligNodeRot::ligNodeRot(ligNode* _lig1, ligNode* _lig2)
{
    itsRotationMin=0;
    itsRotationMax=0;
    
    itsLinkedNodes.resize(0);
    
    itsLinkedNodes.push_back(_lig1);
    itsLinkedNodes.push_back(_lig2); //So, _lig1 will be itsLinkedNodes[0]

}

ligNodeRot::ligNodeRot(ligNode* _lig1, ligNode* _lig2, double _rotationLowBound, double _rotationHighBound)
{
    itsRotationMin=_rotationLowBound;
    itsRotationMax=_rotationHighBound;
    
    itsLinkedNodes.resize(0);
    
    itsLinkedNodes.push_back(_lig1);
    itsLinkedNodes.push_back(_lig2); //So, _lig1 will be itsLinkedNodes[0]

}

ligNodeRot::~ligNodeRot()
{
    howManyLigNodeRot--;
}

