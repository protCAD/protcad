// filename: treeNode.h
// contents: class treeNode defined

#include "assert.h"
#include <vector>
#include <iostream>
#include "typedef.h"

#ifndef LIGNODE_H
class ligNode;
#endif


#ifndef TREENODE_H
#define TREENODE_H

class treeNode
{
public:
// constructors and destructors
	treeNode();
	~treeNode();

// accessors

	void setChild(treeNode* _child);
	void setParent(treeNode* _parent);
	void setNextSib(treeNode* _sib);
	void setPreviousSib(treeNode* _sib);
	bool isHeadNode();

	vector<treeNode*> getImmediateChildren();
        void printImmediateChildren();
	void queryChildren();
	void queryChildrensCoords();
	UInt getNumChildren();
	UInt getNumChildren(UInt _counter);
	void translateChildren(const dblVec& _vec);
	void transformChildren(const dblMat& _mat);

	treeNode* getChild()           {return pItsFirstChild; }
	treeNode* getNextSib()         {return pItsNextSib; }
	treeNode* getParent()          {return pItsParent; }
	treeNode* getPreviousSib()     {return pItsPreviousSib; }

	UInt parentSet;
	UInt firstChildSet;
	UInt nextSibSet;
	UInt prevSibSet;

protected:
	treeNode* pItsFirstChild;
	treeNode* pItsNextSib;
	treeNode* pItsParent;
	treeNode* pItsPreviousSib;
};

#endif
