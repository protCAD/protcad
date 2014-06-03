// filename: treeNode.cpp
// contents: class freeNode implmentation

#include "treeNode.h"
#include "atom.h"
//#include "CMath.h"

treeNode::treeNode()
{	
#ifdef __DEBUG
	cout << "treeNode constructor called" << endl;
#endif
	pItsFirstChild = 0;
	pItsNextSib = 0;
	pItsParent = 0;
	pItsPreviousSib = 0;
	parentSet = 0;
	firstChildSet = 0;
	nextSibSet = 0;
	prevSibSet = 0;
}

treeNode::~treeNode()
{	
#ifdef __DEBUG
	cout << "treeNode destructor called" << endl;
#endif

	if (pItsParent)
	{	pItsParent->setChild(pItsNextSib);
	}
	if (pItsFirstChild)
	{	pItsFirstChild->setParent(0);
	}
	if (pItsPreviousSib)
	{	pItsPreviousSib->setNextSib(pItsNextSib);
	}
	if (pItsNextSib)
	{	pItsNextSib->setPreviousSib(pItsPreviousSib);
	}
}

void treeNode::setChild(treeNode* _child)
{	pItsFirstChild = _child;
	firstChildSet = 1;
	if (_child)
	{	if (!pItsFirstChild->parentSet)
		{	pItsFirstChild->setParent(this);
		}
	}
	firstChildSet = 0;
}

void treeNode::setParent(treeNode* _parent)
{	pItsParent = _parent;
	parentSet = 1;
	if (_parent)
	{	if (!pItsParent->firstChildSet)
		{	pItsParent->setChild(this);
		}
	}
	parentSet = 0;
}

void treeNode::setNextSib(treeNode* _next)
{	pItsNextSib = _next;
	nextSibSet = 1;
	if (_next)
	{	if (!pItsNextSib->prevSibSet)
		{	pItsNextSib->setPreviousSib(this);
		}
	}
	nextSibSet = 0;
}

void treeNode::setPreviousSib(treeNode* _previous)
{	pItsPreviousSib = _previous;
	prevSibSet = 1;
	if (_previous)
	{	if (!pItsPreviousSib->prevSibSet)
		{	pItsPreviousSib->setNextSib(this);
		}
	}
	prevSibSet = 0;
}

vector<treeNode*> treeNode::getImmediateChildren()
{	treeNode* temp = getChild();
	vector<treeNode*> tempVec;
	tempVec.resize(0);

	if( temp )
	{	tempVec.push_back(temp);
		while( (temp = temp->getNextSib()) )
		{	tempVec.push_back(temp);
		}
	}

	return tempVec;
}

void treeNode::printImmediateChildren()
{
    vector<treeNode*> tempVec=getImmediateChildren();
    cout << "# children=" << tempVec.size() << "...";
    
    for(UInt i=0; i<tempVec.size(); i++)
    {
        cout << static_cast<atom*>(tempVec[i])->getName()<< " ";
    }
}

void treeNode::queryChildren()
{	treeNode* temp = getChild();
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		cout << static_cast<atom*>(temp)->getName()<< " ";
		temp->queryChildren();
		temp = temp->getNextSib();
	}
}

void treeNode::queryChildrensCoords()
{	treeNode* temp = getChild();
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		cout << static_cast<atom*>(temp)->getName() << "  ";
		cout << static_cast<atom*>(temp)->getCoords() << endl;
		temp->queryChildrensCoords();
		temp = temp->getNextSib();
	}
}

UInt treeNode::getNumChildren()
{	treeNode* temp = getChild();
	UInt counter = 0;
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		counter += temp->getNumChildren(counter);
		temp = temp->getNextSib();
	}
	return counter;
}

UInt treeNode::getNumChildren(UInt _counter)
{	treeNode* temp = getChild();
	_counter +=1;
	UInt tempcounter = 0;
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		_counter += temp->getNumChildren(tempcounter);
		temp = temp->getNextSib();
	}
	return _counter;
}

void treeNode::translateChildren(const dblVec& _dblVec)
{	treeNode* temp = getChild();
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		static_cast<atom*>(temp)->translate(_dblVec);
		temp->translateChildren(_dblVec);
		temp = temp->getNextSib();
	}
}

void treeNode::transformChildren(const dblMat& _matrix)
{	treeNode* temp = getChild();
	while (temp)
	{	// if we're going to do something with the children,
		// we have to do it here, before the recursive call
		// to collectChildren
		static_cast<atom*>(temp)->transform(_matrix);
		temp->transformChildren(_matrix);
		temp = temp->getNextSib();
	}
}

bool treeNode::isHeadNode()
{
	 if (!getParent() && !getPreviousSib() &&
		( getNextSib() || getChild())     )
		return true;
	 return false;
}
