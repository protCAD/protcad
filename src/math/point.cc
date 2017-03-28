// filename: point.cpp
// contents: class point implementation

#include "point.h"

// static class member initialization

UInt point::howMany = 0;

// Constructors

point::point()
{	
	itsCoords.newsize(3);
	itsCoords[0] = 0.00;
	itsCoords[1] = 0.00;
	itsCoords[2] = 0.00;

#ifdef __POINT_DEBUG
	cout << "Point constructor called: "
	     << "point::point() "
	     << endl;
#endif
	howMany++;
}

point::point(const dblVec&  _dblVec)
{     	
	itsCoords.newsize(3);
#ifdef __POINT_DEBUG
	cout << "Point constructor called: " << endl
	     << "point::point(const dblVec&) "
	     << endl;
#endif
	setCoords(_dblVec);
	howMany++;
}			

point::point(const double _x, const double _y, const double _z)
{      
	itsCoords.newsize(3);
#ifdef __POINT_DEBUG
	cout << "Point constructor called: "
	     << "point::point(const double, const double, const double) "
	     << endl;
#endif
	setCoords(_x,_y,_z);
	howMany++;
}

// Copy constructor

point::point(const point& rhs)
{	
	itsCoords.newsize(3);
#ifdef __POINT_DEBUG
	cout << "Point copy constructor called " << endl;
#endif
	itsCoords = rhs.getCoords();
	//howMany++;
}

// Destructor
point::~point() 
{	
#ifdef __POINT_DEBUG
	cout << "Point destructor called" << endl;
#endif
	howMany--;
}

// Overloaded Operators

point& point::operator= (const point& rhs)
{
	if (this == &rhs)
		return *this;
	else
		itsCoords = rhs.getCoords();
	return *this;
}

dblVec point::operator- (const point& rhs)
{
	return itsCoords - rhs.itsCoords;
}

// Coordinate Related Operations

void point::setCoords(const dblVec& _dblVec)
{	
	double theSize = _dblVec.dim();
	if (theSize!=3)
	{	cout << "Error: setCoords with a vector of extent "
		     << theSize << " !" << endl;
		cout << "       Coordinates Untouched" << endl;
		cout << "Error reported by: point::setCoord(const dblVec&) "
		     << endl;
		return;
	};
        ASSERT(_dblVec[0] < 1e7 && _dblVec[0] > -1e7);
        ASSERT(_dblVec[1] < 1e7 && _dblVec[1] > -1e7);
        ASSERT(_dblVec[2] < 1e7 && _dblVec[2] > -1e7);
	itsCoords = _dblVec;
}

void point::setCoords(const double _x, const double _y, const double _z)
{	
	
        ASSERT(_x < 1e7 && _x > -1e7);
        ASSERT(_y < 1e7 && _y > -1e7);
        ASSERT(_z < 1e7 && _z > -1e7);
	itsCoords[0] = _x;
	itsCoords[1] = _y;
	itsCoords[2] = _z;
}

void point::transform(dblMat* _pDoubleMatrix)
{	dblVec newCoords = CMath::transform(itsCoords,_pDoubleMatrix);
	setCoords(newCoords);
}

void point::transform(const dblMat& _dblMat)
{	
	dblVec newCoords = CMath::transform(itsCoords,_dblMat);
	setCoords(newCoords);
}

void point::translate(dblVec* _pDoubleVector)
{	
	dblVec newCoords = CMath::translate(itsCoords,_pDoubleVector);
	setCoords(newCoords);
}

void point::translate(const dblVec& _dblVec)
{	
	dblVec newCoords = CMath::translate(itsCoords,_dblVec);
	setCoords(newCoords);
}	 

//  Distance Related Operations
 
double point::distance(point* pOtherPoint) const
{	return CMath::distance(itsCoords,pOtherPoint->getCoords());
}

double point::distance(const point& otherPoint) const
{	return CMath::distance(itsCoords,otherPoint.getCoords());
}
