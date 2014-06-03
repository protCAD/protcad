// filename: point.h
// contents: class point is defined

#include "assert.h"
#include <vector>
#include "typedef.h"
//#include "CMath.h"

#ifndef POINT_H
#define POINT_H

class point
{
public:
	// Constructor and Destructor declaration
	point();
	point(const dblVec& _dblVec);
	point(const double _x, const double _y, const double _z);
	point(const point& _rhs);
	virtual ~point();

	// Overloaded operators
	point& operator= (const point& rhs);
	dblVec operator- (const point& rhs);

	// Coordinate Related Operations
	double getX() const {return itsCoords[0]; }
	double getY() const {return itsCoords[1]; }
	double getZ() const {return itsCoords[2]; }
	
	dblVec getCoords() const {return itsCoords; }
	
	void setCoords(const dblVec& _dblVec);
	void setCoords(const double _x, const double _y, const double _z);
	
	void transform(dblMat* _pDoubleMatrix);
	void transform(const dblMat& _dblMat);
	
	void translate(dblVec* _pDoubleVector);
	void translate(const dblVec& _dblVec);   
	
	// Distance Related Operations
	double distance(const point& otherPoint) const;
	double distance(point* pOtherPoint) const;
	
	// Static Variable Accessors
	static UInt getHowMany() {return howMany; }
	
protected:

	//non-static variable declaration
	dblVec itsCoords;

private:

	//static variable declaration
	static UInt howMany;

};
#endif
