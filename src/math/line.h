// filename: line.h
// contents: class line is defined

#include "assert.h"
#include <vector>
#include "typedef.h"
#include "point.h"

#ifndef LINE_H
#define LINE_H

class line
{
public:
	// Constructor and Destructor declaration
	line();
	line(const point& _point1, const point& _point2);
	line(const line& _rhs);
	virtual ~line();

	// Coordinate Related Operations
	point getPoint(UInt _index) const {return itsPoints[_index]; }
	void  setPoint(UInt _index, point& _point)
		{itsPoints[_index] = _point;}
	
	void transform(dblMat* _pDoubleMatrix);
	void transform(const dblMat& _dblMat);
	
	void translate(dblVec* _pDoubleVector);
	void translate(const dblVec& _dblVec);   
	
	// Distance Related Operations
	virtual double distance(const line& otherLine) const;
	virtual double distance(line* pOtherLine) const;
//	virtual double distance(const point& somePoint) const;
//	virtual double distance(point* pSomePoint) const;
	
	// Static Variable Accessors
	static UInt getHowMany() {return howMany; }

protected:
	virtual double calculateDistanceOfClosestApproach(point P0,
		point P1, point Q0, point Q1) const;
	
protected:

	//non-static variable declaration
	vector<point> itsPoints;

	//static variable declaration
	static UInt howMany;

};
#endif
