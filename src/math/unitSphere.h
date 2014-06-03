#include <vector>
#include "point.h"
#include "typedef.h"
#include <string>
#include <iostream>


#ifndef UNITSPHERE_H

#define UNITSPHERE_H

//*** File this one under the catagory of Stupid Code Tricks.
//*** All I wanted to do was just display a bounding sphere in my
//*** engine as a debugging and visualization aid.  Then I found out
//*** how much code is required just to tesselate a sphere.  Here
//*** is a quick and dirty class which will allow you to produce
//*** a tesselated sphere at some point in world space and with some
//*** radius.  The output tesselated sphere is a always a hard-coded
//*** 224 polygon mesh.
//***
//*** To use in your application simply inherte UnitSphereCallback and
//*** provide the virtual method SphereTriangle
//***
//*** To to produce the tesselated sphere simply invoke the static
//*** method:
//***
//*** UnitSphere::Tesselate(*this,center,radius);
//***
//*** You can replace SpherePoint with whatever your own 3d point
//*** representation might be.
//***
//*** Submitted to FlipCode.com Code of the day on October 10, 2000
//*** by John W. Ratcliff jratcliff@verant.com and is released into
//*** the public domain, for what it's worth.  If somebody knows of
//*** a lean mean high speed sphere tesselator please, feel free to
//*** modify this code snippet and replace it with the new implementation.



class spherePoint : public point
{
public:
	point* getSpherePoint(UInt _index);
	void setSphereSize(UInt _size);
	vector <point*> getAllSpherePoints();
	static UInt getSphereSize() { return itsSize; }	
	static double getX(UInt _index);
	static double getY(UInt _index);
	static double getZ(UInt _index);

private:
	static float gVertices[];
	static UInt itsSize;
};

#endif
