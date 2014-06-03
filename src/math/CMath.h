// filename: CMath.h
// contents: all math related functions defined

#include "typedef.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "bestfit.h"

namespace CMath
{	
	double dotProduct(const dblVec& _vec1, const dblVec& _vec2);
	dblVec dotProduct(const dblMat& _mat1, const dblVec& _vec1);
//	dblMat dotProduct(const dblMat& _mat1, const dblMat& _mat2);
	// Non-constant functions
	dblVec translate(const dblVec* _pdv1,const dblVec* _pdv2);
	dblVec translate(const dblVec* _pdv1,const dblVec& _dv2);
	dblVec translate(const dblVec& _dv1,const dblVec* _pdv2);
	dblVec translate(const dblVec& _dv1,const dblVec& _dv2);
	
	dblVec transform(const dblVec* _pdv,const dblMat* _pdm);
	dblVec transform(const dblVec* _pdv,const dblMat& _dm);
	dblVec transform(const dblVec& _dv, const dblMat* _pdm);
	dblVec transform(const dblVec& _dv,const dblMat& _dm);

	// Constant functions
	double distance(const dblVec* _pdv1, const dblVec* _pdv2);
	double distance(const dblVec& _dv1, const dblVec* _pdv2);
	double distance(const dblVec* _pdv1, const dblVec& _dv2);
	double distance(const dblVec& _dv1, const dblVec& _dv2);
	double distanceSquared(const dblVec& _dv1, const dblVec& _dv2);

	double dihedral(const dblVec& _one, const dblVec& _two,
			const dblVec& _three, const dblVec& _four) ;
	double angle(const dblVec& _one, const dblVec& _two,
			const dblVec& _three);
	
	dblMat rotationMatrix(const dblVec& _dblVec, const double& _theta);

	dblVec centroid(const vector<dblVec>* thePoints, const vector<double>* theWeights); 
	dblVec cross(const dblVec& _one, const dblVec& _two) ;
	double determinant(const dblMat& _matrix);
	double determinant(const dblMat* _pMatrix); 
	double linearInterpolate(const double _m, const double _b,
		const double x);
/*
	double rmsd(vector<double>& _theWeights, 
		vector<dblVec>& _coord_fixed, vector, <dblVec>& _coord_mobile,
		vector<int>& _list1, vector<int>& _list2,
		vector<dblVec>& _coord_mobile_new, dblMat& _rotmat,
		dblVec& _centroid1, dblVec& _centroid2);
*/
}
