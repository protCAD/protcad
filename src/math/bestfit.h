// these are linux flags
// FLAGS = -lg2c -lc -lm   //well, the lm is obvious

extern "C" 
{	
	void bestfit_(double *coords1, int *nat1, double *coords2, int *nat2,
		int *nat, double *coords3, int *list1, int* list2, double *rmsd, int*
		ierr, double *rotation_matrix, double *xc1, double *xc2, double *rmsdat);
}
