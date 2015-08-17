#ifndef OBS_H_
#define OBS_H_
#include <petscksp.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "steady.h"
using namespace std;
class cMasterObservables{
private:
	PetscInt DIM2, rhostart, rhoend;
	Mat      RhoMat;

public:
	cMasterObservables(){}
	~cMasterObservables(){}
	void initialize();
	PetscErrorCode destruction();
	PetscErrorCode ReshapeRho();
};
#endif
