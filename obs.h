#ifndef OBS_H_
#define OBS_H_
#include <petscksp.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "steady.h"
class cMasterMatrix;
using namespace std;
class cMasterObservables{
private:
	PetscInt DIM2, rhostart, rhoend;
	Mat      RhoMat;
	PetscErrorCode ierr_obs;
	double    	PhotonNumber, PhotonFluc, tmpRhoDiagonal;
public:
	cMasterObservables(){}
	~cMasterObservables(){}
	void initialize();
	PetscErrorCode destruction();
	PetscErrorCode ReshapeRho();
	PetscErrorCode photon(cMasterMatrix GMatrix);
	PetscErrorCode oscillator(cMasterMatrix GMatrix);
};
#endif
