#ifndef OBS_H_
#define OBS_H_
#include <petscksp.h>
#include <slepceps.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "steady.h"
class cMasterMatrix;
using namespace std;
class cMasterObservables{
private:
	EPS		eps;
	  EPSType		 type;
	PetscInt N,Q, DIM2, rhostart, rhoend,nlocalrow,ROW,nev,its,maxit,nconv;;
	Mat      RhoMat;
	PetscErrorCode ierr;
	double    	PhotonNumber, PhotonFluc, tmpRhoDiagonal;
	double    	OsciNumber, OsciFluc;
	double		ODTNumber_Q,ODTNumber_N;
	PetscReal	tol,error,re,im;
	PetscScalar kr2,ki;
	Vec            xr,xi;          /* RHS, test_exact solutions */
	  PetscMPIInt      rank, size;
		PetscInt m,n,p,q,r,rstart,rend;
public:
	cMasterObservables(){}
	~cMasterObservables(){}
	void initialize(cMasterMatrix);
	PetscErrorCode destruction();
	PetscErrorCode ReshapeRho(cMasterMatrix);
	PetscErrorCode photon(cMasterMatrix);
	PetscErrorCode oscillator(cMasterMatrix);
	PetscErrorCode negativity();
	PetscErrorCode checkODT(cMasterMatrix);
	PetscErrorCode spin_density(cMasterMatrix);
};
#endif
