#include <petscksp.h>
#include "obs.h"
#undef __FUNCT__
#define __FUNCT__ "obs"

void cMasterObservables::initialize(){
//    DIM2 = 2*(N+1)*(Q+1);

}


PetscErrorCode cMasterObservables::ReshapeRho(){

//  ierr = MatCreate(PETSC_COMM_WORLD,&RhoMat);CHKERRQ(ierr);
//  ierr = MatSetSizes(RhoMat,PETSC_DECIDE,PETSC_DECIDE,DIM2,DIM2);CHKERRQ(ierr);
//  ierr = MatSetFromOptions(RhoMat);CHKERRQ(ierr);
//  ierr = MatSetUp(RhoMat);CHKERRQ(ierr);
//  /*
//         Currently, all PETSc parallel matrix formats are partitioned by
//         contiguous chunks of rows across the processors.  Determine which
//         rows of the matrix are locally owned.
//      */
//  ierr = MatGetOwnershipRange(RhoMat,&rhostart,&rhoend);CHKERRQ(ierr);
//
//  for(ROW = rhostart;ROW < rhostart; ROW++){
//
//  }

}
PetscErrorCode cMasterObservables::destruction(){
  /*
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
  */
//	ierr = MatDestroy(&RhoMat);CHKERRQ(ierr);


}
