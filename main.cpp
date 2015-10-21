
static char help[] = "Solves the steady state of Master Equation.\n\n";

/*T
   Concepts: KSP^basic parallel example; Sparse matrix construction;
   Processors: n
   mpirun -n 6 SteadyState -ksp_monitor_short   -pc_type jacobi   -ksp_type gmres -ksp_gmres_restart 200
   or
   ./SteadyState -ksp_monitor_short   -pc_type jacobi   -ksp_type gmres -ksp_gmres_restart 200
T*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/
#include <petscksp.h>
#include <slepceps.h>
#include "steady.h"
#include "obs.h"
#undef __FUNCT__
#define __FUNCT__ "main"
#define root 0
int main(int argc,char **args){
  PetscErrorCode ierr;
  SlepcInitialize(&argc,&args,(char*)0,help);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     "======================================================================\n"
		     "The purpose of this program is to study the steady state solution of\n"
		     " quantum master equation, that governs the laser assisted SOC system.\n"
		     "Motivated, proposed, designed, implemented and researched \n"
		     "by Lin Dong at Rice University. \n"
		     "at "  __TIME__  ", on "  __DATE__  "\n"
		     "Petsc is initialized and program starts from \n"
		     __FILE__  "\n"
		     "======================================================================\n");CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Compute the matrix and right-hand-side vector that define
           the linear system, Gx = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cMasterMatrix GMatrix;
  GMatrix.initialize();
  GMatrix.construction();
  GMatrix.seek_steady_state();
  GMatrix.viewMatrix();
  cMasterObservables DensityOps;
  DensityOps.initialize(GMatrix);
  DensityOps.photon(GMatrix);
  DensityOps.oscillator(GMatrix);
  DensityOps.ReshapeRho(GMatrix);
  DensityOps.negativity();
  DensityOps.checkODT(GMatrix);
  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  DensityOps.destruction();
  GMatrix.destruction();
  ierr = SlepcFinalize();
  return 0;
}
