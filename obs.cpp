#include <petscksp.h>
#include "obs.h"
#undef __FUNCT__
#define __FUNCT__ "obs"

PetscErrorCode cMasterObservables::photon(cMasterMatrix GMatrix){
	int nonzeros = 0;double phtn_n_r, phtn_fluc_r, tmpdiagrho;
	cout.precision(16);
	phtn_n_r = 0;phtn_fluc_r = 0;tmpdiagrho = 0;
	PetscScalar value0;
	for (GMatrix.ROW=GMatrix.rstart;GMatrix.ROW<GMatrix.rend;GMatrix.ROW++){
		GMatrix.block(GMatrix.ROW, GMatrix.r, GMatrix.m, GMatrix.n, GMatrix.p, GMatrix.q);
		if (GMatrix.r==0 || GMatrix.r==3){ // Getting diagonal elements of rho_up_up and rho_dn_dn for all photon and orbital numbers
			if (GMatrix.m==GMatrix.n && GMatrix.p==GMatrix.q){
				ierr_obs = VecGetValues(GMatrix.x,1,&(GMatrix.ROW),&value0);CHKERRQ(ierr_obs);
//				col[nonzeros] = m; // saved for photon number computation
//				cout << value[nonzeros] << '\t' << r+1 << '\t' << m << '\t' << p << endl; //
//				if (PetscImaginaryPart(value[nonzeros]) > 1.e-5) {
//					cerr << "check the convergence, the imaginary part is intolerably large, stopping now... " << endl;
//					exit(1);
//				}
//				tmpdiagrho += PetscRealPart(value0[nonzeros]);
				phtn_n_r += PetscRealPart(value0)*GMatrix.m; // checked the imaginary part is exceedingly small as it should be.
				phtn_fluc_r += PetscRealPart(value0)*GMatrix.m*GMatrix.m;
				nonzeros++;
			}
		}
	}
//	cout << "rank " << rank << " has photon number: " << phtn_n_r << " photon fluc: " << phtn_fluc_r << endl;
	MPI_Reduce(&phtn_n_r, &PhotonNumber, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&phtn_fluc_r, &PhotonFluc, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&tmpdiagrho, &tmpRhoDiagonal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	if (GMatrix.rank == 0) {
		PhotonFluc = (PhotonFluc - PhotonNumber*PhotonNumber)/PhotonNumber; // normalization to 1 for coherent state at root
	cout << "photon number is " << PhotonNumber  << '\t'
			<< "photon number fluctuation is " << PhotonFluc << endl;
//	cout << "sum of diag " << tmpRhoDiagonal << endl;
	}
	return ierr_obs;
}

PetscErrorCode cMasterObservables::oscillator(cMasterMatrix GMatrix){
	int nonzeros = 0;double phtn_n_r, phtn_fluc_r, tmpdiagrho;
	cout.precision(16);
	phtn_n_r = 0;phtn_fluc_r = 0;tmpdiagrho = 0;
	PetscScalar value0;
	for (GMatrix.ROW=GMatrix.rstart;GMatrix.ROW<GMatrix.rend;GMatrix.ROW++){
		GMatrix.block(GMatrix.ROW, GMatrix.r, GMatrix.m, GMatrix.n, GMatrix.p, GMatrix.q);
		if (GMatrix.r==0 || GMatrix.r==3){ // Getting diagonal elements of rho_up_up and rho_dn_dn for all photon and orbital numbers
			if (GMatrix.m==GMatrix.n && GMatrix.p==GMatrix.q){
				ierr_obs = VecGetValues(GMatrix.x,1,&(GMatrix.ROW),&value0);CHKERRQ(ierr_obs);
//				col[nonzeros] = p; // saved for average oscillator number computation
//				cout << value[nonzeros] << '\t' << r+1 << '\t' << m << '\t' << p << endl; //
//				if (PetscImaginaryPart(value[nonzeros]) > 1.e-5) {
//					cerr << "check the convergence, the imaginary part is intolerably large, stopping now... " << endl;
//					exit(1);
//				}
//				tmpdiagrho += PetscRealPart(value[nonzeros]);
				phtn_n_r += PetscRealPart(value0)*GMatrix.p; // checked the imaginary part is exceedingly small as it should be.
				phtn_fluc_r += PetscRealPart(value0)*GMatrix.p*GMatrix.p;
				nonzeros++;
			}
		}
	}
//	cout << "rank " << rank << " has photon number: " << phtn_n_r << " photon fluc: " << phtn_fluc_r << endl;
	MPI_Reduce(&phtn_n_r, &PhotonNumber, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&phtn_fluc_r, &PhotonFluc, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&tmpdiagrho, &tmpRhoDiagonal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	if (GMatrix.rank == 0) {
		PhotonFluc = (PhotonFluc - PhotonNumber*PhotonNumber)/PhotonNumber; // normalization to 1 for coherent state at root
	cout << "orbital number is " << PhotonNumber  << '\t'
			<< "orbital number fluctuation is " << PhotonFluc << endl;
//	cout << "sum of diag " << tmpRhoDiagonal << endl;
	}

//	for (int jtmp = 0; jtmp < size; ++jtmp) {
//		if (jtmp == rank) {
//			for (int itmp = 0; itmp < nonzeros; itmp++) {
////				if (col[itmp] !=0){
//					cout << "rank " << rank << " has " << value[itmp] << '\t' << "and photon number is " << col[itmp] << endl;
//			}
//		}
//	}
	return ierr_obs;
}

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
	return ierr_obs;
}
PetscErrorCode cMasterObservables::destruction(){
  /*
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
  */
//	ierr = MatDestroy(&RhoMat);CHKERRQ(ierr);

	return ierr_obs;
}
