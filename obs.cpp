#include <petscksp.h>
#include "obs.h"
#undef __FUNCT__
#define __FUNCT__ "obs"

void cMasterObservables::initialize(cMasterMatrix GMatrix){
	N=GMatrix.N;
	Q=GMatrix.Q;
	DIM2 = 2*(N+1)*(GMatrix.Q+1);
	rstart=GMatrix.rstart;
	rend=GMatrix.rend;
	rank=GMatrix.rank;
	size=GMatrix.size;
}

PetscErrorCode cMasterObservables::photon(cMasterMatrix GMatrix){
	int nonzeros = 0;double phtn_n_r, phtn_fluc_r, tmpdiagrho;
	cout.precision(16);
	phtn_n_r = 0;phtn_fluc_r = 0;tmpdiagrho = 0;
	PetscScalar value0;
	for (ROW=rstart;ROW<rend;ROW++){
		GMatrix.block(ROW, r, m, n, p, q);
		if (r==0 || r==3){ // Getting diagonal elements of rho_up_up and rho_dn_dn for all photon and orbital numbers
			if (m==n && p==q){
				ierr = VecGetValues(GMatrix.x,1,&ROW,&value0);CHKERRQ(ierr);
//				col[nonzeros] = m; // saved for photon number computation
//				cout << value0 << '\t' << r+1 << '\t' << m << '\t' << p << endl; //
//				cout << value0 << '\t' << r*(N+1)*(Q+1)+m*(Q+1)+p+1 << endl;
//				if (PetscImaginaryPart(value[nonzeros]) > 1.e-5) {
//					cerr << "check the convergence, the imaginary part is intolerably large, stopping now... " << endl;
//					exit(1);
//				}
//				tmpdiagrho += PetscRealPart(value0[nonzeros]);
				phtn_n_r += PetscRealPart(value0)*m; // checked the imaginary part is exceedingly small as it should be.
				phtn_fluc_r += PetscRealPart(value0)*m*m;
				nonzeros++;
			}
		}
	}
//	cout << "rank " << rank << " has photon number: " << phtn_n_r << " photon fluc: " << phtn_fluc_r << endl;
	MPI_Reduce(&phtn_n_r, &PhotonNumber, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&phtn_fluc_r, &PhotonFluc, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&tmpdiagrho, &tmpRhoDiagonal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	if (rank == 0) {
		PhotonFluc = (PhotonFluc - PhotonNumber*PhotonNumber)/PhotonNumber; // normalization to 1 for coherent state at root
	cout << "photon number is " << PhotonNumber  << '\t'
			<< "photon number fluctuation is " << PhotonFluc << endl;
//	cout << "sum of diag " << tmpRhoDiagonal << endl;
	}
	return ierr;
}

PetscErrorCode cMasterObservables::oscillator(cMasterMatrix GMatrix){
	int nonzeros = 0;double phtn_n_r, phtn_fluc_r, tmpdiagrho;
	cout.precision(16);
	phtn_n_r = 0;phtn_fluc_r = 0;tmpdiagrho = 0;
	PetscScalar value0;
	for (ROW=rstart;ROW<rend;ROW++){
		GMatrix.block(ROW, r, m, n, p, q);
		if (r==0 || r==3){ // Getting diagonal elements of rho_up_up and rho_dn_dn for all photon and orbital numbers
			if (m==n && p==q){
				ierr = VecGetValues(GMatrix.x,1,&ROW,&value0);CHKERRQ(ierr);
//				col[nonzeros] = p; // saved for average oscillator number computation
//				cout << value[nonzeros] << '\t' << r+1 << '\t' << m << '\t' << p << endl; //
//				if (PetscImaginaryPart(value[nonzeros]) > 1.e-5) {
//					cerr << "check the convergence, the imaginary part is intolerably large, stopping now... " << endl;
//					exit(1);
//				}
//				tmpdiagrho += PetscRealPart(value[nonzeros]);
				phtn_n_r += PetscRealPart(value0)*p; // checked the imaginary part is exceedingly small as it should be.
				phtn_fluc_r += PetscRealPart(value0)*p*p;
				nonzeros++;
			}
		}
	}
//	cout << "rank " << rank << " has photon number: " << phtn_n_r << " photon fluc: " << phtn_fluc_r << endl;
	MPI_Reduce(&phtn_n_r, &PhotonNumber, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&phtn_fluc_r, &PhotonFluc, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(&tmpdiagrho, &tmpRhoDiagonal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	if (rank == 0) {
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
	return ierr;
}

PetscErrorCode cMasterObservables::ReshapeRho(cMasterMatrix GMatrix){


	int ir, ic, kr, kc;
	ierr = MatCreate(PETSC_COMM_WORLD,&RhoMat);CHKERRQ(ierr);
	ierr = MatSetType(RhoMat,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(RhoMat,PETSC_DECIDE,PETSC_DECIDE,DIM2,DIM2);CHKERRQ(ierr);
	ierr = MatSetFromOptions(RhoMat);CHKERRQ(ierr);
	ierr = MatSetUp(RhoMat);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(RhoMat,&rhostart,&rhoend);CHKERRQ(ierr);
	ierr = MatGetLocalSize(RhoMat,&nlocalrow, NULL);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(RhoMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(RhoMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//	ierr = MatSetOption(RhoMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
	PetscScalar value;
	// --> debug purpose of serial output:
	for (int ig = 0; ig < size; ++ig) {
		if (rank==ig) {
			for (ROW=rstart;ROW<rend;ROW++){
				GMatrix.block(ROW, r, m, n, p, q);
				ierr = VecGetValues(GMatrix.x,1,&ROW,&value);CHKERRQ(ierr);
//				cout << value << endl;
			}
		}
	}
	for (ROW=rstart;ROW<rend;ROW++){
		GMatrix.block(ROW, r, m, n, p, q);
		ierr = VecGetValues(GMatrix.x,1,&ROW,&value);CHKERRQ(ierr);
		if (r==0) {
			ir = 0; ic = 0;
		} else if (r==1) {
			ir = 0; ic = 1;
		} else if (r==2) {
			ir = 1; ic = 0;
		} else {
			ir = 1; ic = 1;
		}

		// TODO: ATTENTION ATTENTION ATTENTION please: This is rho matrix:
//		kr = ir*(N+1)*(Q+1)+m*(Q+1)+p;
//		kc = ic*(N+1)*(Q+1)+n*(Q+1)+q;
		// I don't have to use MPI_Send and MPI_Recv to manually control the communications, and petsc will do it for me.
//		ierr = MPI_Send(&(PetscRealPart(value)), 1, MPI_DOUBLE, int(floor(kr/floor(DIM2/GMatrix.size))), 42,MPI_COMM_WORLD);CHKERRQ(ierr);
//		ierr = MPI_Recv(&ADDRESS_OF_RhoMat_kr_kc, 1, MPI_DOUBLE, GMatrix.rank, 42,MPI_COMM_WORLD);CHKERRQ(ierr);
		//		ierr = MatSetValues(RhoMat,1,&kr,1,&kc,&value,INSERT_VALUES);CHKERRQ(ierr); // <-- rho matrix

		// TODO: ATTENTION ATTENTION ATTENTION please: This is partial transpose of rho matrix:
		kr = ir*(N+1)*(Q+1)+n*(Q+1)+q;
		kc = ic*(N+1)*(Q+1)+m*(Q+1)+p;
		// https://en.wikipedia.org/wiki/Negativity_%28quantum_mechanics%29
		// https://en.wikipedia.org/wiki/Peres%E2%80%93Horodecki_criterion
		ierr = MatSetValues(RhoMat,1,&kr,1,&kc,&value,INSERT_VALUES);CHKERRQ(ierr); // // <-- partial transpose of rho matrix

	}
	ierr = MatAssemblyBegin(RhoMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(RhoMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	/* // --> it would be great if I could check the hermicity of the density matrix. but the matrix format (MATMPIAIJ or MATMPIDENSE) is not yet supported. oops...
	//	PetscBool  *flg_MatIsHermitian;
	//	ierr = MatIsHermitian(RhoMat,1e-9,flg_MatIsHermitian);CHKERRQ(ierr); // a tolerance that I can accept.
	//	cout << "flg_MatIsHermitian value is " << *flg_MatIsHermitian << endl;
	 *
	 */
	/* Well... This is perhaps too strict in checking... no options to set the tolerance level...
	Mat RhoMat_Hermit;
	PetscBool  flg_MatIsHermitian;
	ierr = MatHermitianTranspose(RhoMat,MAT_INITIAL_MATRIX, &RhoMat_Hermit);CHKERRQ(ierr);
	ierr = MatEqual(RhoMat, RhoMat_Hermit, &flg_MatIsHermitian);CHKERRQ(ierr);
	if (!flg_MatIsHermitian) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"RhoMat is not Hermitian");
	ierr = MatDestroy(&RhoMat_Hermit);CHKERRQ(ierr);
	*/
	// Let me try the third route: doing a subtraction.
//	Mat RhoMat_Hermit;
//	ierr = MatHermitianTranspose(RhoMat,MAT_INITIAL_MATRIX, &RhoMat_Hermit);CHKERRQ(ierr);
//	ierr = MatAXPY(RhoMat,-1.0,RhoMat_Hermit,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	//	OKAY: it is hermitian to a degree, say 1e-16... TODO: need to remember to turn the checking off if you still want to proceed from here.
//	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_DENSE  );CHKERRQ(ierr);
//  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_MATLAB  );CHKERRQ(ierr);
//  ierr = MatView(RhoMat,	PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr);
	return ierr;
}

PetscErrorCode cMasterObservables::negativity(){
	ierr = MatCreateVecs(RhoMat,NULL,&xr);CHKERRQ(ierr);
	ierr = MatCreateVecs(RhoMat,NULL,&xi);CHKERRQ(ierr);
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps,RhoMat,NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr); // <-- ATTENTION: partial transpose of density matrix generally is NOT hermitian matrix.
	PetscBool her;
	ierr = EPSIsHermitian(eps, &her); CHKERRQ(ierr); if (rank==0) {cout << endl;cout << "is hermitian? " << her << endl;}
	PetscBool pos;
	ierr = EPSIsPositive(eps, &pos); CHKERRQ(ierr); if (rank==0) {cout << "is positive? " << pos << endl;cout << endl;}
	ierr = EPSSetDimensions(eps,DIM2,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr); //  <-- I am computing all the eigenvalues and let petsc decide the rest.
//	ierr = EPSSetType(eps,EPSPOWER);CHKERRQ(ierr);
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
//	if (nconv>0) {
//		/*
//		Display eigenvalues and relative errors
//		*/
//		ierr = PetscPrintf(PETSC_COMM_WORLD,
//		"           k          ||Ax-kx||/||kx||\n"
//		"   ----------------- ------------------\n");CHKERRQ(ierr);
//		for (int i=0;i<nconv;i++) {
//			ierr = EPSGetEigenpair(eps,i,&kr2,&ki,xr,xi);CHKERRQ(ierr);
//			/*
//			Compute the relative error associated to each eigenpair
//			*/
//			ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
//			#if defined(PETSC_USE_COMPLEX)
//				re = PetscRealPart(kr2);
//				im = PetscImaginaryPart(kr2);
//			#else
//				re = kr2;
//				im = ki;
//			#endif
//			if (im!=0.0) {
//				ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9f j %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
//			} else {
//				ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error);CHKERRQ(ierr);
//			}
////					cout << i <<"th eigenvalue is" << re + im*PETSC_i << endl;
//		}
//		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
//	}
	double negativity = 0.0;
		if (nconv>0) {
			for (int i=0;i<nconv;i++) {
				ierr = EPSGetEigenpair(eps,i,&kr2,&ki,xr,xi);CHKERRQ(ierr);
				re = PetscRealPart(kr2);
				negativity += (abs((double)re)-(double)re)/2.0;
			}
//			ierr = PetscPrintf(PETSC_COMM_WORLD,"negativity is %12f\n", negativity); CHKERRQ(ierr);
			if (rank==0) cout << "negativity is " << negativity << endl;
		} else{
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"None of the eigenparis are converged. Stopping now.");
		}
	return ierr;
}

PetscErrorCode cMasterObservables::destruction(){
  /*
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
  */
	ierr = MatDestroy(&RhoMat);CHKERRQ(ierr);
	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = VecDestroy(&xr);CHKERRQ(ierr);
	ierr = VecDestroy(&xi);CHKERRQ(ierr);
	return ierr;
}
