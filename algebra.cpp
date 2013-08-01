#include "cholmod.h"
#include <cassert>
#include <ctype.h>
#include "global.h"
#include "util.h"
#include "vec.h"
//#include "umfpack.h"
#include "algebra.h"
#include "mpi.h"

// deliver the address of x
void Algebra::solve_CK(Matrix & A, cholmod_factor *L, cholmod_dense *&x, cholmod_dense *b, 
  cholmod_common *cm){
	cm->final_super = false;
	cm->final_asis = false;
	CK_decomp(A, L, cm);
	// then solve
	x = cholmod_solve(CHOLMOD_A, L, b, cm);
}

// doing cholesky decomposition
void Algebra::CK_decomp(Matrix &A, cholmod_factor *&L, cholmod_common *cm){
	cm->final_super = false;
	cm->final_asis = false;
	// doing factorization first
	cholmod_triplet * T;
	size_t n_row = A.get_row();
	size_t n_col = A.get_row();
	size_t nnz = A.size();
	if(nnz==0) return;
	int *Ti;
	int *Tj;
	double *Tx;
	int stype = -1;// lower triangular storage
	T = cholmod_allocate_triplet(n_row, n_col, nnz, stype, 
			CHOLMOD_REAL, cm);
	Ti = static_cast<int *>(T->i);
	Tj = static_cast<int *>(T->j);
	Tx = static_cast<double *>(T->x);
	// copy data into T
	for(size_t k=0;k<nnz;k++){
		Ti[k] = A.Ti[k];
		Tj[k] = A.Tj[k];
		Tx[k] = A.Tx[k];
	}
	T->nnz = nnz;
	cholmod_sparse * A_cholmod;
	A_cholmod = cholmod_triplet_to_sparse(T, nnz, cm);

	// free the triplet pointer
	cholmod_free_triplet(&T, cm);
	L = cholmod_analyze(A_cholmod, cm);
	cholmod_factorize(A_cholmod, L, cm);
	cholmod_free_sparse(&A_cholmod, cm);
}
