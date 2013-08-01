#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include "global.h"
#include "vec.h"
//#include "umfpack.h"
#include "cholmod.h"
class Vec;
class Algebra{
public:
	//static void solve(const Matrix & A, const Vec & b, Vec & x);
	static void solve_CK(Matrix & A, cholmod_factor *L, cholmod_dense *&x, cholmod_dense *b, cholmod_common *cm);
	static void CK_decomp(Matrix &A, cholmod_factor *&L, 
			cholmod_common *cm);
};

#endif
