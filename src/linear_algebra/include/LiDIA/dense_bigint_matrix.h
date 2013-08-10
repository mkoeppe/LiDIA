// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_DENSE_BIGINT_MATRIX_H_GUARD_
#define LIDIA_DENSE_BIGINT_MATRIX_H_GUARD_



#ifndef LIDIA_PRIMELIST_SIZE
#define LIDIA_PRIMELIST_SIZE 300
#endif

#ifndef LIDIA_PRIMELIST_MAXSIZE
#define LIDIA_PRIMELIST_MAXSIZE 86416
#endif


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_HNF_KERNEL_H_GUARD_
# include	"LiDIA/hnf_kernel.h"
#endif

#ifndef LIDIA_DENSE_MATRIX_H_GUARD_
# include	"LiDIA/dense_matrix.h"
#endif
#ifndef LIDIA_DENSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/dense_ring_matrix.h"
#endif

#ifndef LIDIA_DENSE_BIGINT_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/dense_bigint_matrix_kernel.h"
#endif
#ifndef LIDIA_DENSE_BIGINT_MATRIX_MODULES_H_GUARD_
# include	"LiDIA/dense_bigint_matrix_modules.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template<>
class dense_matrix< bigint > : public dense_ring_matrix< bigint >
{

	friend class dense_bigint_matrix_kernel< dense_matrix < bigint > >;

private:

	const dense_bigint_matrix_kernel< dense_matrix < bigint > > dense_modul;

	hnf_kernel< bigint, row_oriented_dense_matrix_modules < bigint, MR < bigint > >,
		nf_conf1< bigint, row_oriented_dense_matrix_modules < bigint, MR < bigint > >, MR< bigint > >,
		MR< bigint > > hnf_modul;

	//
	// constructors
	//

public:

	dense_matrix() : dense_ring_matrix< bigint > () {}
	dense_matrix(lidia_size_t a, lidia_size_t b) : dense_ring_matrix< bigint > (a, b) {}
	dense_matrix(lidia_size_t a, lidia_size_t b, const bigint **v) :
		dense_ring_matrix< bigint > (a, b, v) {}
	dense_matrix(const MR< bigint > &A) : dense_ring_matrix< bigint > (A) {}
	dense_matrix(const dense_base_matrix< bigint > &A) : dense_ring_matrix< bigint > (A) {}

	//
	// destructor
	//

public:

	~dense_matrix() {}

	//
	// Casts
	//

public:

//   operator bigfloat_matrix()
//   {
//     bigfloat_matrix B(rows, columns);
//     for (register lidia_size_t i = 0; i < rows; i++)
//       for (register lidia_size_t j = 0; j < columns; j++)
// 	B.sto(i, j, (bigfloat)value[i][j]);
//     return B;
//   }


	//
	// pseudo-division
	//

public:

	friend void divide(dense_matrix< bigint > &, const dense_matrix< bigint > &, const bigint &);
	friend void compwise_divide(dense_matrix< bigint > &, const dense_matrix< bigint > &,
				    const dense_matrix< bigint > &);

protected:

	void divide(const dense_matrix< bigint > &, const bigint &);
	void compwise_divide(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	//
	// pseudo - division
	//

public:

	friend dense_matrix< bigint > operator / (const dense_matrix< bigint > &, const bigint &);
	friend dense_matrix< bigint > & operator /= (dense_matrix< bigint > &, const bigint &);

	//
	// remainder
	//

public:

	friend dense_matrix< bigint > operator % (const dense_matrix< bigint > &, const bigint &);

	friend dense_matrix< bigint > &operator %= (dense_matrix< bigint > &, const bigint &);

	friend void remainder(dense_matrix< bigint > &, const dense_matrix< bigint > &, const bigint &);

protected:

	void remainder(const dense_matrix< bigint > &, const bigint &);

	//
	// assignments
	//

public:

	dense_matrix< bigint > & operator = (const dense_base_matrix< bigint > & B)
	{
		dense_base_matrix< bigint >::assign(B);
		return *this;
	}

	//
	// norms and bounds
	//

public:

	void max(bigint &) const;
	bigint max() const
	{
		bigint MAX;

		max(MAX);
		return MAX;
	}
	friend bigint max(const dense_matrix< bigint > &);

	void max_abs(bigint &) const;
	bigint max_abs() const
	{
		bigint RES;

		max_abs(RES);
		return RES;
	}
	friend bigint max_abs(const dense_matrix< bigint > &);

	void max_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint max_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MAX;

		max_pos(MAX, x, y);
		return MAX;
	}
	friend bigint max_pos(const dense_matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void max_abs_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint max_abs_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint RES;

		max_abs_pos(RES, x, y);
		return RES;
	}
	friend bigint max_abs_pos(const dense_matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void min(bigint &) const;
	bigint min() const
	{
		bigint MIN;

		min(MIN);
		return MIN;
	}
	friend bigint min(const dense_matrix< bigint > &);

	void min_abs(bigint &) const;
	bigint min_abs() const
	{
		bigint MIN;

		min_abs(MIN);
		return MIN;
	}
	friend bigint min_abs(const dense_matrix< bigint > &);

	void min_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint min_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MIN;

		min_pos(MIN, x, y);
		return MIN;
	}
	friend bigint min_pos(const dense_matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void min_abs_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint min_abs_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MIN;

		min_abs_pos(MIN, x, y);
		return MIN;
	}
	friend bigint min_abs_pos(const dense_matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	bigint hadamard() const
	{
		bigint H;

		hadamard(H);
		return H;
	}
	friend bigint hadamard(const dense_matrix< bigint > &);

	void hadamard(bigint &) const;
	void hadamard2(bigint &) const;

	void row_norm(bigint &, lidia_size_t, long) const;
	bigint row_norm(lidia_size_t i, long n) const
	{
		bigint RES;

		row_norm(RES, i, n);
		return RES;
	}
	friend bigint row_norm(const dense_matrix< bigint > &, lidia_size_t, long);

	void column_norm(bigint &, lidia_size_t, long) const;
	bigint column_norm(lidia_size_t i, long n) const
	{
		bigint RES;

		column_norm(RES, i, n);
		return RES;
	}
	friend bigint column_norm(const dense_matrix< bigint > &, lidia_size_t, long);

	//
	// randomize
	//

public:

	void randomize(const bigint &);
	void randomize_with_det(const bigint &, const bigint &);

	//
	// Chinese remaindering theorem
	//

protected:

	friend dense_matrix< bigint > chinrest(const dense_matrix< bigint > *, const bigint *);
	friend void chinrest(dense_matrix< bigint > &, const dense_matrix< bigint > *, const bigint *);

	///////////////////////////
	// BEGIN: Linear algebra //
	// PART 1                //
	///////////////////////////

	//
	// rank
	//

public:

	lidia_size_t rank() const;
	friend lidia_size_t rank(const dense_matrix< bigint > &);

protected:

	lidia_size_t rank(const bigint &) const;

	//
	// rank and linearly independent rows
	//

public:

	lidia_size_t *lininr() const;
	friend lidia_size_t *lininr(const dense_matrix< bigint > &);

	void lininr(base_vector< lidia_size_t > &) const;
	friend void lininr(base_vector< lidia_size_t > &v, const dense_matrix< bigint > &);

public:

	lidia_size_t *lininr(const bigint &) const;
	lidia_size_t *lininr2(const bigint &) const;
	lidia_size_t *lininr3(const bigint &) const;

	//
	// rank and linearly independent columns
	//

public:

	lidia_size_t *lininc() const;
	friend lidia_size_t *lininc(const dense_matrix< bigint > &);

	void lininc(base_vector< lidia_size_t > &) const;
	friend void lininc(base_vector< lidia_size_t > &, const dense_matrix< bigint > &);

public:

	lidia_size_t *lininc(const bigint &) const;
	lidia_size_t *lininc2(const bigint &) const;
	lidia_size_t *lininc3(const bigint &) const;

	//
	// regular expansion
	//

public:

	void regexpansion(const lidia_size_t *);
	friend void regexpansion(dense_matrix< bigint > &, const lidia_size_t *);

	//
	// adjoint matrix
	//

public:

	void adj(const dense_matrix< bigint > &);
	void adj(const dense_matrix< bigint > &, const bigint &);
	friend dense_matrix< bigint > adj(const dense_matrix< bigint > &);

	//
	// lattice determinant
	//

	///////////////////////////////////
	// BEGIN: INTERFACE - latticedet //
	///////////////////////////////////

public:

	void latticedet(bigint &DET)
	{
		latticedet2(DET);
	}
	bigint latticedet()
	{
		bigint RET;

		latticedet2(RET);
		return RET;
	}

	friend bigint latticedet(const dense_matrix< bigint > &);

	/////////////////////////////////
	// END: INTERFACE - latticedet //
	/////////////////////////////////

public:

	void latticedet1(bigint &) const;
	bigint latticedet1() const
	{
		bigint DET;

		latticedet1(DET);
		return DET;
	}
	friend bigint latticedet1(const dense_matrix< bigint > &);

	void latticedet2(bigint &) const;
	bigint latticedet2() const
	{
		bigint DET;

		latticedet2(DET);
		return DET;
	}
	friend bigint latticedet2(const dense_matrix< bigint > &);

	void latticedet3(bigint &) const;
	bigint latticedet3() const
	{
		bigint DET;

		latticedet3(DET);
		return DET;
	}
	friend bigint latticedet3(const dense_matrix< bigint > &);

public:

	void latticedet1(bigint &, const bigint &) const;
	void latticedet2(bigint &, const bigint &) const;
	void latticedet3(bigint &, const bigint &) const;
	void latticedet4(bigint &, const bigint &) const;
	void latticedet5(bigint &, const bigint &) const;

	//
	// determinant
	//

public:

	void det(bigint &) const;
	bigint det() const
	{
		bigint DET;

		det(DET);
		return DET;
	}
	friend bigint det(const dense_matrix< bigint > &);

	void det(bigint &, const bigint &) const;
	void det(bigint &, const bigint &, int) const;

	//
	// characteristic polynomial
	//

public:

	bigint *charpoly() const;
	friend bigint *charpoly(const dense_matrix< bigint > &);

	void charpoly(base_vector< bigint > &) const;
	friend void charpoly(base_vector< bigint > &, const dense_matrix< bigint > &);

	/////////////////////////
	// END: Linear algebra //
	// PART 1              //
	/////////////////////////

	///////////////////////////
	// BEGIN: Linear algebra //
	// PART 2                //
	///////////////////////////

	//
	// Hermite normal form
	//

	////////////////////////////
	// BEGIN: INTERFACE - HNF //
	////////////////////////////

public:

	void hnf()
	{
		hnf_havas_cont();
	}
	void hnf(dense_matrix< bigint > &A)
	{
		hnf_havas_cont(A);
	}

	friend dense_matrix< bigint > hnf(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > hnf(const dense_matrix< bigint > &, dense_matrix< bigint > &);

	//////////////////////////
	// END: INTERFACE - HNF //
	//////////////////////////

public:

	void hnf_havas();
	void hnf_havas(dense_matrix< bigint > &);

	friend dense_matrix< bigint > hnf_havas(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > hnf_havas(const dense_matrix< bigint > &, dense_matrix< bigint > &);

	void hnf_havas_cont();
	void hnf_havas_cont(dense_matrix< bigint > &);

	friend dense_matrix< bigint > hnf_havas_cont(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > hnf_havas_cont(const dense_matrix< bigint > &, dense_matrix< bigint > &);

	void hnfmod_dkt(const bigint &);
	void hnfmod_dkt()
	{
		hnfmod_dkt(latticedet2());
	}

	friend dense_matrix< bigint > hnfmod_dkt(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > hnfmod_dkt(const dense_matrix< bigint > &, const bigint &);

	void hnfmod_cohen(const bigint &);
	void hnfmod_cohen()
	{
		hnfmod_cohen(latticedet2());
	}

	friend dense_matrix< bigint > hnfmod_cohen(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > hnfmod_cohen(const dense_matrix< bigint > &, const bigint &);

	void hnfmod_mueller(dense_matrix< bigint > &);

	friend dense_matrix< bigint > hnfmod_mueller(const dense_matrix< bigint > &, dense_matrix< bigint > &);

	//
	// Kernel
	//

	///////////////////////////////
	// BEGIN: INTERFACE - kernel //
	///////////////////////////////

public:

	void kernel(const dense_matrix< bigint > &B)
	{
		kernel1(B);
	}

	friend dense_matrix< bigint > kernel(const dense_matrix< bigint > &);

	/////////////////////////////
	// END: INTERFACE - kernel //
	/////////////////////////////

public:

	void kernel1(const dense_matrix< bigint > &);

	friend dense_matrix< bigint > kernel1(const dense_matrix< bigint > &);

	void kernel2(const dense_matrix< bigint > &);

	friend dense_matrix< bigint > kernel2(const dense_matrix< bigint > &);

	//
	// regular Invimage
	//

	////////////////////////////////////
	// BEGIN: INTERFACE - reginvimage //
	////////////////////////////////////

public:

	void reginvimage(const dense_matrix< bigint > &A, const dense_matrix< bigint > &B)
	{
		reginvimage1(A, B);
	}

	friend dense_matrix< bigint > reginvimage(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	//////////////////////////////////
	// END: INTERFACE - reginvimage //
	//////////////////////////////////

public:

	void reginvimage1(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	friend dense_matrix< bigint > reginvimage1(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	void reginvimage2(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	friend dense_matrix< bigint > reginvimage2(const dense_matrix< bigint > &, const dense_matrix< bigint > &);

	//
	// Image
	//

	//////////////////////////////
	// BEGIN: INTERFACE - image //
	//////////////////////////////

public:

	void image(const dense_matrix< bigint > &B)
	{
		image1(B);
	}

	friend dense_matrix< bigint > image(const dense_matrix< bigint > &);

	////////////////////////////
	// END: INTERFACE - image //
	////////////////////////////

public:

	void image1(const dense_matrix< bigint > &);

	friend dense_matrix< bigint > image1(const dense_matrix< bigint > &);

	void image2(const dense_matrix< bigint > &);

	friend dense_matrix< bigint > image2(const dense_matrix< bigint > &);

	//
	// InvImage
	//

public:

	void invimage(const dense_matrix< bigint > &, const bigint *);
	friend dense_matrix< bigint > invimage(const dense_matrix< bigint > &, const bigint *);

	void invimage(const dense_matrix< bigint > &, const math_vector< bigint > &);
	friend dense_matrix< bigint > invimage(const dense_matrix< bigint > &, const math_vector< bigint > &);

	//
	// solve
	//

public:

	void solve(const dense_matrix< bigint > &A, const bigint *b)
	{
		invimage(A, b);
	}
	friend dense_matrix< bigint > solve(const dense_matrix< bigint > &, const bigint *);

	void solve(const dense_matrix< bigint > &A, const math_vector< bigint > &b)
	{
		invimage(A, b);
	}

	friend dense_matrix< bigint > solve(const dense_matrix< bigint > &, const math_vector< bigint > &);

	//
	// Smith normal form
	//

	////////////////////////////
	// BEGIN: INTERFACE - SNF //
	////////////////////////////

public:

	void snf()
	{
		snf_havas();
	}
	void snf(dense_matrix< bigint > &A, dense_matrix< bigint > &B)
	{
		snf_havas(A, B);
	}

	friend dense_matrix< bigint > snf(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snf(const dense_matrix< bigint > &, dense_matrix< bigint > &, dense_matrix< bigint > &);

	//////////////////////////
	// END: INTERFACE - SNF //
	//////////////////////////

public:

	void snf_hartley();
	void snf_hartley(dense_matrix< bigint > &, dense_matrix< bigint > &);

	friend dense_matrix< bigint > snf_hartley(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snf_hartley(const dense_matrix< bigint > &, dense_matrix< bigint > &,
						  dense_matrix< bigint > &);

	void snf_simple();
	void snf_simple(dense_matrix< bigint > &, dense_matrix< bigint > &);

	friend dense_matrix< bigint > snf_simple(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snf_simple(const dense_matrix< bigint > &, dense_matrix< bigint > &,
						 dense_matrix< bigint > &);

	void snf_havas();
	void snf_havas(dense_matrix< bigint > &, dense_matrix< bigint > &);

	friend dense_matrix< bigint > snf_havas(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snf_havas(const dense_matrix< bigint > &, dense_matrix< bigint > &,
						dense_matrix< bigint > &);

	void snf_mult(long art = 1);
	void snf_mult(dense_matrix< bigint > &, dense_matrix< bigint > &, long art = 1);

	friend dense_matrix< bigint > snf_mult(const dense_matrix< bigint > &, long art = 1);
	friend dense_matrix< bigint > snf_mult(const dense_matrix< bigint > &, dense_matrix< bigint > &,
					       dense_matrix< bigint > &, long art = 1);

	void snf_add(long art = 1);
	void snf_add(dense_matrix< bigint > &, dense_matrix< bigint > &, long art = 1);

	friend dense_matrix< bigint > snf_add(const dense_matrix< bigint > &, long art = 1);
	friend dense_matrix< bigint > snf_add(const dense_matrix< bigint > &, dense_matrix< bigint > &,
					      dense_matrix< bigint > &, long art = 1);

	void snf_new(long art = 1);
	void snf_new(dense_matrix< bigint > &, dense_matrix< bigint > &, long art = 1);

	friend dense_matrix< bigint > snf_new(const dense_matrix< bigint > &, long art = 1);
	friend dense_matrix< bigint > snf_new(const dense_matrix< bigint > &, dense_matrix< bigint > &,
					      dense_matrix< bigint > &, long art = 1);

	void snfmod_dkt(const bigint &);
	void snfmod_dkt()
	{
		snfmod_dkt(latticedet2());
	}

	friend dense_matrix< bigint > snfmod_dkt(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snfmod_dkt(const dense_matrix< bigint > &, const bigint &);

	void snfmod_cohen(const bigint &);
	void snfmod_cohen()
	{
		snfmod_cohen(latticedet2());
	}

	friend dense_matrix< bigint > snfmod_cohen(const dense_matrix< bigint > &);
	friend dense_matrix< bigint > snfmod_cohen(const dense_matrix< bigint > &, const bigint &);

	/////////////////////////
	// END: Linear algebra //
	// PART 2              //
	/////////////////////////

protected:

	void gauss();

	//
	// mgcd Computation
	//

public: //protected:

	bigint *mgcd(const bigint *, lidia_size_t);

	friend bigint *mgcd1(const bigint *, lidia_size_t, dense_matrix< bigint > &);
	bigint *mgcd1(const bigint *, lidia_size_t);

	friend bigint *mgcd2(const bigint *, lidia_size_t);
	friend void mgcd2(bigint &, const bigint *, lidia_size_t);

	friend bigint *mgcd2(const bigint *, lidia_size_t, dense_matrix< bigint > &);
	bigint *mgcd2(const bigint *, lidia_size_t);

	friend bigint *mgcd3(const bigint *, lidia_size_t, dense_matrix< bigint > &);
	bigint *mgcd3(const bigint *, lidia_size_t);

	friend bigint *mgcd4(const bigint *, lidia_size_t, dense_matrix< bigint > &);
	bigint *mgcd4(const bigint *, lidia_size_t);
};

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// division
//

inline void
divide(dense_matrix< bigint > &RES, const dense_matrix< bigint > &M, const bigint &a)
{
	RES.divide(M, a);
}



inline void
compwise_divide(dense_matrix< bigint > &RES, const dense_matrix< bigint > &M,
		const dense_matrix< bigint > &N)
{
	RES.compwise_divide(M, N);
}



//
// remainder
//

inline void
remainder(dense_matrix< bigint > &A, const dense_matrix< bigint > &B, const bigint &mod)
{
	A.remainder(B, mod);
}



/////////////////////////////////
// END: arithmetic procedures  //
// BEGIN: arithmetic operators //
/////////////////////////////////

//
// division
//

inline dense_matrix< bigint >
operator / (const dense_matrix< bigint > &A, const bigint &m)
{
	dense_matrix< bigint > RES(A.rows, A.columns);
	RES.divide(A, m);
	return RES;
}



inline dense_matrix< bigint > &
operator /= (dense_matrix< bigint > &A, const bigint &m)
{
	A.divide(A, m);
	return A;
}



//
// remainder
//

inline dense_matrix< bigint >
operator % (const dense_matrix< bigint > &A, const bigint &mod)
{
	dense_matrix< bigint > RES(A.rows, A.columns);
	RES.remainder(A, mod);
	return RES;
}



inline dense_matrix< bigint > &
operator %= (dense_matrix< bigint > &A, const bigint &mod)
{
	A.remainder(A, mod);
	return A;
}



///////////////////////////////
// END: arithmetic operators //
///////////////////////////////

//
// norms and bounds
//

inline bigint
max(const dense_matrix< bigint > &A)
{
	return A.max();
}



inline bigint
max_abs(const dense_matrix< bigint > &A)
{
	return A.max_abs();
}



inline bigint
max_pos(const dense_matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.max_pos(x, y);
}



inline bigint
max_abs_pos(const dense_matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.max_abs_pos(x, y);
}



inline bigint
min(const dense_matrix< bigint > &A)
{
	return A.min();
}



inline bigint
min_abs(const dense_matrix< bigint > &A)
{
	return A.min_abs();
}



inline bigint
min_pos(const dense_matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.min_pos(x, y);
}



inline bigint
min_abs_pos(const dense_matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.min_abs_pos(x, y);
}



inline bigint
hadamard(const dense_matrix< bigint > &A)
{
	bigint H;
	A.hadamard(H);
	return H;
}



inline bigint
hadamard2(const dense_matrix< bigint > &A)
{
	bigint H;
	A.hadamard2(H);
	return H;
}



inline bigint
row_norm(const dense_matrix< bigint > &A, lidia_size_t i, long n)
{
	bigint RES;
	A.row_norm(RES, i, n);
	return RES;
}



inline bigint
column_norm(const dense_matrix< bigint > &A, lidia_size_t i, long n)
{
	bigint RES;
	A.column_norm(RES, i, n);
	return RES;
}



//
// Chinese remaindering theorem
//

inline dense_matrix< bigint >
chinrest(const dense_matrix< bigint > *v, const bigint *prim)
{
	dense_matrix< bigint > A(v[0].rows, v[0].columns);
	chinrest(A, v, prim);
	return A;
}



///////////////////////////
// BEGIN: Linear algebra //
// PART 1                //
///////////////////////////

//
// rank
//

inline lidia_size_t rank(const dense_matrix< bigint > &B)
{
	return B.rank();
}



//
// rank and linearly independent rows
//

inline lidia_size_t *
lininr(const dense_matrix< bigint > &B)
{
	return B.lininr();
}



inline void
lininr(base_vector< lidia_size_t > &v, const dense_matrix< bigint > &B)
{
	B.lininr(v);
}



//
// rank and linearly independent columns
//

inline lidia_size_t *
lininc(const dense_matrix< bigint > &B)
{
	return B.lininc();
}



inline void
lininc(base_vector< lidia_size_t > &v, const dense_matrix< bigint > &B)
{
	B.lininc(v);
}



//
// regular expansion
//

inline void
regexpansion(dense_matrix< bigint > &A, const lidia_size_t *v)
{
	A.regexpansion(v);
}



//
// adjoint matrix
//

inline dense_matrix< bigint >
adj(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B(A.rows, A.columns);
	B.adj(A);
	return B;
}



//
// lattice determinant
//

///////////////////////////////////
// BEGIN: INTERFACE - latticedet //
///////////////////////////////////

inline bigint latticedet(const dense_matrix< bigint > &A)
{
	bigint DET;
	A.latticedet2(DET);
	return DET;
}



/////////////////////////////////
// END: INTERFACE - latticedet //
/////////////////////////////////

inline bigint
latticedet1(const dense_matrix< bigint > &A)
{
	bigint DET;
	A.latticedet1(DET);
	return DET;
}



inline bigint
latticedet2(const dense_matrix< bigint > &A)
{
	bigint DET;
	A.latticedet2(DET);
	return DET;
}



inline bigint
latticedet3(const dense_matrix< bigint > &A)
{
	bigint DET;
	A.latticedet3(DET);
	return DET;
}



//
// determinant
//

inline bigint
det(const dense_matrix< bigint > &A)
{
	bigint DET;
	A.det(DET);
	return DET;
}



//
// characteristic polynomial
//

inline bigint *
charpoly(const dense_matrix< bigint > &A)
{
	return A.charpoly();
}



inline void
charpoly(base_vector< bigint > &v, const dense_matrix< bigint > &A)
{
	A.charpoly(v);
}



/////////////////////////
// END: Linear algebra //
// PART 1              //
/////////////////////////

///////////////////////////
// BEGIN: Linear algebra //
// PART 2                //
///////////////////////////

//
// Hermite normal form
//

////////////////////////////
// BEGIN: INTERFACE - HNF //
////////////////////////////

inline dense_matrix< bigint >
hnf(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas_cont();
	return B;
}



inline dense_matrix< bigint >
hnf(const dense_matrix< bigint > &A, dense_matrix< bigint > &TR)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas_cont(TR);
	return B;
}



//////////////////////////
// END: INTERFACE - HNF //
//////////////////////////

inline dense_matrix< bigint >
hnf_simple(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnf_simple();
	return B;
}



inline dense_matrix< bigint >
hnf_simple(const dense_matrix< bigint > &A, dense_matrix< bigint > &TRANS)
{
	dense_matrix< bigint > B = A;
	B.hnf_simple(TRANS);
	return B;
}



inline dense_matrix< bigint >
hnf_havas(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas();
	return B;
}



inline dense_matrix< bigint >
hnf_havas(const dense_matrix< bigint > &A, dense_matrix< bigint > &TRANS)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas(TRANS);
	return B;
}



inline dense_matrix< bigint >
hnf_havas_cont(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas_cont();
	return B;
}



inline dense_matrix< bigint >
hnf_havas_cont(const dense_matrix< bigint > &A, dense_matrix< bigint > &TRANS)
{
	dense_matrix< bigint > B = A;
	B.hnf_havas_cont(TRANS);
	return B;
}



inline dense_matrix< bigint >
hnfmod_dkt(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnfmod_dkt(B.latticedet2());
	return B;
}



inline dense_matrix< bigint >
hnfmod_dkt(const dense_matrix< bigint > &A, const bigint &h)
{
	dense_matrix< bigint > B = A;
	B.hnfmod_dkt(h);
	return B;
}



inline dense_matrix< bigint >
hnfmod_cohen(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.hnfmod_cohen(B.latticedet2());
	return B;
}



inline dense_matrix< bigint >
hnfmod_cohen(const dense_matrix< bigint > &A, const bigint &h)
{
	dense_matrix< bigint > B = A;
	B.hnfmod_cohen(h);
	return B;
}



inline dense_matrix< bigint >
hnfmod_mueller(const dense_matrix< bigint > &A, dense_matrix< bigint > &TRANS)
{
	dense_matrix< bigint > B = A;
	B.hnfmod_mueller(TRANS);
	return B;
}



//
// Kernel
//

///////////////////////////////
// BEGIN: INTERFACE - kernel //
///////////////////////////////

inline dense_matrix< bigint >
kernel(const dense_matrix< bigint > &B)
{
	dense_matrix< bigint > RET;
	RET.kernel1(B);
	return RET;
}



/////////////////////////////
// END: INTERFACE - kernel //
/////////////////////////////

inline dense_matrix< bigint >
kernel1(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > RET;
	RET.kernel(A);
	return RET;
}



inline dense_matrix< bigint >
kernel2(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > RET;
	RET.kernel2(A);
	return RET;
}



//
// regular Invimage
//

////////////////////////////////////
// BEGIN: INTERFACE - reginvimage //
////////////////////////////////////

inline dense_matrix< bigint >
reginvimage(const dense_matrix< bigint > &A, const dense_matrix< bigint > &B)
{
	dense_matrix< bigint > X;
	X.reginvimage(A, B);
	return X;
}



//////////////////////////////////
// END: INTERFACE - reginvimage //
//////////////////////////////////

inline dense_matrix< bigint >
reginvimage1(const dense_matrix< bigint > &A, const dense_matrix< bigint > &B)
{
	dense_matrix< bigint > X;
	X.reginvimage1(A, B);
	return X;
}



inline dense_matrix< bigint >
reginvimage2(const dense_matrix< bigint > &A, const dense_matrix< bigint > &B)
{
	dense_matrix< bigint > X;
	X.reginvimage(A, B);
	return X;
}



//
// Image
//

//////////////////////////////
// BEGIN: INTERFACE - image //
//////////////////////////////

inline dense_matrix< bigint >
image(const dense_matrix< bigint > &B)
{
	dense_matrix< bigint > RET;
	RET.image1(B);
	return RET;
}



////////////////////////////
// END: INTERFACE - image //
////////////////////////////

inline dense_matrix< bigint >
image1(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > BILD;
	BILD.image1(A);
	return BILD;
}



inline dense_matrix< bigint >
image2(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > BILD;
	BILD.image(A);
	return BILD;
}



//
// InvImage
//

inline dense_matrix< bigint >
invimage(const dense_matrix< bigint > &A, const bigint *b)
{
	dense_matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



inline dense_matrix< bigint >
invimage(const dense_matrix< bigint > &A, const math_vector< bigint > &b)
{
	dense_matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



//
// solve
//

inline dense_matrix< bigint >
solve(const dense_matrix< bigint > &A, const bigint *b)
{
	dense_matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



inline dense_matrix< bigint >
solve(const dense_matrix< bigint > &A, const math_vector< bigint > &b)
{
	dense_matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



//
// Smith normal form
//

////////////////////////////
// BEGIN: INTERFACE - SNF //
////////////////////////////

inline dense_matrix< bigint >
snf(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snf_havas();
	return B;
}



inline dense_matrix< bigint >
snf(const dense_matrix< bigint > &A, dense_matrix< bigint > &B, dense_matrix< bigint > &C)
{
	dense_matrix< bigint > D = A;
	D.snf_havas(B, C);
	return D;
}



//////////////////////////
// END: INTERFACE - SNF //
//////////////////////////

inline dense_matrix< bigint >
snf_hartley(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snf_hartley();
	return B;
}



inline dense_matrix< bigint >
snf_hartley(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	    dense_matrix< bigint > &T2)
{
	dense_matrix< bigint > B = A;
	B.snf_hartley(T1, T2);
	return B;
}



inline dense_matrix< bigint >
snf_simple(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snf_simple();
	return B;
}



inline dense_matrix< bigint >
snf_simple(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	   dense_matrix< bigint > &T2)
{
	dense_matrix< bigint > B = A;
	B.snf_simple(T1, T2);
	return B;
}



inline dense_matrix< bigint >
snf_havas(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snf_havas();
	return B;
}



inline dense_matrix< bigint >
snf_havas(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	  dense_matrix< bigint > &T2)
{
	dense_matrix< bigint > B = A;
	B.snf_havas(T1, T2);
	return B;
}



inline dense_matrix< bigint >
snf_mult(const dense_matrix< bigint > &A, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_mult(art);
	return B;
}



inline dense_matrix< bigint >
snf_mult(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	 dense_matrix< bigint > &T2, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_mult(T1, T2, art);
	return B;
}



inline dense_matrix< bigint >
snf_add(const dense_matrix< bigint > &A, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_add(art);
	return B;
}



inline dense_matrix< bigint >
snf_add(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	dense_matrix< bigint > &T2, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_add(T1, T2, art);
	return B;
}



inline dense_matrix< bigint >
snf_new(const dense_matrix< bigint > &A, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_new(art);
	return B;
}



inline dense_matrix< bigint >
snf_new(const dense_matrix< bigint > &A, dense_matrix< bigint > &T1,
	dense_matrix< bigint > &T2, long art = 1)
{
	dense_matrix< bigint > B = A;
	B.snf_new(T1, T2, art);
	return B;
}



inline dense_matrix< bigint >
snfmod_dkt(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snfmod_dkt(B.latticedet2());
	return B;
}



inline dense_matrix< bigint >
snfmod_dkt(const dense_matrix< bigint > &A, const bigint &h)
{
	dense_matrix< bigint > B = A;
	B.snfmod_dkt(h);
	return B;
}



inline dense_matrix< bigint >
snfmod_cohen(const dense_matrix< bigint > &A)
{
	dense_matrix< bigint > B = A;
	B.snfmod_cohen(B.latticedet2());
	return B;
}



inline dense_matrix< bigint >
snfmod_cohen(const dense_matrix< bigint > &A, const bigint &h)
{
	dense_matrix< bigint > B = A;
	B.snfmod_cohen(h);
	return B;
}



/////////////////////////
// END: Linear algebra //
// PART 2              //
/////////////////////////

//
// mgcd Computation
//

inline bigint *mgcd1(const bigint *v, lidia_size_t l, dense_matrix< bigint > &T)
{
	return T.mgcd1(v, l);
}



inline bigint *mgcd2(const bigint *v, lidia_size_t l, dense_matrix< bigint > &T)
{
	return T.mgcd2(v, l);
}



inline bigint *mgcd3(const bigint *v, lidia_size_t l, dense_matrix< bigint > &T)
{
	return T.mgcd3(v, l);
}



inline bigint *mgcd4(const bigint *v, lidia_size_t l, dense_matrix< bigint > &T)
{
	return T.mgcd4(v, l);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DENSE_BIGINT_MATRIX_H_GUARD_
