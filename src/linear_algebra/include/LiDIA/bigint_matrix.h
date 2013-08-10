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


#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
#define LIDIA_BIGINT_MATRIX_H_GUARD_


#ifndef LIDIA_PRIMELIST_SIZE
#define LIDIA_PRIMELIST_SIZE 300
#endif

#ifndef LIDIA_PRIMELIST_MAXSIZE
#define LIDIA_PRIMELIST_MAXSIZE 86416
#endif



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_PATH_H_GUARD_
# include	"LiDIA/path.h"
#endif
#ifndef LIDIA_MATRIX_INTERN_H_GUARD_
# include	"LiDIA/matrix_intern.h"
#endif
#ifndef LIDIA_LONG_MATRIX_H_GUARD_
# include	"LiDIA/long_matrix.h"
#endif
#ifndef LIDIA_BIGFLOAT_MATRIX_H_GUARD_
# include	"LiDIA/bigfloat_matrix.h"
#endif
#ifndef LIDIA_RING_MATRIX_H_GUARD_
# include	"LiDIA/ring_matrix.h"
#endif
#ifndef LIDIA_DENSE_BIGINT_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_bigint_matrix_kernel.h"
#endif
#ifndef LIDIA_SPARSE_BIGINT_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_bigint_matrix_kernel.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define DRMKex DRMK< bigint >
#define SRMKex SRMK< bigint >

class file_adjoint;
class trans_matrix;



template <>
class matrix< bigint > : public ring_matrix< bigint >
{
	//
	// friend classes
	//

#if defined(_MSC_VER)
	friend template matrix< long >;
	friend template dense_bigint_matrix_kernel< matrix < bigint > >;
	friend template sparse_bigint_matrix_kernel< matrix < bigint > >;
#else
	friend class matrix< long >;
	friend class dense_bigint_matrix_kernel< matrix < bigint > >;
	friend class sparse_bigint_matrix_kernel< matrix < bigint > >;
#endif

private:

	//
	// dense representation
	//

	const dense_bigint_matrix_kernel< matrix < bigint > > D_bigint_modul;

	//
	// sparse representation
	//

	const sparse_bigint_matrix_kernel< matrix < bigint > > S_bigint_modul;

	dense_bigint_matrix_kernel< matrix < bigint > > hnf_ref_modul;

	//
	// constructors
	//

public:

	matrix()
		: ring_matrix< bigint > () {}
	matrix(const matrix_flags &flags)
		: ring_matrix< bigint > (flags) {}
	matrix(lidia_size_t a, lidia_size_t b)
		: ring_matrix< bigint > (a, b) {}
	matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: ring_matrix< bigint > (a, b, flags) {}
	matrix(const base_vector< bigint > &v)
		: ring_matrix< bigint > (v) {}
	matrix(const base_vector< bigint > &v, const matrix_flags &flags)
		: ring_matrix< bigint > (v, flags) {}
	matrix(const base_matrix< bigint > &A)
		: ring_matrix< bigint > (A) {}
	matrix(const base_matrix< bigint > &A, const matrix_flags &flags)
		: ring_matrix< bigint > (A, flags) {}
	matrix(const dense_base_matrix< bigint > &A)
		: ring_matrix< bigint > (A) {}
	matrix(const sparse_base_matrix< bigint > &A)
		: ring_matrix< bigint > (A) {}
	matrix(lidia_size_t a, lidia_size_t b, const bigint **v)
		: ring_matrix< bigint > (a, b, v) {}
	matrix(lidia_size_t a, lidia_size_t b, const bigint **v, const matrix_flags &flags)
		: ring_matrix< bigint > (a, b, v, flags) {}

	matrix(const base_matrix< long > &A);


	//
	// destructor
	//

	~matrix() {}

	//
	// Casts
	//

	operator bigfloat_matrix()
	{
		bigfloat_matrix B(rows, columns, bitfield);
		for (register lidia_size_t i = 0; i < rows; i++)
			for (register lidia_size_t j = 0; j < columns; j++)
				B.sto(i, j, bigfloat(value[i][j]));
		return B;
	}

	//
	// assignment
	//

	matrix< bigint > & operator = (const matrix< bigint > &B)
	{
		assign(B);
		return *this;
	}

	matrix< bigint > & operator = (const base_matrix< bigint > & B)
	{
		base_matrix< bigint >::assign(B);
		return *this;
	}

	//
	// pseudo-division
	//

	friend void divide(matrix< bigint > &, const matrix< bigint > &, const bigint &);
	friend void compwise_divide(matrix< bigint > &, const matrix< bigint > &,
				    const matrix< bigint > &);

protected:

	void divide(const matrix< bigint > &, const bigint &);
	void compwise_divide(const matrix< bigint > &, const matrix< bigint > &);

  //
  // pseudo - division
  //

public:

	friend matrix< bigint > operator / (const matrix< bigint > &, const bigint &);
	friend matrix< bigint > & operator /= (matrix< bigint > &, const bigint &);

	//
	// remainder
	//

	friend matrix< bigint > operator % (const matrix< bigint > &, const bigint &);

	friend matrix< bigint > & operator %= (matrix< bigint > &, const bigint &);

	friend void remainder(matrix< bigint > &, const matrix< bigint > &, const bigint &);
	friend void remainder(matrix< long > &, const matrix< bigint > &, long);

	//
	// norms and bounds
	//

	void max(bigint &) const;
	bigint max() const
	{
		bigint MAX;
		max(MAX);
		return MAX;
	}

	friend bigint max(const matrix< bigint > &);

	void max_abs(bigint &) const;
	bigint max_abs() const
	{
		bigint RES;
		max_abs(RES);
		return RES;
	}

	friend bigint max_abs(const matrix< bigint > &);

	void max_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint max_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MAX;
		max_pos(MAX, x, y);
		return MAX;
	}

	friend bigint max_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void max_abs_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint max_abs_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint RES;
		max_abs_pos(RES, x, y);
		return RES;
	}

	friend bigint max_abs_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void min(bigint &) const;
	bigint min() const
	{
		bigint MIN;
		min(MIN);
		return MIN;
	}

	friend bigint min(const matrix< bigint > &);

	void min_abs(bigint &) const;
	bigint min_abs() const
	{
		bigint MIN;
		min_abs(MIN);
		return MIN;
	}

	friend bigint min_abs(const matrix< bigint > &);

	void min_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint min_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MIN;
		min_pos(MIN, x, y);
		return MIN;
	}

	friend bigint min_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void min_abs_pos(bigint &, lidia_size_t &, lidia_size_t &) const;
	bigint min_abs_pos(lidia_size_t &x, lidia_size_t &y) const
	{
		bigint MIN;
		min_abs_pos(MIN, x, y);
		return MIN;
	}

	friend bigint min_abs_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

	void hadamard(bigint &) const;
	bigint hadamard() const
	{
		bigint H;
		hadamard(H);
		return H;
	}

	friend bigint hadamard(const matrix< bigint > &);

	void binary_hadamard(lidia_size_t &) const;
	lidia_size_t binary_hadamard() const
	{
		lidia_size_t H;
		binary_hadamard(H);
		return H;
	}

	friend lidia_size_t binary_hadamard(const matrix< bigint > &);

	void row_norm(bigint &, lidia_size_t, long) const;
	bigint row_norm(lidia_size_t i, long n) const
	{
		bigint RES;
		row_norm(RES, i, n);
		return RES;
	}

	friend bigint row_norm(const matrix< bigint > &, lidia_size_t, long);

	void column_norm(bigint &, lidia_size_t, long) const;
	bigint column_norm(lidia_size_t i, long n) const
	{
		bigint RES;
		column_norm(RES, i, n);
		return RES;
	}

	friend bigint column_norm(const matrix< bigint > &,
				  lidia_size_t, long);

	//
	// no_of_elements
	//

	bigint no_of_elements()
	{
		bigint RES, SUM = 0;
		for (lidia_size_t i = 0; i < columns; i++) {
			RES = 0;
			column_norm(RES, i, 0);
			SUM += RES;
		}
		return SUM;
	}

	//
	// randomize
	//

	void randomize(const bigint &);
	void randomize_with_det(const bigint &, const bigint &);
	void randomize(const bigint &, const long);

	///////////////////////////
	// BEGIN: Linear algebra //
	// PART 1                //
	///////////////////////////

	//
	// rank
	//

	lidia_size_t rank() const;
	lidia_size_t rank(const bigint &) const;

	friend lidia_size_t rank(const matrix< bigint > &);

	//
	// rank and linearly independent rows
	//

	///////////////////////////////
	// BEGIN: INTERFACE - lininr //
	///////////////////////////////


	lidia_size_t *lininr() const
	{
		return lininr1();
	}

	lidia_size_t *lininr(const bigint &H) const
	{
		return lininr1(H);
	}

	void lininr(base_vector< lidia_size_t > &v) const
	{
		lininr1(v);
	}

	friend lidia_size_t *lininr(const matrix< bigint > &);
	friend void lininr(base_vector< lidia_size_t > &v, const matrix< bigint > &);

	/////////////////////////////
	// END: INTERFACE - lininr //
	/////////////////////////////

private:

	lidia_size_t *lininr1() const;
	lidia_size_t *lininr1(const bigint &) const;
	void lininr1(base_vector< lidia_size_t > &) const;

	lidia_size_t *lininr2() const;
	lidia_size_t *lininr2(const bigint &) const;
	void lininr2(base_vector< lidia_size_t > &) const;

	//
	// rank and linearly independent columns
	//

	///////////////////////////////
	// BEGIN: INTERFACE - lininc //
	///////////////////////////////

public:

	lidia_size_t *lininc() const
	{
		return lininc1();
	}

	lidia_size_t *lininc(const bigint &H) const
	{
		return lininc1(H);
	}

	void lininc(base_vector< lidia_size_t > &v) const
	{
		lininc1(v);
	}

	friend lidia_size_t *lininc(const matrix< bigint > &);
	friend void lininc(base_vector< lidia_size_t > &, const matrix< bigint > &);

	////////////////////////////
	// END: INTERFACE - ininc //
	////////////////////////////

private:

	lidia_size_t *lininc1() const;
	lidia_size_t *lininc1(const bigint &) const;
	void lininc1(base_vector< lidia_size_t > &) const;

	lidia_size_t *lininc2() const;
	lidia_size_t *lininc2(const bigint &) const;
	void lininc2(base_vector< lidia_size_t > &) const;

	//
	// regular expansion
	//

public:

	void regexpansion(const lidia_size_t *);

	friend void regexpansion(matrix< bigint > &, const lidia_size_t *);

	//
	// adjoint matrix
	//

	////////////////////////////
	// BEGIN: INTERFACE - adj //
	////////////////////////////

	void adj(const matrix< bigint > &A)
	{
		adj1(A);
	}

	void adj(const matrix< bigint > &A, const bigint &H)
	{
		adj1(A, H);
	}

	void adj(const matrix< bigint > &A, const bigint &H, const bigint &DET)
	{
		adj1(A, H, DET);
	}

	friend matrix< bigint > adj(const matrix< bigint > &);

	//////////////////////////
	// END: INTERFACE - adj //
	//////////////////////////

private:

	void adj1(const matrix< bigint > &);
	void adj1(const matrix< bigint > &, const bigint &);
	void adj1(const matrix< bigint > &, const bigint &, const bigint &);

	void adj2(const matrix< bigint > &);
	void adj2(const matrix< bigint > &, const bigint &);
	void adj2(const matrix< bigint > &, const bigint &, const bigint &);

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

	void latticedet(bigint &DET, const bigint &H)
	{
		latticedet2(DET, H);
	}

	bigint latticedet()
	{
		bigint RET;
		latticedet2(RET);
		return RET;
	}

	friend bigint latticedet(const matrix< bigint > &);

	/////////////////////////////////
	// END: INTERFACE - latticedet //
	/////////////////////////////////


	void latticedet1(bigint &) const;
	void latticedet2(bigint &) const;
	void latticedet3(bigint &) const;
	void latticedet_special(bigint &) const;

	void latticedet1(bigint &, const bigint &) const;
	void latticedet2(bigint &, const bigint &) const;
	void latticedet3(bigint &, const bigint &) const;

	void real_latticedet(bigint &, const bigint &) const;

	//
	// determinant
	//

	void det(bigint &) const;
	void det(bigint &, const bigint &) const;
	void det(bigint &, const bigint &, int) const;


	bigint det() const
	{
		bigint DET;
		det(DET);
		return DET;
	}

	friend bigint det(const matrix< bigint > &);
	friend bigint det(const matrix< bigint > &, const bigint &);

	//
	// characteristic polynomial
	//

	bigint *charpoly() const;
	bigint *charpoly(const bigint &) const;
	void charpoly(base_vector< bigint > &) const;
	void charpoly(base_vector< bigint > &, const bigint &) const;

	friend bigint *charpoly(const matrix< bigint > &);
	friend void charpoly(base_vector< bigint > &, const matrix< bigint > &);

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


	void hnf()
	{
		hnf_havas();
	}

	void hnf(matrix< bigint > &A)
	{
		hnf_havas(A);
	}

	friend matrix< bigint > hnf(const matrix< bigint > &);
	friend matrix< bigint > hnf(const matrix< bigint > &, matrix< bigint > &);

	//////////////////////////
	// END: INTERFACE - HNF //
	//////////////////////////


	//
	// Group of Algorithms: Hermite
	//

	void hnf_cg(const matrix< long > &, long, const bigint &, int);
	void hnf_cg(const matrix< long > &, matrix< bigint > &, long, const bigint &, int);

	void hnf_ref();
	void hnf_gls_solver();
	void hnf_gls_solver(matrix< bigint > & TR);

	// HNF HAVAS
	void hnf_havas(lidia_size_t KernAlgo = 0, lidia_size_t mgcdModul = 5, lidia_size_t normalizeModul = 1);
	void hnf_havas(matrix< bigint > &, lidia_size_t KernAlgo = 0, lidia_size_t mgcdModul = 5, lidia_size_t normalizeModul = 1);

	void mgcd(lidia_size_t mgcdModul = 0);
	void mgcd(matrix< bigint > &, lidia_size_t mgcdModul = 0);

	void normalize(lidia_size_t normalizeModul = 1);

	// HNF HAVAS & MAJEWSKI
	void hnf_havas_cont()
	{
		hnf_havas();
	}

	void hnf_havas_cont(matrix< bigint > &A)
	{
		hnf_havas(A);
	}

	// HNF KANNAN
	void hnf_kannan(lidia_size_t SW = 0);
	void hnf_kannan(matrix< bigint > &, lidia_size_t SW = 0);

	// HNF STORJOHANN
	void hnf_storjohann();
	void hnf_storjohann(matrix< bigint > &TR, matrix< bigint > &C, matrix< bigint > &Q);

	// HNF Domich Kannan Trotter
	void hnfmod_dkt(const bigint &);
	void hnfmod_dkt(matrix< bigint > &, const bigint &);

	void hnfmod_dkt()
	{
		hnfmod_dkt(latticedet());
	}

	void hnfmod_dkt(matrix< bigint > &TR)
	{
		hnfmod_dkt(TR, latticedet());
	}

	// HNF Cohen
	void hnfmod_cohen(const bigint &);
	void hnfmod_cohen()
	{
		hnfmod_cohen(latticedet());
	}

	// HNF MUELLER (cf. Achim Mueller's Diplom thesis)
	void hnfmod_mueller(matrix< bigint > &);

	//
	// Kernel
	//

	///////////////////////////////
	// BEGIN: INTERFACE - kernel //
	///////////////////////////////

	void kernel(const matrix< bigint > &B)
	{
		kernel1(B);
	}

	friend matrix< bigint > kernel(const matrix< bigint > &);

	/////////////////////////////
	// END: INTERFACE - kernel //
	/////////////////////////////


protected:

	void kernel1(const matrix< bigint > &);
	void kernel2(const matrix< bigint > &);

	//
	// regular Invimage
	//

	////////////////////////////////////
	// BEGIN: INTERFACE - reginvimage //
	////////////////////////////////////

public:

	void reginvimage(const matrix< bigint > &A, const matrix< bigint > &B)
	{
		reginvimage1(A, B);
	}

	friend matrix< bigint > reginvimage(const matrix< bigint > &, const matrix< bigint > &);

	//////////////////////////////////
	// END: INTERFACE - reginvimage //
	//////////////////////////////////


	void reginvimage1(const matrix< bigint > &, const matrix< bigint > &);
	void reginvimage2(const matrix< bigint > &, const matrix< bigint > &);

	void reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b) const;
	void reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b,
			 const bigint &H, const bigint &DET) const;

	void reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b,
			 const bigint &H, const bigint &DET, const bigint &mod) const;

	void reginvimage_ZmZ(math_vector< bigint > &RES, const math_vector< bigint > & b,
			     const bigint &H, bigint &DET) const;

	//
	// Image
	//

	//////////////////////////////
	// BEGIN: INTERFACE - image //
	//////////////////////////////

	void image(const matrix< bigint > &B)
	{
		image1(B);
	}

	friend matrix< bigint > image(const matrix< bigint > &);

	////////////////////////////
	// END: INTERFACE - image //
	////////////////////////////


	void image1(const matrix< bigint > &);
	void image2(const matrix< bigint > &);

	//
	// InvImage
	//


	void invimage(const matrix< bigint > &, const bigint *);
	friend matrix< bigint > invimage(const matrix< bigint > &, const bigint *);

	void invimage(const matrix< bigint > &, const math_vector< bigint > &);
	friend matrix< bigint > invimage(const matrix< bigint > &, const math_vector< bigint > &);

	//
	// solve
	//

	void solve(const matrix< bigint > &A, const bigint *b)
	{
		invimage(A, b);
	}
	friend matrix< bigint > solve(const matrix< bigint > &, const bigint *);

	void solve(const matrix< bigint > &A, const math_vector< bigint > &b)
	{
		invimage(A, b);
	}

	friend matrix< bigint > solve(const matrix< bigint > &, const math_vector< bigint > &);

	//
	// Smith normal form
	//

	////////////////////////////
	// BEGIN: INTERFACE - SNF //
	////////////////////////////

	void snf()
	{
		snf_havas();
	}
	void snf(matrix< bigint > &A, matrix< bigint > &B)
	{
		snf_havas(A, B);
	}

	friend matrix< bigint > snf(const matrix< bigint > &);
	friend matrix< bigint > snf(const matrix< bigint > &, matrix< bigint > &, matrix< bigint > &);

	//////////////////////////
	// END: INTERFACE - SNF //
	//////////////////////////


	void snf_hartley();
	void snf_hartley(matrix< bigint > &, matrix< bigint > &);

	void snf_simple();
	void snf_simple(matrix< bigint > &, matrix< bigint > &);

	void snf_havas();
	void snf_havas(matrix< bigint > &, matrix< bigint > &);

	void snf_mult(long art = 1);
	void snf_mult(matrix< bigint > &, matrix< bigint > &, long art = 1);

	void snf_add(long art = 1);
	void snf_add(matrix< bigint > &, matrix< bigint > &, long art = 1);

	void snf_new(long art = 1);
	void snf_new(matrix< bigint > &, matrix< bigint > &, long art = 1);

	void snfmod_dkt(const bigint &);
	void snfmod_dkt()
	{
		snfmod_dkt(latticedet());
	}

	void snfmod_cohen(const bigint &);
	void snfmod_cohen()
	{
		snfmod_cohen(latticedet());
	}

	/////////////////////////
	// END: Linear algebra //
	// PART 2              //
	/////////////////////////

	//
	// basis completion
	//


	void basis_completion(bigint *, lidia_size_t);
	void simple_basis_completion(bigint *, lidia_size_t);

	//
	// conditioning functions
	//


	lidia_size_t cond_matrix(bigint *, lidia_size_t);

	//
	// gauss elimination
	//


	void gauss();

	//
	// mgcd Computation
	//


	bigint *mgcd(const bigint *, lidia_size_t);
	bigint *mgcd_new(const bigint *, lidia_size_t);
	bigint *mgcd_new2(const bigint *, lidia_size_t);
	bigint *mgcd_new3(const bigint *, lidia_size_t);
	bigint *mgcd_tree(const bigint *, lidia_size_t);

	friend bigint *mgcd1(const bigint *, lidia_size_t, matrix< bigint > &);
	bigint *mgcd1(const bigint *, lidia_size_t);

	friend bigint *mgcd2(const bigint *, lidia_size_t);
	friend void mgcd2(bigint &, const bigint *, lidia_size_t);

	friend bigint *mgcd2(const bigint *, lidia_size_t, matrix< bigint > &);
	bigint *mgcd2(const bigint *, lidia_size_t);

	friend bigint *mgcd3(const bigint *, lidia_size_t, matrix< bigint > &);
	bigint *mgcd3(const bigint *, lidia_size_t);

	friend bigint *mgcd4(const bigint *, lidia_size_t, matrix< bigint > &);
	bigint *mgcd4(const bigint *, lidia_size_t);


	// ADDED BY MJJ - quadratic order hnf stuff

	void adj(file_adjoint &, const bigint &, const bigint &);
	void adj(file_adjoint &, const bigint &, const bigint &, int);

	void size_red_jacobs();
	void size_red_jacobs(trans_matrix &);

	friend void pre_reduction(matrix< long > &, lidia_size_t, lidia_size_t,
				  base_vector< lidia_size_t > &);
	friend void pre_reduction(matrix< long > &, trans_matrix &,
				  lidia_size_t, lidia_size_t,
				  base_vector< lidia_size_t > &);

	void post_reduction(base_vector< lidia_size_t > &,
			    lidia_size_t *RET);

	lidia_size_t *hnf_jacobs0(trans_matrix &, bigint &, int);
	lidia_size_t *hnf_jacobs1(trans_matrix &, bigint &, int, bigint &, timer &);

	lidia_size_t *hnf_cg2(const matrix< long > &, int);
	lidia_size_t *hnf_cg2(const matrix< long > &, int, bool &);
	lidia_size_t *hnf_cg2(const matrix< long > &, trans_matrix &, int, bigint &, int);
	lidia_size_t *hnf_cg2(const matrix< long > &, trans_matrix &, int, bigint &, int, bool &);

	lidia_size_t *hnf_cg3(const matrix< bigint > &);
	lidia_size_t *hnf_cg3_mod(const matrix< bigint > &, bool &);

	lidia_size_t *hnf_cg3(const matrix< bigint > &, trans_matrix &);
	lidia_size_t *hnf_cg3_mod(const matrix< bigint > &, trans_matrix &, int,
				  bool &);

	bool solve_hnf(math_vector< bigint > & b, math_vector< bigint > & x);

	//
	// transpose function for matrxi< bigint >
	//

	matrix< bigint > trans() const;
	void trans(const matrix< bigint > & B);

};


// friend functions

//
// pseudo-division
//

void divide(matrix< bigint > &, const matrix< bigint > &, const bigint &);
void compwise_divide(matrix< bigint > &, const matrix< bigint > &,
			    const matrix< bigint > &);
matrix< bigint > operator / (const matrix< bigint > &, const bigint &);
matrix< bigint > & operator /= (matrix< bigint > &, const bigint &);

//
// remainder
//

matrix< bigint > operator % (const matrix< bigint > &, const bigint &);

matrix< bigint > & operator %= (matrix< bigint > &, const bigint &);

void remainder(matrix< bigint > &, const matrix< bigint > &, const bigint &);
void remainder(matrix< long > &, const matrix< bigint > &, long);

//
// norms and bounds
//

bigint max(const matrix< bigint > &);
bigint max_abs(const matrix< bigint > &);
bigint max_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);
bigint max_abs_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

bigint min(const matrix< bigint > &);
bigint min_abs(const matrix< bigint > &);
bigint min_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);
bigint min_abs_pos(const matrix< bigint > &, lidia_size_t &, lidia_size_t &);

bigint hadamard(const matrix< bigint > &);
lidia_size_t binary_hadamard(const matrix< bigint > &);
bigint row_norm(const matrix< bigint > &, lidia_size_t, long);
bigint column_norm(const matrix< bigint > &,
			  lidia_size_t, long);

///////////////////////////
// BEGIN: Linear algebra //
// PART 1                //
///////////////////////////

//
// rank
//

lidia_size_t rank(const matrix< bigint > &);

//
// rank and linearly independent rows
//

lidia_size_t *lininr(const matrix< bigint > &);
void lininr(base_vector< lidia_size_t > &v, const matrix< bigint > &);

lidia_size_t *lininc(const matrix< bigint > &);
void lininc(base_vector< lidia_size_t > &, const matrix< bigint > &);

//
// regular expansion
//

void regexpansion(matrix< bigint > &, const lidia_size_t *);

//
// adjoint matrix
//

matrix< bigint > adj(const matrix< bigint > &);

//
// lattice determinant
//

bigint latticedet(const matrix< bigint > &);

//
// determinant
//

bigint det(const matrix< bigint > &);
bigint det(const matrix< bigint > &, const bigint &);

//
// characteristic polynomial
//

bigint *charpoly(const matrix< bigint > &);
void charpoly(base_vector< bigint > &, const matrix< bigint > &);

//
// Hermite normal form
//

matrix< bigint > hnf(const matrix< bigint > &);
matrix< bigint > hnf(const matrix< bigint > &, matrix< bigint > &);

//
// Kernel
//

matrix< bigint > kernel(const matrix< bigint > &);

//
// regular Invimage
//

matrix< bigint > reginvimage(const matrix< bigint > &, const matrix< bigint > &);

//
// Image
//

matrix< bigint > image(const matrix< bigint > &);

//
// InvImage
//

matrix< bigint > invimage(const matrix< bigint > &, const bigint *);
matrix< bigint > invimage(const matrix< bigint > &, const math_vector< bigint > &);

//
// solve
//

matrix< bigint > solve(const matrix< bigint > &, const bigint *);
matrix< bigint > solve(const matrix< bigint > &, const math_vector< bigint > &);

//
// Smith normal form
//

matrix< bigint > snf(const matrix< bigint > &);
matrix< bigint > snf(const matrix< bigint > &, matrix< bigint > &, matrix< bigint > &);

//
// mgcd Computation
//

bigint *mgcd1(const bigint *, lidia_size_t, matrix< bigint > &);

bigint *mgcd2(const bigint *, lidia_size_t);
void mgcd2(bigint &, const bigint *, lidia_size_t);
bigint *mgcd2(const bigint *, lidia_size_t, matrix< bigint > &);

bigint *mgcd3(const bigint *, lidia_size_t, matrix< bigint > &);

bigint *mgcd4(const bigint *, lidia_size_t, matrix< bigint > &);

// ADDED BY MJJ - quadratic order hnf stuff
void pre_reduction(matrix< long > &, lidia_size_t, lidia_size_t,
			  base_vector< lidia_size_t > &);
void pre_reduction(matrix< long > &, trans_matrix &,
			  lidia_size_t, lidia_size_t,
			  base_vector< lidia_size_t > &);


//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// division
//

inline void
divide(matrix< bigint > &RES, const matrix< bigint > &M, const bigint &a)
{
	RES.divide(M, a);
}



inline void
compwise_divide(matrix< bigint > &RES, const matrix< bigint > &M,
		const matrix< bigint > &N)
{
	RES.compwise_divide(M, N);
}



/////////////////////////////////
// END: arithmetic procedures  //
// BEGIN: arithmetic operators //
/////////////////////////////////

//
// division
//

inline matrix< bigint >
operator / (const matrix< bigint > &A, const bigint &m)
{
	matrix< bigint > RES(A.rows, A.columns);
	RES.divide(A, m);
	return RES;
}



inline matrix< bigint > &
operator /= (matrix< bigint > &A, const bigint &m)
{
	A.divide(A, m);
	return A;
}



//
// remainder
//

inline matrix< bigint >
operator % (const matrix< bigint > &A, const bigint &mod)
{
	matrix< bigint > RES(A.rows, A.columns);
	remainder(RES, A, mod);
	return RES;
}



inline matrix< bigint > &
operator %= (matrix< bigint > &A, const bigint &mod)
{
	remainder(A, A, mod);
	return A;
}



///////////////////////////////
// END: arithmetic operators //
///////////////////////////////

//
// norms and bounds
//

inline bigint
max(const matrix< bigint > &A)
{
	return A.max();
}



inline bigint
max_abs(const matrix< bigint > &A)
{
	return A.max_abs();
}



inline bigint
max_pos(const matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.max_pos(x, y);
}



inline bigint
max_abs_pos(const matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.max_abs_pos(x, y);
}



inline bigint
min(const matrix< bigint > &A)
{
	return A.min();
}



inline bigint
min_abs(const matrix< bigint > &A)
{
	return A.min_abs();
}



inline bigint
min_pos(const matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.min_pos(x, y);
}



inline bigint
min_abs_pos(const matrix< bigint > &A, lidia_size_t &x, lidia_size_t &y)
{
	return A.min_abs_pos(x, y);
}



inline bigint
hadamard(const matrix< bigint > &A)
{
	bigint H;
	A.hadamard(H);
	return H;
}



inline lidia_size_t
binary_hadamard(const matrix< bigint > &A)
{
	lidia_size_t H;
	A.binary_hadamard(H);
	return H;
}



inline bigint
row_norm(const matrix< bigint > &A, lidia_size_t i, long n)
{
	bigint RES;
	A.row_norm(RES, i, n);
	return RES;
}



inline bigint
column_norm(const matrix< bigint > &A, lidia_size_t i, long n)
{
	bigint RES;
	A.column_norm(RES, i, n);
	return RES;
}



///////////////////////////
// BEGIN: Linear algebra //
// PART 1                //
///////////////////////////

//
// rank
//

inline lidia_size_t rank(const matrix< bigint > &B)
{
	return B.rank();
}



//
// rank and linearly independent rows
//

inline lidia_size_t *
lininr(const matrix< bigint > &B)
{
	return B.lininr();
}



inline void
lininr(base_vector< lidia_size_t > &v, const matrix< bigint > &B)
{
	B.lininr(v);
}



//
// rank and linearly independent columns
//

inline lidia_size_t *
lininc(const matrix< bigint > &B)
{
	return B.lininc();
}



inline void
lininc(base_vector< lidia_size_t > &v, const matrix< bigint > &B)
{
	B.lininc(v);
}



//
// regular expansion
//

inline void
regexpansion(matrix< bigint > &A, const lidia_size_t *v)
{
	A.regexpansion(v);
}



//
// adjoint matrix
//

inline matrix< bigint >
adj(const matrix< bigint > &A)
{
	matrix< bigint > B(A.rows, A.columns);
	B.adj(A);
	return B;
}



//
// lattice determinant
//

inline bigint latticedet(const matrix< bigint > &A)
{
	bigint DET;
	A.latticedet2(DET);
	return DET;
}



//
// determinant
//

inline bigint
det(const matrix< bigint > &A)
{
	bigint DET;
	A.det(DET);
	return DET;
}



//
// characteristic polynomial
//

inline bigint *
charpoly(const matrix< bigint > &A)
{
	return A.charpoly();
}



inline void
charpoly(base_vector< bigint > &v, const matrix< bigint > &A)
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

inline matrix< bigint >
hnf(const matrix< bigint > &A)
{
	matrix< bigint > B = A;
	B.hnf_havas();
	return B;
}



inline matrix< bigint >
hnf(const matrix< bigint > &A, matrix< bigint > &TR)
{
	matrix< bigint > B = A;
	B.hnf_havas(TR);
	return B;
}



//
// Kernel
//

inline matrix< bigint >
kernel(const matrix< bigint > &B)
{
	matrix< bigint > RET;
	RET.kernel1(B);
	return RET;
}



//
// regular Invimage
//

inline matrix< bigint >
reginvimage(const matrix< bigint > &A, const matrix< bigint > &B)
{
	matrix< bigint > X;
	X.reginvimage(A, B);
	return X;
}



//
// Image
//

inline matrix< bigint >
image(const matrix< bigint > &B)
{
	matrix< bigint > RET;
	RET.image1(B);
	return RET;
}



//
// InvImage
//

inline matrix< bigint >
invimage(const matrix< bigint > &A, const bigint *b)
{
	matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



inline matrix< bigint >
invimage(const matrix< bigint > &A, const math_vector< bigint > &b)
{
	matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



//
// solve
//

inline matrix< bigint >
solve(const matrix< bigint > &A, const bigint *b)
{
	matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



inline matrix< bigint >
solve(const matrix< bigint > &A, const math_vector< bigint > &b)
{
	matrix< bigint > B;
	B.invimage(A, b);
	return B;
}



//
// Smith normal form
//

inline matrix< bigint >
snf(const matrix< bigint > &A)
{
	matrix< bigint > B = A;
	B.snf_havas();
	return B;
}



inline matrix< bigint >
snf(const matrix< bigint > &A, matrix< bigint > &B, matrix< bigint > &C)
{
	matrix< bigint > D = A;
	D.snf_havas(B, C);
	return D;
}



/////////////////////////
// END: Linear algebra //
// PART 2              //
/////////////////////////

//
// mgcd Computation
//

inline bigint *mgcd1(const bigint *v, lidia_size_t l, matrix< bigint > &T)
{
	return T.mgcd1(v, l);
}



inline bigint *mgcd2(const bigint *v, lidia_size_t l, matrix< bigint > &T)
{
	return T.mgcd2(v, l);
}



inline bigint *mgcd3(const bigint *v, lidia_size_t l, matrix< bigint > &T)
{
	return T.mgcd3(v, l);
}



inline bigint *mgcd4(const bigint *v, lidia_size_t l, matrix< bigint > &T)
{
	return T.mgcd4(v, l);
}



//
// transpose function for ring_matrix< bigint >
//

inline void
matrix< bigint >::trans(const matrix< bigint > & B)
{
	base_matrix< bigint >::trans(B);
}



inline matrix< bigint >
matrix< bigint >::trans() const
{
	return base_matrix< bigint >::trans();
}



inline matrix< bigint >
trans(const matrix< bigint > &A)
{
	return A.trans();
}



#undef DRMKex
#undef SRMKex

typedef matrix< bigint > bigint_matrix;



// ADDED BY MJJ for transformation matrix computation

//
// Class: file_adjoint
//
// This class represents an adjoint matrix using several files to store rows
//      of the matrix.  The rows are stored in vector form modulo small primes
//      p.  The original rows are recovered using Chinese remaindering.
//

class file_adjoint
{
protected:

	int index; // which instance is this?
	lidia_size_t rows; // number of rows in files
	lidia_size_t columns; // number of columns in files;
	lidia_size_t actual_rows; // number of rows (actual)
	lidia_size_t actual_columns; // number of columns (actual)
	lidia_size_t num_files; // number of files to use
	lidia_size_t rows_per_file; // number of rows stored per file
	base_vector< bigint > mods; // vector of modulii used



private:

	void combine_rows(matrix< bigint > & RES, lidia_size_t idx);



public:

	file_adjoint();
	~file_adjoint();

	void kill();

	void init(lidia_size_t nrows, lidia_size_t ncols, lidia_size_t nfiles);

	lidia_size_t get_no_of_files();
	lidia_size_t get_no_of_rows();
	lidia_size_t get_no_of_columns();

	void touch_files();

	void resize(lidia_size_t nrows, lidia_size_t ncols);

	void add_new_prime(matrix< bigint > & A, bigint & p);

	bool test_adjoint(matrix< bigint > & A, const bigint & DET);
	bool test_adjoint_diag(matrix< bigint > & A, const bigint & DET);

	void multiply(math_vector< bigint > & RES, math_vector< bigint > & x);
	void multiply(matrix< bigint > & RES, matrix< bigint > & X);
};



//
// Class: trans_matrix
//
// This class represents a power product of unimodular matrices corresponding
//      to a Hermite normal form transformation matrix.  The individual
//      matrices are either stored in files or in a vector.  This class allows
//      individual columns to be accessed transparently, and matrix-vector
//      products to be computed.
//

class trans_matrix
{
protected:

	int size; // number of matrices
	bool files; // are we storing in files?
	base_vector< matrix < bigint > > TR; // vector of matrices
	int index; // which instance is this?
	lidia_size_t cols; // # of columns
	lidia_size_t last_cols; // # of columns in last matrix

	bigint diag_element; // special value for diagonal
	lidia_size_t notone_idx; // index of matrix with
	//   special value on diagonal
	lidia_size_t notone_rows; // # or rows in the matrix
	//   we need to consider
	bool adj_files; // is the adjoint in files?
	lidia_size_t adj_index; // which matrix is the adjoint?

	bigint DET_div; // determinante to divide by
	lidia_size_t det_idx; // corresponding mat idx


public:

	file_adjoint ADJ; // adjoint representation


	trans_matrix();
	~trans_matrix();

	void kill();

	void set_mode(int mode);

	void set_adjoint_mode(int mode);
	bool get_adjoint_mode();

	void touch_files();

	void remove_columns(lidia_size_t *rcols);

	void set_no_of_columns(lidia_size_t ncols);
	lidia_size_t get_no_of_columns();

	lidia_size_t get_size();
	void set_size(lidia_size_t nsize);

	void get_matrix(matrix< bigint > & tran, lidia_size_t idx);

	void store_matrix(matrix< bigint > & tran);
	void store_matrix(matrix< bigint > & tran, bigint & DET);
	void store_adjoint();
	void store_det(bigint & DET);

	void get_column_vector(math_vector< bigint > & vec, lidia_size_t idx);
	void get_submatrix(matrix< bigint > & V, lidia_size_t *idx);

	math_vector< bigint > multiply(math_vector< bigint > & x);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGINT_MATRIX_H_GUARD_
