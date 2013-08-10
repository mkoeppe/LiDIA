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
//	Author	: Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGFLOAT_MATRIX_H_GUARD_
#define LIDIA_BIGFLOAT_MATRIX_H_GUARD_


#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BIGFLOAT_INT_H_GUARD_
# include	"LiDIA/bigfloat_int.h"
#endif
#ifndef LIDIA_MATRIX_INTERN_H_GUARD_
# include	"LiDIA/matrix_intern.h"
#endif
#ifndef LIDIA_FIELD_MATRIX_H_GUARD_
# include	"LiDIA/field_matrix.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// MM: included bigfloat_int, because it is used as a template
//     argument in bigfloat_matrix.c and the HP CC compiler
//     only examines the header file for class definitions.

template<>
class matrix< bigfloat > : public field_matrix< bigfloat >
{

	//
	// constructors
	//

public:

	matrix()
		: field_matrix< bigfloat > () {}
	matrix(const matrix_flags &flags)
		: field_matrix< bigfloat > (flags) {}
	matrix(lidia_size_t a, lidia_size_t b)
		: field_matrix< bigfloat > (a, b) {}
	matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: field_matrix< bigfloat > (a, b, flags) {}
	matrix(const base_matrix< bigfloat > & A)
		: field_matrix< bigfloat > (A) {}
	matrix(const base_matrix< bigfloat > &A, const matrix_flags &flags)
		: field_matrix< bigfloat > (A, flags) {}
	matrix(const dense_base_matrix< bigfloat > &A)
		: field_matrix< bigfloat > (A) {}
	matrix(const sparse_base_matrix< bigfloat > &A)
		: field_matrix< bigfloat > (A) {}
	matrix(const base_vector< bigfloat > &v)
		: field_matrix< bigfloat > (v) {}
	matrix(const base_vector< bigfloat > &v, const matrix_flags &flags)
		: field_matrix< bigfloat > (v, flags) {}
	matrix(lidia_size_t a, lidia_size_t b, const bigfloat **v)
		: field_matrix< bigfloat > (a, b, v) {}
	matrix(lidia_size_t a, lidia_size_t b, const bigfloat **v, const matrix_flags &flags)
		: field_matrix< bigfloat > (a, b, v, flags) {}

	//
	// destructor
	//

	~matrix() {}

	//
	// public functions
	//

	void randomize();

	int equal(const matrix< bigfloat > &) const;
	friend bool operator == (const matrix< bigfloat > &, const matrix< bigfloat > &);
	friend bool operator != (const matrix< bigfloat > &, const matrix< bigfloat > &);

	void invert(matrix< bigfloat > &) const;

	void gauss_jordan(base_vector< bigfloat > &, const base_vector< bigfloat > &) const;
	base_vector< bigfloat > gauss_jordan(const base_vector< bigfloat > &) const;

	void cholesky(matrix< bigfloat > &) const;
	void mod_cholesky(matrix< bigfloat > &) const;

	void LR(matrix< bigfloat > &, matrix< bigfloat > &) const;
	void QR(matrix< bigfloat > &, matrix< bigfloat > &) const;

	friend void gauss_jordan(base_vector< bigfloat > & x, const matrix< bigfloat > &, const base_vector< bigfloat > &);
	friend base_vector< bigfloat > gauss_jordan(const matrix< bigfloat > &, const base_vector< bigfloat > &);
	friend void invert(base_matrix< bigfloat > &, const base_matrix< bigfloat > &);
};



inline bool
operator == (const matrix< bigfloat > & A, const matrix< bigfloat > & B)
{
	return A.equal(B);
}



inline bool
operator != (const matrix< bigfloat > & A, const matrix< bigfloat > & B)
{
	return !A.equal(B);
}



inline void
gauss_jordan (base_vector< bigfloat > & x, const matrix< bigfloat > & A, const base_vector< bigfloat > & b)
{
	A.gauss_jordan(x, b);
}



inline base_vector< bigfloat >
gauss_jordan (const matrix< bigfloat > & A, const base_vector< bigfloat > & b)
{
	return A.gauss_jordan(b);
}



inline void
invert (base_matrix< bigfloat > & B, const base_matrix< bigfloat > & A)
{
	matrix< bigfloat > C(A), D;
	C.invert(D);
	B.assign(D);
}



typedef matrix< bigfloat > bigfloat_matrix;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGFLOAT_MATRIX_H_GUARD_
