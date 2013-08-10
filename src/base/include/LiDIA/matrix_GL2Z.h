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
//	Author	: Thomas Papanikolaou (TP), Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MATRIX_GL2Z_H_GUARD_
#define LIDIA_MATRIX_GL2Z_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BASE_MATRIX_H_GUARD_
# include	"LiDIA/base_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class matrix_GL2Z
{
	friend class quadratic_form;
private:

	// (s u)
	// (t v)
	//
	// determinant = sv - tu = +-1

	bigint s, u;
	bigint t, v;
	mutable int determinant;


public:

	//
	// c'tors and d'tor
	//

	matrix_GL2Z();
	matrix_GL2Z(const bigint & a, const bigint & c,
		    const bigint & b, const bigint & d);
	matrix_GL2Z(const matrix_GL2Z & a);
	matrix_GL2Z(const base_matrix< bigint > & U);
	~matrix_GL2Z();


	//
	// accessors
	//

	const bigint &  get_s() const;
	const bigint &  get_t() const;
	const bigint &  get_u() const;
	const bigint &  get_v() const;

	int det() const;

	const bigint & at(int, int) const;
	const bigint & operator () (int, int) const;



	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(const matrix_GL2Z & A);
	void assign(const base_matrix< bigint > & A);

	matrix_GL2Z & operator = (const matrix_GL2Z & a);
	matrix_GL2Z & operator = (const base_matrix< bigint > & U);



	//
	// comparators
	//

	bool is_zero() const;
	bool is_one() const;

	bool is_equal(const matrix_GL2Z & A) const;



	//
	// modifiers
	//

	void invert();

	void es(const bigint &);
	void te();
	void randomize(const bigint &);



	friend void multiply(matrix_GL2Z & c,
			     const matrix_GL2Z & a,
			     const matrix_GL2Z & b);

	friend matrix_GL2Z operator *(const matrix_GL2Z & a, const matrix_GL2Z & b);
	friend matrix_GL2Z operator /(const matrix_GL2Z & a, const matrix_GL2Z & b);

	friend matrix_GL2Z inverse(const matrix_GL2Z & a);

	friend std::ostream & operator << (std::ostream & out, const matrix_GL2Z & a);
	friend std::istream & operator >> (std::istream & in, matrix_GL2Z & A);

};



//
// c'tors and d'tor
//

inline
matrix_GL2Z::matrix_GL2Z ()
	: s(1L),
	  u(0L),
	  t(0L),
	  v(1L),
	  determinant(1)
{
	// nothing to do
}



inline
matrix_GL2Z::matrix_GL2Z (const matrix_GL2Z & A)
	: s(A.s),
	  u(A.u),
	  t(A.t),
	  v(A.v),
	  determinant(A.determinant)
{
	// nothing to do
}



inline
matrix_GL2Z::~matrix_GL2Z()
{
	// nothing to do
}



//
// accessors
//

inline const bigint &
matrix_GL2Z::get_s () const
{
	return s;
}



inline const bigint &
matrix_GL2Z::get_t () const
{
	return t;
}



inline const bigint &
matrix_GL2Z::get_u () const
{
	return u;
}



inline const bigint &
matrix_GL2Z::get_v () const
{
	return v;
}



inline const bigint &
matrix_GL2Z::operator () (int i, int j) const
{
	return at(i, j);
}



//
// assigners
//

inline void
matrix_GL2Z::assign_zero ()
{
	s.assign_zero();
	u.assign_zero();
	t.assign_zero();
	v.assign_zero();
	determinant = 0;
}



inline void
matrix_GL2Z::assign_one ()
{
	s.assign_one();
	u.assign_zero();
	t.assign_zero();
	v.assign_one();
	determinant = 1;
}



inline void
matrix_GL2Z::assign (const matrix_GL2Z & A)
{
	s.assign(A.s);
	t.assign(A.t);
	u.assign(A.u);
	v.assign(A.v);
	determinant = A.determinant;
}



inline matrix_GL2Z &
matrix_GL2Z::operator = (const matrix_GL2Z & A)
{
	assign(A);
	return *this;
}



inline matrix_GL2Z &
matrix_GL2Z::operator = (const base_matrix< bigint > & A)
{
	assign(A);
	return *this;
}



//
// comparators
//

inline bool
matrix_GL2Z::is_zero () const
{
	return (s.is_zero() &&
		u.is_zero() &&
		t.is_zero() &&
		v.is_zero());
}



inline bool
matrix_GL2Z::is_one () const
{
	return (s.is_one() &&
		u.is_zero() &&
		t.is_zero() &&
		v.is_one());
}



inline bool
matrix_GL2Z::is_equal (const matrix_GL2Z & A) const
{
	return (s.compare(A.s) == 0 &&
		u.compare(A.u) == 0 &&
		t.compare(A.t) == 0 &&
		v.compare(A.v) == 0);
}



inline bool
operator == (const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	return A.is_equal(B);
}



inline bool
operator != (const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	return !A.is_equal(B);
}



//
// modifiers
//

inline void
matrix_GL2Z::es (const bigint & a)
{
	u += a * s;
	v += a * t;
}



inline void
matrix_GL2Z::te ()
{
	swap(s, u);
	u.negate();
	swap(t, v);
	v.negate();
}



//
// procedural arithmetic
//

void multiply(matrix_GL2Z & c, const matrix_GL2Z & a, const matrix_GL2Z & b);



inline void
divide (matrix_GL2Z & C, const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	matrix_GL2Z inv_B(B);

	inv_B.invert();
	multiply(C, A, inv_B);
}



//
// operational arithmetic
//

inline matrix_GL2Z
operator * (const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	matrix_GL2Z C;

	multiply(C, A, B);
	return C;
}



inline matrix_GL2Z
operator / (const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	matrix_GL2Z C;

	divide(C, A, B);
	return C;
}



inline matrix_GL2Z &
operator *= (matrix_GL2Z & A, const matrix_GL2Z & B)
{
	multiply(A, A, B);
	return A;
}



inline matrix_GL2Z &
operator /= (matrix_GL2Z & A, const matrix_GL2Z & B)
{
	divide(A, A, B);
	return A;
}



inline matrix_GL2Z
inverse (const matrix_GL2Z & A)
{
	matrix_GL2Z inv_A(A);

	inv_A.invert();
	return inv_A;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MATRIX_GL2Z_H_GUARD_
