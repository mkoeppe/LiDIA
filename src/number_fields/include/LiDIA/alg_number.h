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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ALG_NUMBER_H_GUARD_
#define LIDIA_ALG_NUMBER_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif

#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif

#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif
#ifndef LIDIA_BIGMOD_MATRIX_H_GUARD_
# include	"LiDIA/bigmod_matrix.h"
#endif
#ifndef LIDIA_BIGFLOAT_MATRIX_H_GUARD_
# include	"LiDIA/bigfloat_matrix.h"
#endif

#ifndef LIDIA_BIGINT_POLYNOMIAL_H_GUARD_
# include	"LiDIA/bigint_polynomial.h"
#endif
#ifndef LIDIA_BIGCOMPLEX_POLYNOMIAL_H_GUARD_
# include	"LiDIA/bigcomplex_polynomial.h"
#endif



#ifdef LIDIA_IMPLICIT_CAST_EXPLICIT
#define alg_ideal_cast(O) alg_ideal(O)
#else
#define alg_ideal_cast(O) O
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class number_field;
class order;
class alg_number;
class module;
class alg_ideal;
class prime_ideal;



class nf_base {
	// The classes allowed to change bases get immediate access.
	// They need to do all the work, since otherwise reference counting
	// gets difficult. So everything is private, since nobody is allowed to
	// explicitly use this class !!

	friend class number_field; // Those need access to all
	friend class order; // information of class nf_base.

	friend class alg_number; // Those need access only to
	friend class module; // the variable "current_base".
	friend class alg_ideal;
	friend class prime_ideal;

	// Additional friend functions:
	friend bool operator == (const number_field &, const number_field &);
	friend std::ostream& operator << (std::ostream&, const number_field&);
	friend std::istream& operator >> (std::istream&, number_field&);
	friend const bigint & disc(const order & O); // discriminant
	friend std::ostream& operator << (std::ostream &, const order &);
	friend std::istream& operator >> (std::istream &, order &);
	friend void multiply(alg_number &, const alg_number &, const alg_number &);
	friend void square(alg_number &, const alg_number &);
	friend void multiply2(module &, const module &, const module &);
	friend bigint_matrix rep_matrix(const alg_number &);
	friend std::istream& operator >> (std::istream &, alg_number &);
	friend std::istream& operator >> (std::istream &, module &);

	static nf_base * current_base; // Default nf_base
	static nf_base * dummy_base; // Dummy nf_base
	int references; // How many objects use this base.

	mutable polynomial< bigint > f; // the irreducible polynomial
	mutable lidia_size_t real_roots; // number of real roots

	mutable bigint_matrix base;
	mutable bigint den;
	mutable math_matrix< bigint > table; // multiplication table

	// Internal information, that is frequently needed and so stored
	// with the order, as soon, as it is computed:
	mutable math_vector< bigint > One; // Representation of 1
	bigfloat_matrix  conjugates; // column i contains the conjugates
	// of w_i (0 <= i < degree().

	// Functions for checking, what is already computed:
	bool table_computed() const
	{
		return (table.get_no_of_columns() > 1);
	}

	void compute_table() const;

	bool base_computed() const
	{
		return (base.get_no_of_columns() > 1);
	}

	void compute_base() const;

	void compute_conjugates();

	// whether to use MT or polynomial for multiplying.
	bool using_necessary() const
	{
		return (f.is_zero() || base_computed());
	}

	const bigfloat & get_conjugate(lidia_size_t, lidia_size_t);
	const bigfloat_matrix & get_conjugates();

	const polynomial< bigint > & which_polynomial() const;

	lidia_size_t degree() const;

	// Constructor and Destructor

	nf_base():references(0), f(), real_roots(0), base(), den(1), table(),
		  One(), conjugates()
	{}			// Do nothing.
	~nf_base() {}		// Just call member destructors.

	//  Default copy constructor and assignment are OK.

	void assign(const nf_base&);
	// Input and output are used by number_field and order!
	friend std::ostream& operator << (std::ostream&, const nf_base&); // jeweils Polynom
	friend std::istream& operator >> (std::istream&, nf_base&); // ein-, ausgeben

public:
	const math_vector< bigint > & get_one() const; // 1
	void dec_ref()
	{
		references--;
		if (!references && this != nf_base::dummy_base) {
			if (nf_base::current_base == this)
				nf_base::current_base = nf_base::dummy_base;
			delete this;
		}
	}
	void inc_ref()
	{
		references++;
	}
};

// friend functions of nf_base

bool operator == (const number_field &, const number_field &);
std::ostream& operator << (std::ostream&, const number_field&);
std::istream& operator >> (std::istream&, number_field&);
const bigint & disc(const order & O); // discriminant
std::ostream& operator << (std::ostream &, const order &);
std::istream& operator >> (std::istream &, order &);
void multiply(alg_number &, const alg_number &, const alg_number &);
void square(alg_number &, const alg_number &);
void multiply2(module &, const module &, const module &);
bigint_matrix rep_matrix(const alg_number &);
std::istream& operator >> (std::istream &, alg_number &);
std::istream& operator >> (std::istream &, module &);

    // Input and output are used by number_field and order!
std::ostream& operator << (std::ostream&, const nf_base&); // jeweils Polynom
std::istream& operator >> (std::istream&, nf_base&); // ein-, ausgeben



class number_field {
#ifdef LIDIA_DEBUG
	static int count;
#endif
	nf_base * base; // pointer to the base of O

public:
	// Contructors & destructor:
	number_field();
	number_field(const polynomial< bigint > & p);
	number_field(const bigint *v, lidia_size_t deg);
	number_field(const order &);
	number_field(const number_field &);
	~number_field();

	// Cast operator
	operator nf_base *() const
	{
		return base;
	}

	// member functions:
	const bigfloat & get_conjugate(lidia_size_t, lidia_size_t);
	const bigfloat_matrix & get_conjugates();

	const polynomial< bigint > & which_polynomial() const
	{
		return base->which_polynomial();
	}

	nf_base * which_base() const
	{
		return base;
	}

	lidia_size_t degree() const
	{
		return base->degree();
	}

	lidia_size_t no_of_real_embeddings() const;

#ifndef HEADBANGER
	friend number_field overfield(const number_field &, const number_field &);
	void assign(const number_field & K);
#endif

	number_field & operator = (const number_field & F)
	{
		assign(F);
		return *this;
	}

	// Comparisions:
	// Since deciding whether fields are isomorphic is difficult,
	// we consider fields to be equal only, if they are generated by
	// the same polynomial.
	friend bool operator == (const number_field &, const number_field &);
	friend bool operator != (const number_field &, const number_field &);
	friend bool operator <= (const number_field &, const number_field &);
	friend bool operator < (const number_field &, const number_field &);
	friend bool operator >= (const number_field &, const number_field &);
	friend bool operator > (const number_field &, const number_field &);

	// friend functions:
	friend const polynomial< bigint > & which_polynomial(const number_field &);
	friend nf_base * which_base(const number_field &);
	friend lidia_size_t degree(const number_field &);
	friend lidia_size_t no_of_real_embeddings(const number_field &);
	// In-/Output:
	friend std::ostream& operator << (std::ostream&, const number_field&); // jeweils Polynom
	friend std::istream& operator >> (std::istream&, number_field&); // ein-, ausgeben
};

// friend functions of number_field

#ifndef HEADBANGER
number_field overfield(const number_field &, const number_field &);
#endif

    // Comparisions:
    // Since deciding whether fields are isomorphic is difficult,
    // we consider fields to be equal only, if they are generated by
    // the same polynomial.
bool operator == (const number_field &, const number_field &);
bool operator != (const number_field &, const number_field &);
bool operator <= (const number_field &, const number_field &);
bool operator < (const number_field &, const number_field &);
bool operator >= (const number_field &, const number_field &);
bool operator > (const number_field &, const number_field &);

const polynomial< bigint > & which_polynomial(const number_field &);
nf_base * which_base(const number_field &);
lidia_size_t degree(const number_field &);
lidia_size_t no_of_real_embeddings(const number_field &);
    // In-/Output:
std::ostream& operator << (std::ostream&, const number_field&); // jeweils Polynom
std::istream& operator >> (std::istream&, number_field&); // ein-, ausgeben


inline bool
operator != (const number_field & A, const number_field & B)
{
	return !(A == B);
}



inline bool
operator < (const number_field & A, const number_field & B)
{
	return (!(A == B) && (A <= B));
}



inline bool
operator >= (const number_field & A, const number_field & B)
{
	return B <= A;
}



inline bool
operator > (const number_field & A, const number_field & B)
{
	return B < A;
}



inline const polynomial< bigint > &
which_polynomial (const number_field &F)
{
	return F.which_polynomial();
}



inline nf_base *
which_base (const number_field & F)
{
	return F.base;
}



inline lidia_size_t
degree (const number_field &F)
{
	return F.degree();
}



inline lidia_size_t
no_of_real_embeddings (const number_field &F)
{
	return F.no_of_real_embeddings();
}



class order {
#ifdef LIDIA_DEBUG
	static int count;
#endif
#ifndef HEADBANGER
	// Components of the data - type:
	mutable bigint discriminant;
	nf_base * base; // pointer to the base of O

	// Internal routine for the work involved in comparisons
	void compare(const order&, bool &, bool &) const;
#endif

public:
	// Contructors & destructor:
	order(const nf_base * base1 = nf_base::current_base);
	order(const polynomial< bigint > & p, // initialize with equation order
	      const matrix< bigint > & = bigint_matrix(), // (default)
	      const bigint & = 1); // or with transformed eq. order.
	order(const matrix< bigint > &); // table uebergeben
	order(const matrix< bigint > &, const bigint &,
	      const nf_base * base1 = nf_base::current_base); // TRAFO uebergeben.
	order(const order &); // copy constructor
	order(const number_field &); // Maximal order !
	~order();

	order & operator = (const order & O)
	{
		assign(O);
		return *this;
	}

	// member functions:
	// Accessing Information:

	const bigint_matrix & base_numerator() const
	{
		return base->base;
	}

	const bigint & base_denominator() const
	{
		return base->den;
	}

	const polynomial< bigint > & which_polynomial() const
	{
		return base->which_polynomial();
	}

	nf_base * which_base() const
	{
		return base;
	}

	const bigint & MT(lidia_size_t, lidia_size_t, lidia_size_t); // accessing MT

	const bigfloat & get_conjugate(lidia_size_t, lidia_size_t);
	const bigfloat_matrix & get_conjugates();

	lidia_size_t degree() const			// degree of extension
	{
		return base->degree();
	}

	lidia_size_t no_of_real_embeddings() const;

	lidia_size_t no_of_roots_of_unity();

#ifndef HEADBANGER
	void assign(const order &);
#endif

	// Cast operator:
	operator alg_ideal() const;

	operator nf_base *() const
	{
		return base;
	}

	// Comparisions:
	// We are interested in comparing orders only if they are
	// over the same field! So nothing else is supported!!
	friend bool operator == (const order&, const order&);
	friend bool operator != (const order&, const order&);
	friend bool operator <= (const order&, const order&);
	friend bool operator < (const order&, const order&);
	friend bool operator >= (const order&, const order&);
	friend bool operator > (const order&, const order&);

	// friend functions:
	friend const bigint_matrix & base_numerator(const order &);
	friend const bigint & base_denominator(const order &);
	friend const polynomial< bigint > & which_polynomial(const order &);
	friend nf_base * which_base(const order &);
	friend void swap(order &, order &);

	// Number-theoretic functions:
	friend lidia_size_t degree(const order & O); // degree of field extension
	friend lidia_size_t no_of_real_embeddings(const order & O);
	friend lidia_size_t no_of_roots_of_unity(order & O);
	friend const bigint & disc(const order & O); // discriminant
	module pseudo_radical(const bigint & p) const;
	// computes pseudo-radical
	bool dedekind(const bigint & p, polynomial< bigint > & h2) const;
	// uses Dedekind criterion for prime p; returns true,
	// if p is an index divisor and sets extended to the new order.

	order maximize(const bigint & p) const;
	// maximizes order at p.
	order maximize() const; // full maximization (round2)
	order maximize2() const; // full maximization (round2) with a
				// sligtly different strategy.

	// In-/Output:
	friend std::ostream& operator << (std::ostream &, const order &); // table ausgeben
	friend std::istream& operator >> (std::istream &, order &);
	// tries to interpret input either as multiplication table
	// or as base transformation according to dimension.
};

// friend functions of class order

// Comparisions:
// We are interested in comparing orders only if they are
// over the same field! So nothing else is supported!!
bool operator == (const order&, const order&);
bool operator != (const order&, const order&);
bool operator <= (const order&, const order&);
bool operator < (const order&, const order&);
bool operator >= (const order&, const order&);
bool operator > (const order&, const order&);

const bigint_matrix & base_numerator(const order &);
const bigint & base_denominator(const order &);
const polynomial< bigint > & which_polynomial(const order &);
nf_base * which_base(const order &);
void swap(order &, order &);

// Number-theoretic functions:
lidia_size_t degree(const order & O); // degree of field extension
lidia_size_t no_of_real_embeddings(const order & O);
lidia_size_t no_of_roots_of_unity(order & O);
const bigint & disc(const order & O); // discriminant

// In-/Output:
std::ostream& operator << (std::ostream &, const order &); // table ausgeben
std::istream& operator >> (std::istream &, order &);
// tries to interpret input either as multiplication table
// or as base transformation according to dimension.



inline const bigint_matrix &
base_numerator (const order & O)
{
	return O.base_numerator();
}



inline const bigint &
base_denominator (const order & O)
{
	return O.base_denominator();
}



inline const polynomial< bigint > &
which_polynomialc(const order & O)
{
	return O.which_polynomial();
}



inline nf_base *
which_base (const order & O)
{
	return O.base;
}



inline lidia_size_t
degree (const order & O)
{
	return O.degree();
}



inline lidia_size_t
no_of_real_embeddings (const order & O)
{
	return O.no_of_real_embeddings();
}



inline lidia_size_t
no_of_roots_of_unity (order & O)
{
	return O.no_of_roots_of_unity();
}



class alg_number {
#ifdef LIDIA_DEBUG
	static int count;
#endif
	bigint den;
	math_vector< bigint > coeff;
	nf_base * O;

public:
	// Constructors & destructor:
	alg_number(const nf_base * O1 = nf_base::current_base);
	alg_number(const bigint &, const nf_base * O1 = nf_base::current_base);
	alg_number(const base_vector< bigint > &, const bigint & i,
		   const nf_base * O1 = nf_base::current_base);
	alg_number(const bigint *, const bigint & i = 1,
		   const nf_base * O1 = nf_base::current_base);
	alg_number(const alg_number & a);
	~alg_number();

	alg_number & operator = (const alg_number & a)
	{
		assign(a);
		return *this;
	}

	// member-functions
	const bigint & denominator() const
	{
		return den;
	}

	const math_vector< bigint > & coeff_vector() const
	{
		return coeff;
	}

	alg_number numerator() const
	{
		return alg_number(coeff, bigint(1), O);
	}

	nf_base * which_base() const
	{
		return O;
	}

	lidia_size_t degree() const
	{
		return O->degree();
	}

	bigfloat get_conjugate(lidia_size_t) const;
	math_vector< bigfloat > get_conjugates() const;

#ifndef HEADBANGER
	bool is_zero() const;
	bool is_one() const;

	void normalize();
	void negate();
	void invert();
	void multiply_by_2();
	void divide_by_2();

	void assign_zero(); // Remains member of same order
	void assign_one(); // Remains member of same order
	void assign(const bigint &); // Remains member of same order
	void assign(const alg_number &); // Becomes member of same order as a!!

	// Procedural versions of arithmetic operations:
	friend void add(alg_number &, const alg_number &, const alg_number &);
	friend void subtract(alg_number &, const alg_number &, const alg_number &);
	friend void multiply(alg_number &, const alg_number &, const alg_number &);
	friend void divide(alg_number &, const alg_number &, const alg_number &);

	friend void add(alg_number &, const alg_number &, const bigint &);
	friend void subtract(alg_number &, const alg_number &, const bigint &);
	friend void multiply(alg_number &, const alg_number &, const bigint &);
	friend void divide(alg_number &, const alg_number &, const bigint &);

	friend void add(alg_number &, const bigint &, const alg_number &);
	friend void subtract(alg_number &, const bigint &, const alg_number &);
	friend void multiply(alg_number &, const bigint &, const alg_number &);
	friend void divide(alg_number &, const bigint &, const alg_number &);

	friend void power(alg_number &, const alg_number &, const bigint &);
#endif
	friend void power_mod_p(alg_number &, const alg_number &,
				const bigint &, const bigint &);

	// arithmetic operators:
	friend alg_number operator -(const alg_number &);
	friend alg_number operator +(const alg_number &, const alg_number &);
	friend alg_number operator +(const alg_number &, const bigint &);
	friend alg_number operator +(const bigint &, const alg_number &);
	friend alg_number operator -(const alg_number &, const alg_number &);
	friend alg_number operator -(const alg_number &, const bigint &);
	friend alg_number operator -(const bigint &, const alg_number &);
	friend alg_number operator *(const alg_number &, const alg_number &);
	friend alg_number operator *(const alg_number &, const bigint &);
	friend alg_number operator *(const bigint &, const alg_number &);
	friend alg_number operator /(const alg_number &, const alg_number &);
	friend alg_number operator /(const alg_number &, const bigint &);
	friend alg_number operator /(const bigint &, const alg_number &);

	alg_number & operator += (const alg_number & a)
	{
		add(*this, *this, a);
		return *this;
	}

	alg_number & operator += (const bigint & a)
	{
		add(*this, *this, a);
		return *this;
	}

	alg_number & operator -= (const alg_number & a)
	{
		subtract(*this, *this, a);
		return *this;
	}

	alg_number & operator -= (const bigint & a)
	{
		subtract(*this, *this, a);
		return *this;
	}

	alg_number & operator *= (const alg_number & a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	alg_number & operator *= (const bigint & a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	alg_number & operator /= (const alg_number & a)
	{
		divide(*this, *this, a);
		return *this;
	}

	alg_number & operator /= (const bigint & a)
	{
		divide(*this, *this, a);
		return *this;
	}

	// Comparisions:
	// By now, only comparision of numbers over the same order is implemented.
	friend bool operator == (const alg_number&, const alg_number&);
	friend bool operator != (const alg_number&, const alg_number&);

	bool operator ! () const
	{
		return is_zero();
	}

	// Some number-theoretic functions:
	friend lidia_size_t degree(const alg_number &);
	friend bigint_matrix rep_matrix(const alg_number &);
	friend bigrational norm(const alg_number &); // Norm
	friend bigrational trace(const alg_number &); // Trace
	friend polynomial< bigint > charpoly(const alg_number &);
	// characteristic polynomial

	// other functions:
#ifndef HEADBANGER
	friend void negate(alg_number &, const alg_number &);
	friend void invert(alg_number &, const alg_number &);
	friend alg_number inverse(const alg_number &);
#endif
	friend const bigint & denominator(const alg_number &);
	friend const math_vector< bigint > & coeff_vector(const alg_number &);
	friend alg_number numerator(const alg_number &);
	friend nf_base * which_base(const alg_number &);
#ifndef HEADBANGER
	friend void square(alg_number &, const alg_number &);
#endif
	friend void swap(alg_number &, alg_number &);

	// random numbers
	void randomize(const bigint &);

	// In-/Output:
	friend std::ostream& operator << (std::ostream&, const alg_number&);
	friend std::istream& operator >> (std::istream&, alg_number&);

	// friends:
	friend class module;
	//  friend module operator *(const module &, const module &);
	friend void multiply(alg_ideal &, const alg_ideal &, const alg_number &);
	friend void multiply(alg_ideal &, const alg_number &, const alg_ideal &);
	friend void divide(alg_ideal &, const alg_ideal &, const alg_number &);

};

// friends of class alg_number

#ifndef HEADBANGER
// Procedural versions of arithmetic operations:
void add(alg_number &, const alg_number &, const alg_number &);
void subtract(alg_number &, const alg_number &, const alg_number &);
void multiply(alg_number &, const alg_number &, const alg_number &);
void divide(alg_number &, const alg_number &, const alg_number &);

void add(alg_number &, const alg_number &, const bigint &);
void subtract(alg_number &, const alg_number &, const bigint &);
void multiply(alg_number &, const alg_number &, const bigint &);
void divide(alg_number &, const alg_number &, const bigint &);

void add(alg_number &, const bigint &, const alg_number &);
void subtract(alg_number &, const bigint &, const alg_number &);
void multiply(alg_number &, const bigint &, const alg_number &);
void divide(alg_number &, const bigint &, const alg_number &);

void power(alg_number &, const alg_number &, const bigint &);
#endif
void power_mod_p(alg_number &, const alg_number &,
		 const bigint &, const bigint &);

// arithmetic operators:
alg_number operator -(const alg_number &);
alg_number operator +(const alg_number &, const alg_number &);
alg_number operator +(const alg_number &, const bigint &);
alg_number operator +(const bigint &, const alg_number &);
alg_number operator -(const alg_number &, const alg_number &);
alg_number operator -(const alg_number &, const bigint &);
alg_number operator -(const bigint &, const alg_number &);
alg_number operator *(const alg_number &, const alg_number &);
alg_number operator *(const alg_number &, const bigint &);
alg_number operator *(const bigint &, const alg_number &);
alg_number operator /(const alg_number &, const alg_number &);
alg_number operator /(const alg_number &, const bigint &);
alg_number operator /(const bigint &, const alg_number &);

// Comparisions:
// By now, only comparision of numbers over the same order is implemented.
bool operator == (const alg_number&, const alg_number&);
bool operator != (const alg_number&, const alg_number&);

// Some number-theoretic functions:
lidia_size_t degree(const alg_number &);
bigint_matrix rep_matrix(const alg_number &);
bigrational norm(const alg_number &); // Norm
bigrational trace(const alg_number &); // Trace
polynomial< bigint > charpoly(const alg_number &);
// characteristic polynomial

// other functions:
#ifndef HEADBANGER
void negate(alg_number &, const alg_number &);
void invert(alg_number &, const alg_number &);
alg_number inverse(const alg_number &);
#endif
const bigint & denominator(const alg_number &);
const math_vector< bigint > & coeff_vector(const alg_number &);
alg_number numerator(const alg_number &);
nf_base * which_base(const alg_number &);
#ifndef HEADBANGER
void square(alg_number &, const alg_number &);
#endif
void swap(alg_number &, alg_number &);

// In-/Output:
std::ostream& operator << (std::ostream&, const alg_number&);
std::istream& operator >> (std::istream&, alg_number&);

void multiply(alg_ideal &, const alg_ideal &, const alg_number &);
void multiply(alg_ideal &, const alg_number &, const alg_ideal &);
void divide(alg_ideal &, const alg_ideal &, const alg_number &);



inline alg_number
operator - (const alg_number & a)
{
	alg_number c(-a.coeff, a.den, a.O);


	return c;
}



inline alg_number
operator +(const alg_number & a, const alg_number & b)
{
	alg_number c(a.O);


	add(c, a, b);
	return c;
}



inline alg_number
operator + (const alg_number & a, const bigint & b)
{
	alg_number c(a.O);


	add(c, a, b);
	return c;
}



inline alg_number
operator + (const bigint & b, const alg_number & a)
{
	alg_number c(a.O);


	add(c, a, b);
	return c;
}



inline alg_number
operator - (const alg_number & a, const alg_number & b)
{
	alg_number c(a.O);


	subtract(c, a, b);
	return c;
}



inline alg_number
operator - (const alg_number & a, const bigint & b)
{
	alg_number c(a.O);


	subtract(c, a, b);
	return c;
}



inline alg_number
operator - (const bigint & b, const alg_number & a)
{
	alg_number c(a.O);


	subtract(c, b, a);
	return c;
}



inline alg_number
operator * (const alg_number & a, const alg_number & b)
{
	alg_number c(a.O);
	multiply(c, a, b);
	return c;
}



inline alg_number
operator * (const alg_number & a, const bigint & b)
{
	alg_number c(a.O);


	multiply(c, a, b);
	return c;
}



inline alg_number
operator * (const bigint & b, const alg_number & a)
{
	alg_number c(a.O);


	multiply(c, a, b);
	return c;
}



inline alg_number
operator / (const alg_number & a, const alg_number & b)
{
	alg_number c(a.O);


	divide(c, a, b);
	return c;
}



inline alg_number
operator / (const alg_number & a, const bigint & b)
{
	alg_number c(a.O);


	divide(c, a, b);
	return c;
}



inline alg_number
operator / (const bigint & a, const alg_number & b)
{
	alg_number c(b.O);


	divide(c, a, b);
	return c;
}



inline bool
operator != (const alg_number & a, const alg_number & b)
{
	return !(a == b);
}



inline lidia_size_t
degree(const alg_number & a)
{
	return a.degree();
}



inline const bigint &
denominator(const alg_number & a)
{
	return a.den;
}



inline const math_vector< bigint > &
coeff_vector (const alg_number & a)
{
	return a.coeff;
}



inline alg_number
numerator (const alg_number & a)
{
	return alg_number(a.coeff, bigint(1), a.O);
}



inline nf_base *
which_base (const alg_number & a)
{
	return a.O;
}



class module {
protected:
#ifdef LIDIA_DEBUG
	static int count;
#endif
	mutable bigmod_matrix base;
	bigint den;
	nf_base * O;
	mutable bool is_exp;

	void compare(const module&, bool &, bool &) const;
	// internal routine for comparisons

public:
	// Contructors & destructor:
	module(const nf_base * O1 = nf_base::current_base); // zero module
	module(const alg_number & a,
	       const alg_number & b = alg_number(bigint(0)));
	module(const matrix< bigint > &, const bigint & d = 1,
	       const nf_base * O1 = nf_base::current_base);
	module(const bigmod_matrix &, const bigint & d = 1,
	       const nf_base * O1 = nf_base::current_base);
	module(const module &);
	virtual ~module();

	virtual module & operator = (const module & A)
	{
		assign(A);
		return *this;
	}

	// member-functions
	const bigint & denominator() const
	{
		return den;
	}

	const bigmod_matrix & coeff_matrix() const
	{
		return base;
	}

	bigint_matrix z_basis() const
	{
		bigmod_matrix tmp(base);

		tmp.lift(0);
		return tmp;
	}

	module numerator() const
	{
		return module(base, bigint(1), O);
	}

	nf_base * which_base() const
	{
		return O;
	}

	lidia_size_t degree() const
	{
		return O->degree();
	}

#ifndef HEADBANGER
	bool is_zero() const; // In the obvious sense
#endif

	bool is_whole_order() const;

#ifndef HEADBANGER
	bool is_one() const		// i.e. is_whole_order
	{
		return is_whole_order();
	}
#endif

	void normalize();

#ifndef HEADBANGER
	void invert();

	void assign_zero(); // Remains in the same order

	inline void assign_whole_order(); // Remains in the same order

	void assign_one()		// Remains in the same order, is the same as:
	{
		assign_whole_order();
	}

	void assign(const bigint &); // Remains in the same order
	void assign(const alg_number &); // Becomes subset of same order as a !
	virtual void assign(const module &); // Becomes subset of same order as a !

	// Procedural versions of arithmetic operations:
	friend void add(module &, const module &, const module &);
	friend void intersect(module &, const module &, const module &);
	friend void multiply(module &, const module &, const module &);
	//  friend void multiply2(module &, const module &, const module &);
	friend void multiply(module &, const module &, const bigint &);
	friend void multiply(module &, const bigint &, const module &);
	friend void divide(module &, const module &, const module &); // Warning:
	friend void divide(alg_ideal &, const alg_ideal &, const module &); // Warning:
	friend void divide(module &, const module &, const bigint &);
	friend void remainder(module &, const module &, const bigint &);

	friend void power(module &, const module &, const bigint &);
#endif

	// arithmetic operators:
	friend module operator +(const module &, const module &); // sum or union
	friend module operator &(const module &, const module &); // intersection
	friend module operator *(const module &, const module &); // product
	friend module operator *(const module &, const bigint &); // product
	friend module operator *(const bigint &, const module &); // product
	friend module operator /(const module &, const bigint &); // quotient
	friend module operator /(const module &, const module &); // quotient (??):
	// Division by an order is only guaranteed to produce correct results, if
	// you are in the maximal order!!
	friend module operator %(const module &, const bigint &); // Reduce mod pO

	module& operator += (const module & a)
	{
		add(*this, *this, a);
		return *this;
	}

	module& operator &= (const module & a)
	{
		intersect(*this, *this, a);
		return *this;
	}

	module& operator *= (const module & a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	module& operator *= (const bigint & a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	module& operator /= (const bigint & a)
	{
		divide(*this, *this, a);
		return *this;
	}

	module& operator /= (const module & a)	// Warning:
	{
		divide(*this, *this, a);
		return *this;
	}

	module& operator %= (const bigint & a)
	{
		remainder(*this, *this, a);
		return *this;
	}

	// Comparisions:
	// By now, only comparision of modules over the same order is implemented.
	friend bool operator == (const module&, const module&);
	friend bool operator != (const module&, const module&);
	friend bool operator <= (const module&, const module&);
	friend bool operator < (const module&, const module&);
	friend bool operator >= (const module&, const module&);
	friend bool operator > (const module&, const module&);

	bool operator ! () const
	{
		return is_zero();
	}

	// Some number-theoretic function:
	friend lidia_size_t degree(const module &);
	friend bigrational norm(const module &); // Norm
	friend bigrational exponent(const module &); // exponent
	order ring_of_multipliers(const bigint &p) const;

	// other functions:
#ifndef HEADBANGER
	friend void invert(module &, const module &);
	friend module inverse (const module &);
#endif

	friend const bigint & denominator(const module &);
	friend const bigmod_matrix & coeff_matrix(const module &);
	friend bigint_matrix z_basis(const module &);
	friend module numerator(const module &);
	friend nf_base * which_base(const module &);
#ifndef HEADBANGER
	friend void square(module &, const module &);
#endif
	friend void swap(module &, module &);

	// random numbers
	void randomize(const bigint &);

	// In-/Output:
	friend std::ostream& operator << (std::ostream &, const module &);
	friend std::istream& operator >> (std::istream &, module &);

	// friends:
	friend void multiply(alg_ideal &, const alg_ideal &, const alg_ideal &);
};

// friend functions of class module

#ifndef HEADBANGER
// Procedural versions of arithmetic operations:
void add(module &, const module &, const module &);
void intersect(module &, const module &, const module &);
void multiply(module &, const module &, const module &);
//void multiply2(module &, const module &, const module &);
void multiply(module &, const module &, const bigint &);
void multiply(module &, const bigint &, const module &);
void divide(module &, const module &, const module &); // Warning:
void divide(alg_ideal &, const alg_ideal &, const module &); // Warning:
void divide(module &, const module &, const bigint &);
void remainder(module &, const module &, const bigint &);

void power(module &, const module &, const bigint &);
#endif

// arithmetic operators:
module operator +(const module &, const module &); // sum or union
module operator &(const module &, const module &); // intersection
module operator *(const module &, const module &); // product
module operator *(const module &, const bigint &); // product
module operator *(const bigint &, const module &); // product
module operator /(const module &, const bigint &); // quotient
module operator /(const module &, const module &); // quotient (??):
// Division by an order is only guaranteed to produce correct results, if
// you are in the maximal order!!
module operator %(const module &, const bigint &); // Reduce mod pO

// Comparisions:
// By now, only comparision of modules over the same order is implemented.
bool operator == (const module&, const module&);
bool operator != (const module&, const module&);
bool operator <= (const module&, const module&);
bool operator < (const module&, const module&);
bool operator >= (const module&, const module&);
bool operator > (const module&, const module&);

// Some number-theoretic function:
lidia_size_t degree(const module &);
bigrational norm(const module &); // Norm
bigrational exponent(const module &); // exponent

// other functions:
#ifndef HEADBANGER
void invert(module &, const module &);
module inverse (const module &);
#endif

const bigint & denominator(const module &);
const bigmod_matrix & coeff_matrix(const module &);
bigint_matrix z_basis(const module &);
module numerator(const module &);
nf_base * which_base(const module &);
#ifndef HEADBANGER
void square(module &, const module &);
#endif
void swap(module &, module &);

// In-/Output:
std::ostream& operator << (std::ostream &, const module &);
std::istream& operator >> (std::istream &, module &);

void multiply(alg_ideal &, const alg_ideal &, const alg_ideal &);


inline module
operator + (const module & a, const module & b)
{
	module c(a.O);


	add(c, a, b);
	return c;
}



inline module
operator & (const module & a, const module & b)
{
	module c(a.O);


	intersect(c, a, b);
	return c;
}



inline module
operator * (const module & a, const module & b)
{
	module c(a.O);


	multiply(c, a, b);
	return c;
}



inline module
operator * (const module &a, const bigint &b)
{
	module c(a.O);


	multiply(c, a, b);
	return c;
}



inline module
operator * (const bigint &b, const module &a)
{
	module c(a.O);


	multiply(c, a, b);
	return c;
}



inline module
operator / (const module & a, const module & b)
{
	module c(a.O);


	divide(c, a, b);
	return c;
}



inline module
operator / (const module & a, const bigint & b)
{
	module c(a.O);


	divide(c, a, b);
	return c;
}



inline module
operator % (const module &a, const bigint &p)
{
	module c(a.O);


	remainder(c, a, p);
	return c;
}



inline lidia_size_t
degree (const module & a)
{
	return a.degree();
}



inline const bigint &
denominator (const module & a)
{
	return a.den;
}



inline const bigmod_matrix &
coeff_matrix (const module & a)
{
	return a.base;
}



inline bigint_matrix
z_basis (const module & a)
{
	return a.z_basis();
}



inline module
numerator (const module & a)
{
	return module(a.base, bigint(1), a.O);
}



inline nf_base *
which_base (const module & a)
{
	return a.O;
}



class alg_ideal : public module {
public:
	// Contructors & destructor:
	alg_ideal(const nf_base * O1 = nf_base::current_base); // zero ideal
	alg_ideal(const bigint & a,
		  const alg_number & b = alg_number(bigint(0)));
	alg_ideal(const alg_number & a,
		  const alg_number & b = alg_number(bigint(0)));
	alg_ideal(const matrix< bigint > &, const bigint & d = 1,
		  const nf_base * O1 = nf_base::current_base);
	alg_ideal(const bigmod_matrix &, const bigint & d = 1,
		  const nf_base * O1 = nf_base::current_base);
	alg_ideal(const alg_ideal &);
	virtual ~alg_ideal();

	module & operator = (const module & A)
	{
		multiply(*this, A, alg_ideal_cast(order(A.which_base())));
		return *this;
	}

	void assign(const module & A)
	{
		multiply(*this, A, alg_ideal_cast(order(A.which_base())));
	}

	alg_ideal & operator = (const alg_ideal & A)
	{
		assign(A);
		return *this;
	}

	// member-functions

	alg_ideal numerator() const
	{
		return alg_ideal(base, bigint(1), O);
	}

#ifndef HEADBANGER
	// Procedural versions of arithmetic operations:
	friend void multiply(alg_ideal &, const alg_ideal &, const alg_ideal &);
	friend void square(alg_ideal &, const alg_ideal &);
	friend void multiply(alg_ideal &, const alg_ideal &, const alg_number &);
	friend void multiply(alg_ideal &, const alg_number &, const alg_ideal &);
	friend void divide(alg_ideal &, const alg_ideal &, const module &); // Warning:
	friend void divide(alg_ideal &, const alg_ideal &, const alg_number &);
	friend long ord(const prime_ideal &, const alg_ideal &);
#endif

	// arithmetic operators:
	friend alg_ideal operator +(const alg_ideal &, const alg_ideal &); // sum or union
	friend alg_ideal operator &(const alg_ideal &, const alg_ideal &); // intersection
	friend alg_ideal operator *(const alg_ideal &, const alg_ideal &); // product
	friend alg_ideal operator *(const alg_ideal &, const bigint &); // product
	friend alg_ideal operator *(const bigint &, const alg_ideal &); // product
	friend alg_ideal operator *(const alg_ideal &, const alg_number &); // product
	friend alg_ideal operator *(const alg_number &, const alg_ideal &); // product
	friend alg_ideal operator /(const alg_ideal &, const bigint &); // quotient
	friend alg_ideal operator /(const alg_ideal &, const alg_number &); // quotient
	friend alg_ideal operator /(const alg_ideal &, const module &); // quotient (??):
	// Division by an order is only guaranteed to produce correct results, if
	// you are in the maximal order!!

	alg_ideal& operator += (const alg_ideal & a)
	{
		add(*this, *this, a);
		return *this;
	}

	alg_ideal& operator &= (const alg_ideal & a)
	{
		intersect(*this, *this, a);
		return *this;
	}

	alg_ideal& operator *= (const alg_ideal & a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	alg_ideal& operator *= (const bigint & a)
	{
		multiply(*this, *(static_cast<module *>(this)), a);
		return *this;
	}

	alg_ideal& operator /= (const bigint & a)
	{
		divide(*this, *(static_cast<module *>(this)), a);
		return *this;
	}

	alg_ideal& operator /= (const module & a)	// Warning:
	{
		divide(*this, *this, a);
		return *this;
	}

	alg_ideal& operator %= (const bigint & a)
	{
		remainder(*this, *this, a);
		return *this;
	}

	// reduction:
	void reduce(alg_number & divisor);
	alg_number reduce()
	{
		alg_number div;

		reduce(div);
		return div;
	}

	// other functions:
#ifndef HEADBANGER
	friend alg_ideal inverse (const alg_ideal &);
#endif

	friend alg_ideal numerator(const alg_ideal &);

	// random numbers
	void randomize(const bigint &);
};

// friend functions of class alg_ideal

#ifndef HEADBANGER
// Procedural versions of arithmetic operations:
void multiply(alg_ideal &, const alg_ideal &, const alg_ideal &);
void square(alg_ideal &, const alg_ideal &);
void multiply(alg_ideal &, const alg_ideal &, const alg_number &);
void multiply(alg_ideal &, const alg_number &, const alg_ideal &);
void divide(alg_ideal &, const alg_ideal &, const module &); // Warning:
void divide(alg_ideal &, const alg_ideal &, const alg_number &);
long ord(const prime_ideal &, const alg_ideal &);
#endif

// arithmetic operators:
alg_ideal operator +(const alg_ideal &, const alg_ideal &); // sum or union
alg_ideal operator &(const alg_ideal &, const alg_ideal &); // intersection
alg_ideal operator *(const alg_ideal &, const alg_ideal &); // product
alg_ideal operator *(const alg_ideal &, const bigint &); // product
alg_ideal operator *(const bigint &, const alg_ideal &); // product
alg_ideal operator *(const alg_ideal &, const alg_number &); // product
alg_ideal operator *(const alg_number &, const alg_ideal &); // product
alg_ideal operator /(const alg_ideal &, const bigint &); // quotient
alg_ideal operator /(const alg_ideal &, const alg_number &); // quotient
alg_ideal operator /(const alg_ideal &, const module &); // quotient (??):
// Division by an order is only guaranteed to produce correct results, if
// you are in the maximal order!!

// other functions:
#ifndef HEADBANGER
alg_ideal inverse (const alg_ideal &);
#endif
alg_ideal numerator(const alg_ideal &);



inline alg_number
reduce(alg_ideal & c, const alg_ideal & a)
{
	alg_number d;

	c.assign(a);
	c.reduce(d);
	return d;
}



inline alg_ideal
operator + (const alg_ideal & a, const alg_ideal & b)
{
	alg_ideal c(a.which_base());


	add(c, a, b);
	return c;
}



inline alg_ideal
operator & (const alg_ideal & a, const alg_ideal & b)
{
	alg_ideal c(a.which_base());


	intersect(c, a, b);
	return c;
}



inline alg_ideal
operator * (const alg_ideal & a, const alg_ideal & b)
{
	alg_ideal c(a.which_base());


	multiply(c, a, b);
	return c;
}



inline alg_ideal
operator * (const alg_ideal &a, const bigint &b)
{
	alg_ideal c(a.which_base());


	multiply(c, static_cast<const module &>(a), b);
	return c;
}



inline alg_ideal
operator * (const bigint &b, const alg_ideal &a)
{
	alg_ideal c(a.which_base());


	multiply(c, static_cast<const module&>(a), b);
	return c;
}



inline alg_ideal
operator * (const alg_ideal &a, const alg_number &b)
{
	alg_ideal c(a.which_base());


	multiply(c, a, b);
	return c;
}



inline alg_ideal
operator * (const alg_number &b, const alg_ideal &a)
{
	alg_ideal c(a.which_base());


	multiply(c, a, b);
	return c;
}



inline alg_ideal
operator / (const alg_ideal & a, const module & b)
{
	alg_ideal c(a.which_base());


	divide(c, a, b);
	return c;
}



inline alg_ideal
operator / (const alg_ideal & a, const bigint & b)
{
	alg_ideal c(a.which_base());


	divide(c, static_cast<const module &>(a), b);
	return c;
}



inline alg_ideal
operator / (const alg_ideal & a, const alg_number & b)
{
	alg_ideal c(a.which_base());


	divide(c, a, b);
	return c;
}



inline alg_ideal
operator % (const alg_ideal &a, const bigint &p)
{
	alg_ideal c(a.which_base());


	remainder(c, a, p);
	return c;
}



inline alg_ideal
numerator(const alg_ideal & a)
{
	return alg_ideal(a.coeff_matrix(), bigint(1), a.which_base());
}


inline void
module::assign_whole_order ()	// Remains in the same order
{
	assign(alg_ideal_cast(order(O)));
}



#undef alg_ideal_cast


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_ALG_IDEAL

#include	"LiDIA/specialization/alg_ideal.special"



#endif	// LIDIA_ALG_NUMBER_H_GUARD_
