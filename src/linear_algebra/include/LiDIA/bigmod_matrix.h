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
//	Author	: Patrick Theobald (PT), Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGMOD_MATRIX_H_GUARD_
#define LIDIA_BIGMOD_MATRIX_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigmod_matrix : public bigint_matrix
{
public:

	//********************************************
	//** Constructors ****************************
	//********************************************

	bigmod_matrix();
	bigmod_matrix(lidia_size_t, lidia_size_t);
	bigmod_matrix(lidia_size_t, lidia_size_t, const bigint &);
	bigmod_matrix(const bigmod_matrix &);
	bigmod_matrix(const base_matrix< bigint > &, const bigint & mod = 0);

	//********************************************
	//** destructor ******************************
	//********************************************

	~bigmod_matrix();

	//
	// Input / Output
	//

	friend std::ostream & operator << (std::ostream &, const bigmod_matrix &);
	friend std::istream & operator >> (std::istream &, bigmod_matrix &);

	//
	// stream handling
	//

	void write_to_stream(std::ostream &) const;

	void read_from_stream(std::istream &);

	//
	// handling the modulus
	//

	const bigint & get_modulus() const
	{
		return p;
	}
	const bigint & get_modulus(const bigmod_matrix & A)
	{
		return A.get_modulus();
	}

	void set_modulus(const bigint & mod)
	{
		p.assign(mod);
	}

	void set_modulus(bigmod_matrix & A, const bigint & mod)
	{
		A.set_modulus(mod);
	}

	void reduce(const bigint &);
	void reduce(bigmod_matrix & A, const bigint & mod)
	{
		A.reduce(mod);
	}

	void lift(const bigint &);
	void lift(bigmod_matrix & A, const bigint & mod)
	{
		A.lift(mod);
	}

	//
	// access functions
	//

	void sto(lidia_size_t, lidia_size_t, const bigint &);
	void sto_column(const bigint *, lidia_size_t,
			lidia_size_t, lidia_size_t from = 0);
	void sto_row(const bigint *, lidia_size_t,
		     lidia_size_t, lidia_size_t from = 0);

	//
	// access functions using vectors
	//

	void sto_column_vector(const base_vector< bigint > &, lidia_size_t,
			       lidia_size_t, lidia_size_t from = 0);
	void sto_row_vector(const base_vector< bigint > &, lidia_size_t,
			    lidia_size_t, lidia_size_t from = 0);


	//
	// split functions
	//

	void split(bigmod_matrix &, bigmod_matrix &, bigmod_matrix &, bigmod_matrix &) const;
	void split_h(bigmod_matrix &, bigmod_matrix &) const;
	void split_v(bigmod_matrix &, bigmod_matrix &) const;

	//
	// compose functions
	//

	void compose(const bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
	void compose_h(const bigmod_matrix &, const bigmod_matrix &);
	void compose_v(const bigmod_matrix &, const bigmod_matrix &);

	//
	// exchange functions / swap functions
	//

	friend void swap(bigmod_matrix &, bigmod_matrix &);

	//
	// BEGIN: matrix arithmetic
	//

	//
	// procedures
	//

	friend void add(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
	friend void add(bigmod_matrix &, const bigmod_matrix &, const bigint &);
	friend void subtract(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
	friend void subtract(bigmod_matrix &, const bigmod_matrix &, const bigint &);
	friend void negate(bigmod_matrix &, const bigmod_matrix &);
	friend void multiply(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
	friend void multiply(bigmod_matrix &, const bigmod_matrix &, const bigint &);
	friend void multiply(bigint *, const bigmod_matrix &, const bigint *);
	friend bool divide(bigmod_matrix &, const bigmod_matrix &, const bigint &);


	//
	// operators
	//

	bigmod_matrix operator - () const;
	bigmod_matrix operator + (const bigmod_matrix &) const;
	bigmod_matrix operator - (const bigmod_matrix &) const;
	bigmod_matrix operator * (const bigmod_matrix &) const;
	bigmod_matrix operator + (const bigint &) const;
	bigmod_matrix operator - (const bigint &) const;
	bigmod_matrix operator * (const bigint &) const;
	bigint * operator * (const bigint *) const;

	bigmod_matrix & operator += (const bigmod_matrix &);
	bigmod_matrix & operator += (const bigint &);
	bigmod_matrix & operator -= (const bigmod_matrix &);
	bigmod_matrix & operator -= (const bigint &);
	bigmod_matrix & operator *= (const bigmod_matrix &);
	bigmod_matrix & operator *= (const bigint &);


	//********************************************
	//** assign operator *************************
	//********************************************

	bigmod_matrix & operator = (const bigmod_matrix &);
	void assign(const bigmod_matrix &);
	friend void assign(bigmod_matrix &, const bigmod_matrix &);

	void assign(const base_matrix< bigint > &, const bigint & mod = 0);
	friend void assign(bigmod_matrix &,
			   const base_matrix< bigint > &, const bigint & mod);

	//********************************************
	//** boolean operators ***********************
	//********************************************

	bool operator == (const bigmod_matrix &) const;
	bool equal(const bigmod_matrix &) const;
	friend bool equal(const bigmod_matrix &, const bigmod_matrix &);

	bool operator != (const bigmod_matrix &) const;
	bool unequal(const bigmod_matrix &) const;
	friend bool unequal(const bigmod_matrix &, const bigmod_matrix &);

	//
	// randomize
	//

	void randomize(const bigint &);

	//********************************************
	//** diag function ***************************
	//********************************************

	void diag(const bigint &, const bigint &);
	friend void diag(bigmod_matrix &, const bigint &, const bigint &);


	//********************************************
	//** transpose function **********************
	//********************************************

	bigmod_matrix trans() const;
	friend bigmod_matrix trans(const bigmod_matrix &);


	//********************************************
	//** regular expansion ***********************
	//********************************************

	void regexpansion(const lidia_size_t *);
	friend void regexpansion(bigmod_matrix &, const lidia_size_t *);

	// LINEAR ALGEBRA - METHOD I:
	// We are satisfied with finding a factor.

	//
	// rank
	// (not implemented)
	//

	//lidia_size_t rank(bigint &) const;
	//friend lidia_size_t rank(const bigmod_matrix &, bigint &);

	//
	// rank and linearly independent rows
	//

	lidia_size_t *lininr(bigint &) const;
	friend lidia_size_t *lininr(const bigmod_matrix &, bigint &);

	//
	// rank and linearly independent columns
	// (not implemented)
	//

	//lidia_size_t *lininc(bigint &) const;
	//friend lidia_size_t *lininc(const bigmod_matrix &, bigint &);

	//
	// inverse and adjoint matrix
	//

	void inv(const bigmod_matrix &, bigint &);
	friend bigmod_matrix inv(const bigmod_matrix &, bigint &);

	void adj(const bigmod_matrix &, bigint &);
	friend bigmod_matrix adj(const bigmod_matrix &, bigint &);

	//
	// determinant
	//

	void det(bigint &, bigint &) const;
	bigint det(bigint &) const;
	friend bigint det(const bigmod_matrix &, bigint &);

	//
	// characteristic polynomial
	//

	bigint *charpoly(bigint &) const;
	friend bigint *charpoly(const bigmod_matrix &, bigint &);

	//********************************************
	//** special forms ***************************
	//********************************************

	friend int stf(bigmod_matrix &, bigmod_matrix &, bigint &);
	void stf(int, bigint &);
	friend bigint stf(bigmod_matrix &, bigint &);

	friend void hbf(bigmod_matrix &, bigint &);

	//
	// Kernel
	//

	void kernel(const bigmod_matrix &, bigint &);
	friend bigmod_matrix kernel(const bigmod_matrix &, bigint &);

	//
	// regular Invimage
	// (not implemented!!)
	//

	//void reginvimage(const bigmod_matrix &, const bigmod_matrix &, bigint &);
	//friend bigmod_matrix reginvimage(const bigmod_matrix &, const bigmod_matrix &, bigint &);

	//
	// Image
	//

	void image(const bigmod_matrix &, bigint &);
	friend bigmod_matrix image(const bigmod_matrix &, bigint &);

	//
	// solve
	//

	void solve(const bigmod_matrix &, const bigint *, bigint &);
	friend bigmod_matrix solve(const bigmod_matrix &, const bigint *, bigint &);
	void solve(const bigmod_matrix &, const base_vector< bigint > &, bigint &);
	friend bigmod_matrix solve(const bigmod_matrix &,
				   const base_vector< bigint > &, bigint &);


	// LINEAR ALGEBRA - METHOD II:
	// We are not satisfied with finding a factor.

#ifdef MNF_DEBUG
	friend void mnf_debug_info(const bigmod_matrix &, const bigmod_matrix &,
				   const bigmod_matrix &, const bigmod_matrix &,
				   const bigmod_matrix &, lidia_size_t,
				   const bigint &);
#endif

	//
	// Exponent
	//

	void exponent(bigint &) const;
	bigint exponent() const;
	friend bigint exponent(const bigmod_matrix &);

	//
	// Kernel
	//

	void kernel(const bigmod_matrix &);
	friend bigmod_matrix kernel(const bigmod_matrix &);

	//
	// Image
	//

	void image(const bigmod_matrix &);
	friend bigmod_matrix image(const bigmod_matrix &);

	void unique_image(const bigmod_matrix &);
	friend bigmod_matrix unique_image(const bigmod_matrix &);

private:
	bigint p; // Internal modulus for matrix

	//********************************************
	//** base ************************************
	//********************************************

	void add_mod(bigint &, const bigint &, const bigint &) const;
	void sub_mod(bigint &, const bigint &, const bigint &) const;
	void div_mod(bigint &, bigint &, const bigint &, const bigint &) const;
	void mult_mod(bigint &, const bigint &, const bigint &) const;
	void inv_mod(bigint &, bigint &, const bigint &) const;
};

// friend functions

//
// Input / Output
//

std::ostream & operator << (std::ostream &, const bigmod_matrix &);
std::istream & operator >> (std::istream &, bigmod_matrix &);

//
// exchange functions / swap functions
//

void swap(bigmod_matrix &, bigmod_matrix &);

//
// BEGIN: matrix arithmetic
//

//
// procedures
//

void add(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
void add(bigmod_matrix &, const bigmod_matrix &, const bigint &);
void subtract(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
void subtract(bigmod_matrix &, const bigmod_matrix &, const bigint &);
void negate(bigmod_matrix &, const bigmod_matrix &);
void multiply(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &);
void multiply(bigmod_matrix &, const bigmod_matrix &, const bigint &);
void multiply(bigint *, const bigmod_matrix &, const bigint *);
bool divide(bigmod_matrix &, const bigmod_matrix &, const bigint &);


//********************************************
//** assign operator *************************
//********************************************

void assign(bigmod_matrix &, const bigmod_matrix &);
void assign(bigmod_matrix &,
		const base_matrix< bigint > &, const bigint & mod = 0);

//********************************************
//** boolean operators ***********************
//********************************************

bool equal(const bigmod_matrix &, const bigmod_matrix &);
bool unequal(const bigmod_matrix &, const bigmod_matrix &);

//********************************************
//** diag function ***************************
//********************************************

void diag(bigmod_matrix &, const bigint &, const bigint &);


//********************************************
//** transpose function **********************
//********************************************

bigmod_matrix trans(const bigmod_matrix &);


//********************************************
//** regular expansion ***********************
//********************************************

void regexpansion(bigmod_matrix &, const lidia_size_t *);

//
// rank and linearly independent rows
//

lidia_size_t *lininr(const bigmod_matrix &, bigint &);

//
// inverse and adjoint matrix
//

bigmod_matrix inv(const bigmod_matrix &, bigint &);
bigmod_matrix adj(const bigmod_matrix &, bigint &);

//
// determinant
//

bigint det(const bigmod_matrix &, bigint &);

//
// characteristic polynomial
//

bigint *charpoly(const bigmod_matrix &, bigint &);

//********************************************
//** special forms ***************************
//********************************************

int stf(bigmod_matrix &, bigmod_matrix &, bigint &);
bigint stf(bigmod_matrix &, bigint &);
void hbf(bigmod_matrix &, bigint &);

//
// Kernel
//

bigmod_matrix kernel(const bigmod_matrix &, bigint &);

//
// Image
//

bigmod_matrix image(const bigmod_matrix &, bigint &);

//
// solve
//

bigmod_matrix solve(const bigmod_matrix &, const bigint *, bigint &);
bigmod_matrix solve(const bigmod_matrix &,
			const base_vector< bigint > &, bigint &);


// LINEAR ALGEBRA - METHOD II:
// We are not satisfied with finding a factor.

#ifdef MNF_DEBUG
void mnf_debug_info(const bigmod_matrix &, const bigmod_matrix &,
			const bigmod_matrix &, const bigmod_matrix &,
			const bigmod_matrix &, lidia_size_t,
			const bigint &);
#endif

//
// Exponent
//

bigint exponent(const bigmod_matrix &);

//
// Kernel
//

bigmod_matrix kernel(const bigmod_matrix &);

//
// Image
//

bigmod_matrix image(const bigmod_matrix &);
bigmod_matrix unique_image(const bigmod_matrix &);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGMOD_MATRIX_H_GUARD_
