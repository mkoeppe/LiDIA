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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


// Description  : Implementation of multivariate polynomials over GF(2^n)


#ifndef LIDIA_MV_POLY_H_GUARD_
#define LIDIA_MV_POLY_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_GF2N_H_GUARD_
# include	"LiDIA/gf2n.h"
#endif
#ifndef LIDIA_MV_TERM_H_GUARD_
# include	"LiDIA/mv_term.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



extern "C" {

typedef struct mvle
{
	mv_term term;
	struct mvle * next;
} mv_list_entry;

}



class mv_poly
{
	//
	// this class represents a multivariate polynomial
	// as a list of objects of type mv_term
	//
	// a freelist is maintained to optimize memory-allocation
	//
	// a term has the form coeff * var
	// where var is of form x_0^i_0 ... x_(i.bit_length-1)^(i.bit_length-1)
	// The polynomial 0 is stored as empty list.
	//

private:

	mv_list_entry *first_term; // pointer to the first element of list
	mv_list_entry *last_term; // pointer to the last element of list
	int len; // number of terms

	static mv_list_entry *free_list; // the free list pointer
	static int size_free_list;
	static int refc; // reference counter

public:

	//
	// constructors and destructor
	//

	mv_poly()    // = 0
	{
		first_term = last_term = NULL;
		len = 0;
		mv_poly::refc ++;
	}

	mv_poly(const gf2n & g)
	{
		first_term = mv_poly::get_new();
		first_term->next = NULL;
		first_term->term.assign(mv_term(g));
		last_term = first_term;
		mv_poly::refc ++;
		len = 1;
	}

	mv_poly(const mv_term & t)
	{
		first_term = mv_poly::get_new();
		first_term->next = NULL;
		first_term->term.assign(t);
		last_term = first_term;
		mv_poly::refc ++;
		len = 1;
	}


	mv_poly(const mv_poly & p)
	{
		mv_poly::refc ++;
		first_term = last_term = NULL;
		*this = p;
	}

	~mv_poly()
	{
		mv_poly::put_in_free(*this);

		if (mv_poly::refc == 1)
			mv_poly::delete_free_list();

		mv_poly::refc --;
	}

	//
	// assignments
	//

	mv_poly & operator = (const mv_poly & p);
	mv_poly & operator = (const mv_term & t);

	friend void swap(mv_poly &, mv_poly &);

	void assign_zero()
	{
		mv_poly::put_in_free(*this);
		len = 0;
	}

	void assign_one();

	void assign(const mv_poly & mp)
	{
		*this = mp;
	}

	void assign(const mv_term & t)
	{
		*this = t;
	}

	const mv_list_entry *  get_first() const
	{
		return first_term;
	}

	const mv_list_entry *  get_last() const
	{
		return last_term;
	}

	int max_index_variable() const;


	//
	// comparisons
	//

	friend bool operator == (const mv_poly & a, const mv_poly & b);

	bool is_zero() const
	{
		return (first_term == NULL);
	}

	bool is_one() const;

	// returns iff mv_poly is non-zero (!!) constant

	bool is_const() const;

	//
	// operator overloading
	//

	int length() const
	{
		return len;
	}

	mv_poly & operator += (const mv_poly & a);

	mv_poly & operator += (const mv_term & t)
	{
		mv_poly temp(t);

		(*this) += temp;
		return *this;
	}

	mv_poly & operator *= (const mv_term & t);

	mv_poly & operator *= (const mv_poly & a);

	//
	// Procedural versions
	//

	friend void add(mv_poly & c, const mv_poly & a, const mv_poly & b);

	friend void multiply(mv_poly & c, const mv_poly & a, const mv_poly & b);

	friend void multiply(mv_poly & c, const gf2n & a, const mv_poly & b);

	friend void square(mv_poly & b, const mv_poly & a);

	friend void sqrt(mv_poly & b, const mv_poly & a);

	//
	// general functions
	//


	// solve x_k=q, where q is a polynomial, which contains no x_k;
	// the function determines and returns a suitable k

	friend bool solve_x_k (mv_poly & q, lidia_size_t & k);

	// solve x_k=q, where q contains no x_k; k is chosen by the user

	friend bool solve_x_k_fixed(mv_poly & q, lidia_size_t k);

	// the k-th Variable is substituted by the polynomial p and stored in q

	friend void substitute(mv_poly & q, const mv_poly & p, lidia_size_t k);

	// the k-th Variable of q is substituted by the polynomial n/d, the numerator
	// of the result is returned in q

	friend void substitute(mv_poly & q, const mv_poly & n, const mv_poly & d,
			       lidia_size_t k);

	friend void substitute(mv_poly & qn, mv_poly & qd, const mv_poly & n,
			       const mv_poly & d, lidia_size_t k);

	// all terms of t are evaluated, X_i is set to the i-th bit of c

	friend void evaluate(gf2n & a, const mv_poly & t, const bigint & c);

	friend void evaluate(mv_poly &, const bigint &, const bigint &);

	// split *this in two parts: terms not containing X_k are kept in *this,
	// terms with X_k are removed from *this and returned in q_k, X_k is removed
	// from all these terms.

	void split_x_k (mv_poly & q_k, lidia_size_t k);

	// the function returns k if *this has the from coeff_1 + coeff_2 * X_k
	// where coeff_1 is allowed to be zero. Otherwise it returns -1.

	lidia_size_t linear_poly_in_only_one_var() const;

	// function returns true iff there is some term with variable X_k

	bool has_var_k(lidia_size_t k) const;

	// checks whether *this has a linear term X_k and X_k
	// does not occur in any other list term;
	// if k=-1, the function determines some value for k.
	// If no such term exists, the function returns NULL;
	// otherwise it returns a pointer to the term before the
	// linear term, or a pointer to the first term of the list

	mv_list_entry* has_one_linear_term(lidia_size_t & k);

	// T is a mv_poly whose terms are exactly the trace of the corresponding
	// terms in *this (variables are not changed)

	void trace_computation(mv_poly & T);

	//
	//   output
	//

	friend std::ostream & operator << (std::ostream & out, const mv_poly & a);



private:

	//
	// the following functions are for the handling of the free_list
	//

	static void put_in_free (mv_poly& p);
	static void put_in_free (mv_list_entry& t);
	static mv_list_entry * get_new(unsigned int anz = 1);

public:

	static void delete_free_list(unsigned int anz = 2147483647);

	//
	// Compare-Method needed for qsort
	//

private:

	// sort a polynomial with counter many terms, by using the library-function
	// qsort()

	void sort(lidia_size_t);

	// collect terms with the same variables in one term

	void clean();

};



inline void add(mv_poly & c, const mv_term t, const mv_poly & b)
{
	mv_poly temp(t);

	add(c, b, temp);
}



inline void multiply(mv_poly & c, const mv_poly & b, const gf2n & a)
{
	multiply(c, a, b);
}



inline void multiply(mv_poly & c, const mv_term & t, const mv_poly & b)
{
	mv_poly temp(t);

	multiply(c, b, temp);
}



inline void multiply(mv_poly & c, const mv_poly & b, const mv_term & t)
{
	mv_poly temp(t);

	multiply(c, b, temp);
}



inline mv_poly operator - (const mv_poly & a)
{
	mv_poly c(a);

	return c;
}



inline mv_poly operator + (const mv_poly & a, const mv_poly & b)
{
	mv_poly c;

	add(c, a, b);
	return c;
}



inline mv_poly operator - (const mv_poly & a, const mv_poly & b)
{
	mv_poly c;

	add(c, a, b);
	return c;
}



inline mv_poly operator * (const mv_poly & a, const mv_poly & b)
{
	mv_poly c;

	multiply(c, a, b);
	return c;
}



inline bool operator != (const mv_poly & a, const mv_poly & b)
{
	return ! (a == b);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MV_POLY_H_GUARD_
