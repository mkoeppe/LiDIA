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
//	Author	: Andrea Rau (AR), Robert Carls (RC), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


// Description  : Class that stores one term of a multi-variate
//                binary polynomial.


#ifndef LIDIA_MV_TERM_H_GUARD_
#define LIDIA_MV_TERM_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_GF2N_H_GUARD_
# include	"LiDIA/gf2n.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class mv_term
{
	//
	// this class represents a term from a multivariate polynomial with
	// binary variables
	//
	// a term has the form coeff * var
	// where var is of the form x_0^i_0 ... x_(i.bit_length-1)^(i.bit_length-1)
	//

	gf2n coeff; // Coefficient
	bigint var; // Variable-list, will be stored in dense representation

public:

	//
	// constructors and destructor
	//

	mv_term()
	{
		coeff.assign_zero();
		var.assign_zero();
	}

	mv_term(const gf2n & c)
	{
		coeff.assign(c);
		var.assign_zero();
	}

	mv_term(const gf2n & c, const bigint & v)
	{
		coeff.assign(c);
		var.assign(v);
	}

	mv_term(const mv_term & t)  : coeff(t.coeff), var(t.var)
	{ }

	~mv_term()
	{ }

	//
	// access
	//

	const gf2n & get_coeff() const
	{
		return coeff;
	}

	const bigint & get_var() const
	{
		return var;
	}

	//
	// assignments
	//

	void assign_zero()
	{
		coeff.assign_zero();
		var.assign_zero();
	}

	void assign_one()
	{
		coeff.assign_one();
		var.assign_zero();
	}

	void assign(const gf2n & c)
	{
		coeff.assign(c);
		var.assign_zero();
	}

	void assign(const mv_term & c)
	{
		coeff.assign(c.coeff);
		var.assign(c.var);
	}

	void assign_coeff(const gf2n & c)
	{
		coeff.assign(c);
	}

	void assign_var(const bigint & v)
	{
		var.assign(v);
	}

	void assign(const gf2n & c, const bigint & v)
	{
		coeff.assign(c);
		var.assign(v);
	}

	mv_term & operator = (const gf2n & c)
	{
		assign(c);
		return *this;
	}

	mv_term & operator = (const mv_term & t)
	{
		assign(t.coeff, t.var);
		return *this;
	}

	// variable and coefficient are equal

	bool all_equal(const mv_term & a) const
	{
		return (get_var() == a.get_var()  &&
			get_coeff() == a.get_coeff());
	}

	bool is_zero() const
	{
		return coeff.is_zero();
	}

	bool is_one() const
	{
		return coeff.is_one();
	}

	bool is_const() const
	{
		return var.is_zero();
	}

	// returns true if X_k is variable of term

	bool has_var_k(unsigned int k) const
	{
		return var.bit(k);
	}


	bool is_linear() const; // returns true if term is of form C*X_k
	// for some gf2n C and integer k

	bool is_linear(lidia_size_t & k) const; // same as bool is_linear() plus
	// sets k if term is of form C*X_k
	// for a certain k, otherwise  k = -1


	friend void add(mv_term & c, const mv_term & a, const mv_term & b);
	friend void multiply(mv_term & c, const gf2n & a, const mv_term & b);
	friend void multiply(mv_term & c, const mv_term & a, const mv_term & b);
	friend void square(mv_term & b, const mv_term & a);
	friend void sqrt(mv_term & b, const mv_term & a);

	//
	// functions
	//

	// in term t the variable X_i is substituted by term s

	friend void substitute(mv_term & t, const mv_term & s, unsigned int i);

	// in term t all the variables with index 1 bits in 'bits' are
	// substituted by the corresponding bits in c

	friend void substitute(mv_term & t, const bigint & c, const bigint & bits);

};



//
// comparisons
//

// only the variables are compared (we use lexicographical order)

inline bool
operator == (const mv_term & a, const mv_term & b)
{
	return (a.get_var() == b.get_var());
}



inline bool
operator != (const mv_term & a, const mv_term & b)
{
	return (a.get_var() != b.get_var());
}



inline bool
operator > (const mv_term & a, const mv_term & b)
{
	return (a.get_var() > b.get_var());
}



inline bool
operator >= (const mv_term & a, const mv_term & b)
{
	return (a.get_var() >= b.get_var());
}



inline bool
operator < (const mv_term & a, const mv_term & b)
{
	return (a.get_var() < b.get_var());
}



inline bool
operator <= (const mv_term & a, const mv_term & b)
{
	return (a.get_var() <= b.get_var());
}



//
// Procedural arithmetic
//

inline void
add(mv_term & c, const mv_term & a, const mv_term & b)
{
	if (a.var != b.var)
		lidia_error_handler("mv_term",
				    "void add::terms have different variables.");

	add(c.coeff, a.coeff, b.coeff);
	c.var.assign(a.var);
}



inline void
subtract(mv_term & c, const mv_term & a, const mv_term & b)
{
	add(c, a, b);
}



inline void
multiply(mv_term & c, const gf2n & a, const mv_term & b)
{
	if (a.is_zero() || b.is_zero())
		c.assign_zero();
	else {
		multiply(c.coeff, a, b.coeff);
		c.var.assign(b.var);
	}
}



inline void
multiply(mv_term & c, const mv_term & a, const gf2n & b)
{
	multiply(c, b, a);
}



inline void
multiply(mv_term & c, const mv_term & a, const mv_term & b)
{
	if (a.is_zero() || b.is_zero())
		c.assign_zero();
	else {
		multiply(c.coeff, a.coeff, b.coeff);
		bitwise_or(c.var, a.var, b.var);
	}
}



inline void
square(mv_term & b, const mv_term & a)
{
	square(b.coeff, a.coeff);
	b.var.assign(a.var);
}



inline void
sqrt(mv_term & b, const mv_term & a)
{
	sqrt(b.coeff, a.coeff);
	b.var.assign(a.var);
}



//
// operator overloading
//

inline mv_term
operator + (const mv_term & a, const mv_term & b)
{
	mv_term c;

	add(c, a, b);
	return c;
}



inline mv_term
operator - (const mv_term & a, const mv_term & b)
{
	mv_term c;

	add(c, a, b);
	return c;
}



inline mv_term
operator * (const gf2n & a, const mv_term & b)
{
	mv_term c;

	multiply(c, a, b);
	return c;
}



inline mv_term
operator * (const mv_term & a, const mv_term & b)
{
	mv_term c;

	multiply(c, a, b);
	return c;
}



inline mv_term &
operator += (mv_term & a, const mv_term & b)
{
	add(a, a, b);
	return a;
}



inline mv_term &
operator -= (mv_term & a, const mv_term & b)
{
	add(a, a, b);
	return a;
}



inline mv_term &
operator *= (mv_term & a, const mv_term & b)
{
	multiply(a, a, b);
	return a;
}



// term t will be evaluated by the values for the variables that are in c

inline void
evaluate(gf2n & a, const mv_term & t, const bigint & c)
{
	const bigint& t_var(t.get_var());

	if ((t_var & c) == t_var)
		a.assign(t.get_coeff());
	else
		a.assign_zero();
}



//
// input / output
//

std::ostream & operator << (std::ostream & out, const mv_term & a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MV_TERM_H_GUARD_
