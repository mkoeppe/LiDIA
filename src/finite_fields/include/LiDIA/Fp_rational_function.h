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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================

//	Description: This class realizes rational functions and
//              rational functions modulo some modulus which can
//              be either given as Fp_polynomial or as Fp_poly_modulus.


#ifndef LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_
#define LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_FP_POLY_MODULUS_H_GUARD_
# include	"LiDIA/Fp_poly_modulus.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class Fp_rational_function
{
private:

	Fp_polynomial *num;
	Fp_polynomial *den;

public:

	Fp_rational_function();
	Fp_rational_function(const bigint & p);
	Fp_rational_function(const Fp_polynomial &);
	Fp_rational_function(const Fp_polynomial & , const Fp_polynomial &);
	Fp_rational_function(const Fp_rational_function &);
	~Fp_rational_function();


	//******************************************************************
	// some simple assignment functions   etc.
	//******************************************************************

	void kill();

	void set_modulus(const bigint & p);
	const bigint & modulus() const;

	lidia_size_t degree_numerator() const;
	lidia_size_t degree_denominator() const;

	void get_coefficient_numerator(bigint & a, lidia_size_t i) const;
	void get_coefficient_denominator(bigint & a, lidia_size_t i) const;

	void set_coefficient_numerator(const bigint& a, lidia_size_t i);
	void set_coefficient_denominator(const bigint& a, lidia_size_t i);

	void set_coefficient_numerator(lidia_size_t i);
	void set_coefficient_denominator(lidia_size_t i);

	const bigint& lead_coefficient_numerator() const;
	const bigint& lead_coefficient_denominator() const;

	const bigint& const_term_numerator() const;
	const bigint& const_term_denominator() const;


	//*************************************************************
	// Assignments
	//*************************************************************

	Fp_rational_function & operator = (const Fp_rational_function &);
	Fp_rational_function & operator = (const Fp_polynomial &);

	void assign(const Fp_rational_function & f);

	void assign(const Fp_polynomial & f);

	void assign(const Fp_polynomial & f, const Fp_polynomial & g);

	void assign_numerator(const Fp_polynomial &a);

	void assign_denominator (const Fp_polynomial &a);

	void assign_zero();
	void assign_one();
	void assign_x();

	void randomize(lidia_size_t deg_num, lidia_size_t deg_denom = 0);

	Fp_polynomial & numerator();
	Fp_polynomial & denominator();

	const  Fp_polynomial & numerator() const;
	const  Fp_polynomial & denominator() const;

	bigint operator()(const bigint & a) const;



	//******* Comparisons ***************************************

	friend bool operator == (const Fp_rational_function &,
				 const Fp_rational_function &);
	friend bool operator != (const Fp_rational_function & a,
				 const Fp_rational_function & b);

	friend bool equal_mod(const Fp_rational_function &, const Fp_rational_function &,
			      const Fp_polynomial &);

	friend bool equal(const Fp_rational_function &, const Fp_rational_function &,
			  const Fp_poly_modulus &);

	bool is_zero() const;
	bool is_one() const;


	//***************************************************************
	// procedural versions for arithmetics
	//***************************************************************


#ifndef HEADBANGER

	void negate()
	{
		num->negate();
	}

	friend void negate(Fp_rational_function & c, const Fp_rational_function & a)
	{
		c.assign(a);
		c.negate();
	}

	friend void add(Fp_rational_function &, const Fp_rational_function &,
			const Fp_rational_function &);

	friend void subtract(Fp_rational_function &, const Fp_rational_function &,
			     const Fp_rational_function &);

	friend void multiply(Fp_rational_function & c, const Fp_rational_function & a,
			     const Fp_rational_function & b)
	{
		multiply(*c.num, *a.num, *b.num);
		multiply(*c.den, *a.den, *b.den);
	}

	friend void multiply (Fp_rational_function &x, const Fp_rational_function &a,
			      const Fp_polynomial & b)
	{
		multiply(*x.num, *a.num, b);
		x.den->assign(*a.den);
	}

	friend void multiply (Fp_rational_function &x, const Fp_polynomial &a,
			      const Fp_rational_function & b)
	{
		multiply(*x.num, a, *b.num);
		x.den->assign(*b.den);
	}

	friend void multiply(Fp_rational_function & x, const bigint & a,
			     const Fp_rational_function & b)
	{
		multiply_by_scalar(*x.num, *b.num, a);
		x.den->assign(*b.den);
	}

	friend void multiply(Fp_rational_function & c, const Fp_rational_function & a,
			     const bigint & b)
	{
		multiply_by_scalar(*c.num, *a.num, b);
		c.den->assign(*a.den);
	}

	void multiply_by_2()
	{
		add(*num, *num, *num);
	}

	void normalize()
	{
		Fp_polynomial c(gcd(*num, *den)); divide(*num, *num, c);
		divide(*den, *den, c);
	}

	friend void square(Fp_rational_function & x, const Fp_rational_function & a)
	{
		square(*x.num, *a.num);
		square(*x.den, *a.den);
	}

#ifndef HEADBANGER
	friend void div_rem(Fp_polynomial & q, Fp_rational_function & f)
	{
		Fp_polynomial rr;

		div_rem(q, rr, *f.num, *f.den);
		f.num->assign(rr);
	}

	friend void divide(Fp_rational_function & q, const Fp_rational_function & a,
			   const Fp_rational_function & b);

	friend void divide(Fp_rational_function & q, const Fp_polynomial & a,
			   const Fp_rational_function& b);

	friend void divide(Fp_rational_function& q, const Fp_rational_function& a,
			   const Fp_polynomial & b)
	{
		q.den->assign(*a.den);
		multiply(*q.num, *a.num, b);
	}
#endif

	void invert();

	friend void invert(Fp_rational_function & c, const Fp_rational_function & a);

	void reduce(const Fp_polynomial & f)
	{
		remainder(*num, *num, f);
		remainder(*den, *den, f);
	}

	void reduce(const Fp_poly_modulus & f)
	{
		remainder(*num, *num, f.modulus());
		remainder(*den, *den, f.modulus());
	}

	//***************************************************************
	// operators
	//***************************************************************

	Fp_rational_function & operator += (const Fp_rational_function &a)
	{
		add(*this, *this, a);
		return *this;
	}

	Fp_rational_function & operator -= (const Fp_rational_function &a)
	{
		subtract(*this, *this, a);
		return *this;
	}

	Fp_rational_function & operator *= (const Fp_rational_function &a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	Fp_rational_function & operator /= (const Fp_rational_function &a)
	{
		divide(*this, *this, a);
		return *this;
	}


	//***************************************************************
	// Miscellaneaous
	//***************************************************************

	friend void swap(Fp_rational_function& x, Fp_rational_function& y)
	{
		swap(*x.num, *y.num);
		swap(*x.den, *y.den);
	}

	friend void shift(Fp_rational_function & c, const Fp_rational_function & a,
			  lidia_size_t n)   // c = a * x^n
	{
		if (n == 0)
			c.assign(a);
		else
			if (n < 0) {
				c.num->assign(*a.num); shift_left(*c.den, *a.den, -n);
			}
			else {
				c.den->assign(*a.den); shift_left(*c.num, *a.num, n);
			}
	}

	friend void derivative(Fp_rational_function & c, const Fp_rational_function & a);
	// c = derivative of a

	friend Fp_rational_function derivative(const Fp_rational_function & a)
	{
		Fp_rational_function x;

		derivative(x, a);
		return x;
	}


	//***************************************************************
	// operators
	//***************************************************************

	friend Fp_rational_function operator - (const Fp_rational_function &a)
	{
		Fp_rational_function c(a);

		c.negate();
		return c;
	}

	friend Fp_rational_function operator + (const Fp_rational_function &a,
						const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		add(c, a, b);
		return c;
	}

	friend Fp_rational_function  operator + (const Fp_polynomial &a,
						 const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		add(c, Fp_rational_function(a), b);
		return c;
	}

	friend Fp_rational_function operator + (const Fp_rational_function &a,
						const Fp_polynomial &b)
	{
		Fp_rational_function c(a.modulus());

		add(c, Fp_rational_function(b), a);
		return c;
	}

	friend Fp_rational_function operator - (const Fp_rational_function &a,
						const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		subtract(c, a, b);
		return c;
	}

	friend Fp_rational_function  operator - (const Fp_polynomial &a,
						 const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		subtract(c, Fp_rational_function(a), b);
		return c;
	}

	friend Fp_rational_function operator - (const Fp_rational_function &a,
						const Fp_polynomial &b)
	{
		Fp_rational_function c(a.modulus());

		subtract(c, a, Fp_rational_function(b));
		return c;
	}

	friend Fp_rational_function  operator * (const Fp_rational_function &a,
						 const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		multiply(c, b, a);
		return c;
	}

	friend Fp_rational_function operator * (const Fp_polynomial &a,
						const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		multiply(c, b, a);
		return c;
	}

	friend Fp_rational_function operator * (const Fp_rational_function &a,
						const Fp_polynomial &b)
	{
		Fp_rational_function c(a.modulus());

		multiply(c, b, a);
		return c;
	}

	friend Fp_rational_function  operator * (const bigint & a,
						 const Fp_rational_function &b)
	{
		Fp_rational_function c(b.modulus());

		multiply(c, a, b);
		return c;
	}

	friend Fp_rational_function operator / (const Fp_rational_function &a,
						const Fp_rational_function &b)
	{
		Fp_rational_function c(a.modulus());

		divide(c, a, b);
		return c;
	}

	friend Fp_rational_function operator / (const Fp_polynomial &a,
						const Fp_rational_function &b)
	{
		Fp_rational_function c(b.modulus());

		divide(c, a, b);
		return c;
	}

	friend Fp_rational_function operator / (const Fp_rational_function &a,
						const Fp_polynomial &b)
	{
		Fp_rational_function c(a.modulus());

		divide(c, a, b);
		return c;
	}


	//***************************************************************
	//
	//	Modular Arithmetic without pre-conditioning
	//
	//***************************************************************

#ifndef HEADBANGER
	// arithmetic mod f.
	// all inputs and outputs are polynomials of degree less than deg(f).
	// Note that numerator and denominator are reduced modulo f, but
	// no conversion to Fp_polynomial is done.
	// ASSUMPTION: f is assumed monic, and deg(f) > 0.
	// NOTE: if you want to do many computations with a fixed f,
	//       use the Fp_poly_modulus data structure and associated
	//       routines below.

	friend void add_mod(Fp_rational_function& c, const Fp_rational_function & a,
			    const Fp_rational_function & b, const Fp_polynomial & f);

	friend void subtract_mod(Fp_rational_function& c, const Fp_rational_function & a,
				 const Fp_rational_function & b, const Fp_polynomial & f);

	friend void multiply_mod(Fp_rational_function& c, const Fp_rational_function & a,
				 const Fp_rational_function & b, const Fp_polynomial & f)
	{
		multiply_mod(*c.num, *a.num, *b.num, f);
		multiply_mod(*c.den, *a.den, *b.den, f);
	}

	friend void square_mod(Fp_rational_function & c, const Fp_rational_function & a,
			       const Fp_polynomial & f)
	{
		square_mod(*c.num, *a.num, f);
		square_mod(*c.den, *a.den, f);
	}


	friend void convert_mod(Fp_polynomial & c, const Fp_rational_function & a,
				const Fp_polynomial & f)
		// c = a % f, error if *a.den is not invertible
	{
		invert_mod(c, *a.den, f);
		multiply_mod(c, c, *a.num, f);
	}

	friend bool convert_mod_status(Fp_polynomial & c, const Fp_rational_function & a,
				       const Fp_polynomial& f)
		// if (a.den, f) = 1, returns 1 and sets c = a % f
		// otherwise, returns 0 and sets c = (a.den, f)
	{
		bool t = invert_mod_status(c, *a.den, f);
		if (t)
			multiply_mod(c, c, *a.num, f);
		return t;
	}

	friend void divide_mod(Fp_rational_function & q, const Fp_rational_function & a,
			       const Fp_rational_function & b, const Fp_polynomial & f)
	{
		Fp_rational_function h(b);

		h.invert();
		multiply_mod(q, a, h, f);
	}
#endif

	//***************************************************************
	// I/O
	//
	// two I/O formats:
	//
	// 1)   [a_n ... a_1 a_0] / [b_m ... b_0] mod p
	//
	// 2)   a_n*x^n + ... + a_1*x + a_0 / b_m*x^m + ... b_0 mod p
	// is self-explanatory.
	// The default is 2).
	//
	// On output, all coefficients will be integers between 0 and p-1,
	// and a_n, b_m not zero (the zero polynomial is [ ]).
	// On input, the coefficients are arbitrary integers which are
	// then reduced modulo p, and leading zeros stripped. The function
	// Fp_polynomial::read(std::istream&) can handle both formats; the decision
	// is made according to the first character being '[' or not.
	//
	//*********************************************************************

	friend std::istream& operator >> (std::istream& s, Fp_rational_function & a)
	{
		a.read(s);
		return s;
	}

	friend std::ostream& operator << (std::ostream & s, const Fp_rational_function & a)
	{
		a.pretty_print(s);
		return s;
	}

private:

	void pretty_read(std::istream &s); //reads format 2)
	void pretty_print(std::ostream &s) const; //writes format 2)

public:

	void read(std::istream &s = std::cin); // reads any format
	void print(std::ostream &s = std::cout) const; //writes format 1)

	//***************************************************************
	// now the Fp_poly_modulus functions for faster arithmetic
	//***************************************************************

	// If you need to do a lot of arithmetic modulo a fixed f,
	// build Fp_poly_modulus F for f.  This pre-computes information about f
	// that speeds up the computation a great deal.
	// f should be monic, and deg(f) > 0.

	friend void add(Fp_rational_function & x, const Fp_rational_function & a,
			const Fp_rational_function & b, const Fp_poly_modulus & F);

	friend void subtract(Fp_rational_function & x, const Fp_rational_function & a,
			     const Fp_rational_function & b, const Fp_poly_modulus & F);

	friend void multiply(Fp_rational_function & x, const Fp_rational_function & a,
			     const Fp_rational_function & b, const Fp_poly_modulus & F)
		// x = (a * b) % f, deg(a), deg(b) < f.degree()
	{
		multiply(*x.num, *a.num, *b.num, F);
		multiply(*x.den, *a.den, *b.den, F);
	}

	friend void square(Fp_rational_function& x, const Fp_rational_function& a,
			   const Fp_poly_modulus& F)
		// x = a^2 % f, deg(a) < f.degree()
	{
		square (*x.num, *a.num, F);
		square(*x.den, *a.den, F);
	}

	friend void divide(Fp_rational_function& x, const Fp_rational_function& a,
			   const Fp_rational_function & b, const Fp_poly_modulus& F)
	{
		Fp_rational_function h(b);

		h.invert();
		multiply(x, a, h, F);
	}

#endif	// HEADBANGER
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_
