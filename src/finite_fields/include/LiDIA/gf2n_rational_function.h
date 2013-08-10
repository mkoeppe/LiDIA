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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


// Description: This class realizes rational functions and
//              rational functions modulo some modulus which can
//              be either given as gf2n_polynomial or as poly_modulus.


#ifndef LIDIA_GF2N_RATIONAL_FUNCTION_H_GUARD_
#define LIDIA_GF2N_RATIONAL_FUNCTION_H_GUARD_


#ifndef LIDIA_GF2N_POLYNOMIAL_H_GUARD_
# include	"LiDIA/gf2n_polynomial.h"
#endif
#ifndef LIDIA_GF2N_POLY_MODULUS_H_GUARD_
# include	"LiDIA/gf2n_poly_modulus.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf2n_rational_function
{
private:

	gf2n_polynomial *num;
	gf2n_polynomial *den;

public:

	gf2n_rational_function();
	gf2n_rational_function(const gf2n_polynomial &);
	gf2n_rational_function(const gf2n_polynomial & , const gf2n_polynomial &);
	gf2n_rational_function(const gf2n_rational_function &);
	~gf2n_rational_function();


	//******************************************************************
	// some simple assignment functions   etc.
	//******************************************************************

	void kill()
	{
		num->kill();
		den->kill();
	}

	lidia_size_t degree_numerator() const;
	lidia_size_t degree_denominator() const;

	void get_coefficient_numerator(gf2n & a, lidia_size_t i) const;
	void get_coefficient_denominator(gf2n & a, lidia_size_t i) const;

	void set_coefficient_numerator(const gf2n& a, lidia_size_t i);
	void set_coefficient_denominator(const gf2n& a, lidia_size_t i);

	void set_coefficient_numerator(lidia_size_t i);
	void set_coefficient_denominator(lidia_size_t i);

	gf2n lead_coefficient_numerator() const;
	gf2n lead_coefficient_denominator() const;

	gf2n const_term_numerator() const;
	gf2n const_term_denominator() const;


	//************************************************************
	// Assignments
	//************************************************************

	gf2n_rational_function & operator = (const gf2n_rational_function &);
	gf2n_rational_function & operator = (const gf2n_polynomial &);

	void assign(const gf2n_rational_function & f);

	void assign(const gf2n_polynomial & f);

	void assign(const gf2n_polynomial & f, const gf2n_polynomial & g);

	void assign_numerator(const gf2n_polynomial &a);

	void assign_denominator (const gf2n_polynomial &a);

	void assign_zero();
	void assign_one();
	void assign_x();

	void randomize(lidia_size_t deg_num, lidia_size_t deg_denom = 0);

	gf2n_polynomial & numerator();
	gf2n_polynomial & denominator();

	const  gf2n_polynomial & numerator() const;
	const  gf2n_polynomial & denominator() const;

	gf2n operator()(const gf2n & a) const;



	// ******* Comparisons **************************************

	friend bool operator == (const gf2n_rational_function &,
				 const gf2n_rational_function &);
	friend bool operator != (const gf2n_rational_function & a,
				 const gf2n_rational_function & b);

	friend bool equal_mod(const gf2n_rational_function &,
			      const gf2n_rational_function &,
			      const gf2n_polynomial &);

	friend bool equal(const gf2n_rational_function &,
			  const gf2n_rational_function &,
			  gf2n_poly_modulus &);

	bool is_zero() const;
	bool is_one() const;


	//***************************************************************
	// procedural versions for arithmetics
	//***************************************************************


#ifndef HEADBANGER

	void negate()
	{  }

	friend void negate(gf2n_rational_function & c,
			   const gf2n_rational_function & a)
	{
		c.assign(a);
	}

	friend void add(gf2n_rational_function &, const gf2n_rational_function &,
			const gf2n_rational_function &);

	friend void add(gf2n_rational_function & c, const gf2n_rational_function & a,
			const gf2n & b)
	{
		add(*c.num, *a.num, b * (*a.den));
		c.den->assign(*a.den);
	}

	friend void add(gf2n_rational_function & c, const gf2n & b,
			const gf2n_rational_function & a)
	{
		add(*c.num, *a.num, b * (*a.den));
		c.den->assign(*a.den);
	}

	friend void subtract(gf2n_rational_function & c,
			     const gf2n_rational_function & a,
			     const gf2n_rational_function & b)
	{
		add(c, a, b);
	}

	friend void multiply(gf2n_rational_function & c,
			     const gf2n_rational_function & a,
			     const gf2n_rational_function & b)
	{
		multiply(*c.num, *a.num, *b.num);
		multiply(*c.den, *a.den, *b.den);
	}

	friend void multiply (gf2n_rational_function &x,
			      const gf2n_rational_function &a,
			      const gf2n_polynomial & b)
	{
		multiply(*x.num, *a.num, b);
		x.den->assign(*a.den);
	}

	friend void multiply (gf2n_rational_function &x, const gf2n_polynomial &a,
			      const gf2n_rational_function & b)
	{
		multiply(*x.num, a, *b.num);
		x.den->assign(*b.den);
	}

	friend void multiply(gf2n_rational_function & x, const gf2n & a,
			     const gf2n_rational_function & b)
	{
		multiply_by_scalar(*x.num, *b.num, a);
		x.den->assign(*b.den);
	}

	friend void multiply(gf2n_rational_function & c,
			     const gf2n_rational_function & a,
			     const gf2n & b)
	{
		multiply_by_scalar(*c.num, *a.num, b);
		c.den->assign(*a.den);
	}

	void multiply_by_2()
	{
		assign_zero();
	}

	void normalize()
	{
		gf2n_polynomial c;
		gcd(c, *num, *den); divide(*num, *num, c);
		divide(*den, *den, c);
	}

	friend void square(gf2n_rational_function & x,
			   const gf2n_rational_function & a)
	{
		square(*x.num, *a.num);
		square(*x.den, *a.den);
	}

#ifndef HEADBANGER
	friend void div_rem(gf2n_polynomial & q, gf2n_rational_function & f)
	{
		gf2n_polynomial rr;

		div_rem(q, rr, *f.num, *f.den);
		f.num->assign(rr);
	}

	friend void divide(gf2n_rational_function & q,
			   const gf2n_rational_function & a,
			   const gf2n_rational_function & b);

	friend void divide(gf2n_rational_function & q, const gf2n_polynomial & a,
			   const gf2n_rational_function& b);

	friend void divide(gf2n_rational_function& q, const gf2n_rational_function& a,
			   const gf2n_polynomial & b)
	{
		q.den->assign(*a.den);
		multiply(*q.num, *a.num, b);
	}
#endif	// HEADBANGER

	void invert();

	friend void invert(gf2n_rational_function & c,
			   const gf2n_rational_function & a);

	//***************************************************************
	// operators
	//***************************************************************

	gf2n_rational_function & operator += (const gf2n_rational_function &a)
	{
		add(*this, *this, a);
		return *this;
	}

	gf2n_rational_function & operator -= (const gf2n_rational_function &a)
	{
		add (*this, *this, a);
		return *this;
	}

	gf2n_rational_function & operator *= (const gf2n_rational_function &a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	gf2n_rational_function & operator /= (const gf2n_rational_function &a)
	{
		divide(*this, *this, a);
		return *this;
	}




	void reduce(const gf2n_polynomial & f)
	{
		remainder(*num, *num, f);
		remainder(*den, *den, f);
	}

	void reduce(gf2n_poly_modulus & f)
	{
		remainder(*num, *num, f.modulus());
		remainder(*den, *den, f.modulus());
	}


	//***************************************************************
	// Miscellaneaous
	//***************************************************************

	friend void swap(gf2n_rational_function& x, gf2n_rational_function& y)
	{
		swap(*x.num, *y.num);
		swap(*x.den, *y.den);
	}

	friend void compose(gf2n_rational_function & f,
			    const gf2n_rational_function & g,
			    const gf2n_polynomial & p);

	friend void shift(gf2n_rational_function & c, const gf2n_rational_function & a,
			  lidia_size_t n)   // c = a * x^n
	{
		if (n == 0)
			c.assign(a);
		else
			if (n < 0) {
				c.num->assign(*a.num);
				shift_left(*c.den, *a.den, -n);
			}
			else {
				c.den->assign(*a.den);
				shift_left(*c.num, *a.num, n);
			}
	}



	//***************************************************************
	// operators
	//***************************************************************

	friend gf2n_rational_function operator - (const gf2n_rational_function &a)
	{
		gf2n_rational_function c(a);

		return c;
	}

	friend gf2n_rational_function operator + (const gf2n_rational_function &a,
						  const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		add(c, a, b);
		return c;
	}

	friend gf2n_rational_function  operator + (const gf2n_polynomial &a,
						   const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		add(c, gf2n_rational_function(a), b);
		return c;
	}

	friend gf2n_rational_function operator + (const gf2n_rational_function &a,
						  const gf2n_polynomial &b)
	{
		gf2n_rational_function c;

		add(c, gf2n_rational_function(b), a);
		return c;
	}

	friend gf2n_rational_function operator - (const gf2n_rational_function &a,
						  const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		add(c, a, b);
		return c;
	}

	friend gf2n_rational_function  operator - (const gf2n_polynomial &a,
						   const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		add(c, gf2n_rational_function(a), b);
		return c;
	}

	friend gf2n_rational_function operator - (const gf2n_rational_function &a,
						  const gf2n_polynomial &b)
	{
		gf2n_rational_function c;

		add(c, a, gf2n_rational_function(b));
		return c;
	}

	friend gf2n_rational_function  operator * (const gf2n_rational_function &a,
						   const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		multiply(c, b, a);
		return c;
	}

	friend gf2n_rational_function operator * (const gf2n_polynomial &a,
						  const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		multiply(c, b, a);
		return c;
	}

	friend gf2n_rational_function operator * (const gf2n_rational_function &a,
						  const gf2n_polynomial &b)
	{
		gf2n_rational_function c;

		multiply(c, b, a);
		return c;
	}

	friend gf2n_rational_function  operator * (const gf2n & a,
						   const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		multiply(c, a, b);
		return c;
	}

	friend gf2n_rational_function operator / (const gf2n_rational_function &a,
						  const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		divide(c, a, b);
		return c;
	}

	friend gf2n_rational_function operator / (const gf2n_polynomial &a,
						  const gf2n_rational_function &b)
	{
		gf2n_rational_function c;

		divide(c, a, b);
		return c;
	}

	friend gf2n_rational_function operator / (const gf2n_rational_function &a,
						  const gf2n_polynomial &b)
	{
		gf2n_rational_function c;

		divide(c, a, b);
		return c;
	}


	//***************************************************************
	//
	// Modular Arithmetic without pre-conditioning
	//
	//***************************************************************

#ifndef HEADBANGER
	// arithmetic mod f.
	// all inputs and outputs are polynomials of degree less than deg(f).
	// Note that numerator and denominator are reduced modulo f, but
	// no conversion to gf2n_polynomial is done.
	// ASSUMPTION: f is assumed monic, and deg(f) > 0.
	// NOTE: if you want to do many computations with a fixed f,
	//       use the gf2n_poly_modulus data structure and associated
	//       routines below.

	friend void add_mod(gf2n_rational_function& c,
			    const gf2n_rational_function & a,
			    const gf2n_rational_function & b,
			    const gf2n_polynomial & f);

	friend void subtract_mod(gf2n_rational_function& c,
				 const gf2n_rational_function & a,
				 const gf2n_rational_function & b,
				 const gf2n_polynomial & f)
	{
		add_mod(c, a, b, f);
	}

	friend void multiply_mod(gf2n_rational_function& c,
				 const gf2n_rational_function & a,
				 const gf2n_rational_function & b,
				 const gf2n_polynomial & f)
	{
		multiply_mod(*c.num, *a.num, *b.num, f);
		multiply_mod(*c.den, *a.den, *b.den, f);
	}

	friend void square_mod(gf2n_rational_function & c,
			       const gf2n_rational_function & a,
			       const gf2n_polynomial & f)
	{
		square_mod(*c.num, *a.num, f);
		square_mod(*c.den, *a.den, f);
	}


	friend void convert_mod(gf2n_polynomial & c, const gf2n_rational_function & a,
				const gf2n_polynomial & f)
		// c = a % f, error if *a.den is not invertible
	{
		invert_mod(c, *a.den, f);
		multiply_mod(c, c, *a.num, f);
	}

	friend bool convert_mod_status(gf2n_polynomial & c,
				       const gf2n_rational_function & a,
				       const gf2n_polynomial& f)
		// if (a.den, f) = 1, returns 1 and sets c = a % f
		// otherwise, returns 0 and sets c = (a.den, f)
	{
		bool t = invert_mod_status(c, *a.den, f);
		if (t)
			multiply_mod(c, c, *a.num, f);
		return t;
	}

	friend void divide_mod(gf2n_rational_function & q,
			       const gf2n_rational_function & a,
			       const gf2n_rational_function & b,
			       const gf2n_polynomial & f)
	{
		gf2n_rational_function h(b);

		h.invert();
		multiply_mod(q, a, h, f);
	}
#endif	// HEADBANGER

	//*********************************************************************
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
	// gf2n_polynomial::read(std::istream&) can handle both formats; the decision
	// is made according to the first character being '[' or not.
	//
	//*********************************************************************

	friend std::istream& operator >> (std::istream& s, gf2n_rational_function & a)
	{
		a.read(s);
		return s;
	}

	friend std::ostream& operator << (std::ostream & s, const gf2n_rational_function & a)
	{
		a.pretty_print(s);
		return s;
	}

private:

#if 0
	void pretty_read(std::istream &s); //reads format 2)
#endif
	void pretty_print(std::ostream &s) const; //writes format 2)

public:

	void read(std::istream &s = std::cin); // reads any format
	void print(std::ostream &s = std::cout) const; //writes format 1)

	//***************************************************************
	// now the gf2n_poly_modulus functions for faster arithmetic
	//***************************************************************

	// If you need to do a lot of arithmetic modulo a fixed f,
	// build gf2n_poly_modulus F for f.  This pre-computes information about f
	// that speeds up the computation a great deal.
	// f should be monic, and deg(f) > 0.

	friend void add(gf2n_rational_function & x, const gf2n_rational_function & a,
			const gf2n_rational_function & b, gf2n_poly_modulus & F);

	friend void subtract(gf2n_rational_function & x,
			     const gf2n_rational_function & a,
			     const gf2n_rational_function & b,
			     gf2n_poly_modulus & F)
	{
		add (x, a, b, F);
	}

	friend void multiply(gf2n_rational_function & x,
			     const gf2n_rational_function & a,
			     const gf2n_rational_function & b,
			     gf2n_poly_modulus & F)
		// x = (a * b) % f, deg(a), deg(b) < f.degree()
	{
		multiply(*x.num, *a.num, *b.num, F);
		multiply(*x.den, *a.den, *b.den, F);
	}

	friend void square(gf2n_rational_function& x, const gf2n_rational_function& a,
			   gf2n_poly_modulus& F)
		// x = a^2 % f, deg(a) < f.degree()
	{
		square (*x.num, *a.num, F);
		square(*x.den, *a.den, F);
	}

	friend void divide(gf2n_rational_function& x,
			   const gf2n_rational_function& a,
			   const gf2n_rational_function & b,
			   gf2n_poly_modulus & F)
	{
		gf2n_rational_function h(b);

		h.invert();
		multiply(x, a, h, F);
	}
#endif	// HEADBANGER
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF2N_RATIONAL_FUNCTION_H_GUARD_
