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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
#define LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class qo_node;
class quadratic_order;



// Remaining includes after class definition.
// Necessary, because quadratic_number_standard .h needs
// quadratic order.h and vice versa.

//
// A quadratic number z \in Q(sqrt(Delta)) is represented as
//
//        z = (a + b sqrt(Delta))/d
//
// where a, b, d \in \Z, d > 0, and Delta is a quadratic discriminant.
// Delta must not be a square in Q.
//
// gcd(a, b, d) == 1.
//

class quadratic_number_standard
{
private:

	bigint a, b, d;
	qo_node *QO;

	static char output_mode;

	static const char WITHOUT_DISCRIMINANT;
	static const char WITH_DISCRIMINANT;

	void check_and_normalize(char *s);


public:

	//
	//  constructors / destructor
	//

	quadratic_number_standard ();
	quadratic_number_standard (const quadratic_number_standard & x);
	~quadratic_number_standard ();

	//
	//  access
	//

	bigint get_discriminant () const;
	const quadratic_order & which_order() const;
	const quadratic_order & get_order() const;

	void get (bigint & pa, bigint & pb, bigint & pd);

	const bigint & get_a () const;
	const bigint & get_b () const;
	const bigint & get_d () const;

	//
	//  assignments
	//

	quadratic_number_standard & operator = (const quadratic_number_standard & x);

	void assign_order(const quadratic_order & QO2);
	void set_order(const quadratic_order & QO2);

	void set (const bigint & pa, const bigint & pb, const bigint & pd);
	void assign (const bigint & pa, const bigint & pb, const bigint & pd);

	void assign_one ();
	void assign_one (const quadratic_order & QO2);
	void assign_zero();
	void assign_zero (const quadratic_order & QO2);

	void assign (const quadratic_number_standard & x);

	void randomize();

	friend void swap (quadratic_number_standard & x, quadratic_number_standard & y);

	//
	// comparisons
	//

	friend bool operator == (const quadratic_number_standard & x,
				 const quadratic_number_standard & y);

	friend bool operator == (const bigint & n,
				 const quadratic_number_standard & x);

	friend bool operator == (const quadratic_number_standard & x,
				 const bigint & n);

	friend bool operator == (const bigrational & n,
				 const quadratic_number_standard & x);

	friend bool operator == (const quadratic_number_standard & x,
				 const bigrational & n);

	friend bool operator != (const quadratic_number_standard & x,
				 const quadratic_number_standard & y);

	friend bool operator != (const bigint & n,
				 const quadratic_number_standard & x);

	friend bool operator != (const quadratic_number_standard & x,
				 const bigint & n);

	friend bool operator != (const bigrational & n,
				 const quadratic_number_standard & x);

	friend bool operator != (const quadratic_number_standard & x,
				 const bigrational & n);

	bool is_one() const;
	bool is_zero() const;

	bool is_positive() const;
	bool is_negative() const;

	bool is_rational_number () const;
	bool is_integer () const;
	bool is_unit () const;

	//
	//  arithmetic
	//

	void negate ();
	void invert ();

	friend void negate   (quadratic_number_standard & c,
			      const quadratic_number_standard & a);

	friend void add      (quadratic_number_standard & c,
			      const quadratic_number_standard & a,
			      const quadratic_number_standard & b);

	friend void subtract (quadratic_number_standard & c,
			      const quadratic_number_standard & a,
			      const quadratic_number_standard & b);

	friend void multiply (quadratic_number_standard & c,
			      const quadratic_number_standard & a,
			      const quadratic_number_standard & b);

	friend void multiply (quadratic_number_standard & c,
			      const bigint            & n,
			      const quadratic_number_standard & a);

	friend void multiply (quadratic_number_standard & c,
			      const quadratic_number_standard & a,
			      const bigint            & n);

	friend void multiply (quadratic_number_standard & c,
			      const bigrational      & n,
			      const quadratic_number_standard & a);

	friend void multiply (quadratic_number_standard & c,
			      const quadratic_number_standard & a,
			      const bigrational      & n);

	void multiply_by_denominator();

	friend void inverse (quadratic_number_standard & x,
			     const quadratic_number_standard & y);

	friend void divide (quadratic_number_standard & c,
			    const quadratic_number_standard & a,
			    const quadratic_number_standard & b);

	friend void divide (quadratic_number_standard         & c,
			    const bigrational      & n,
			    const quadratic_number_standard & a);

	friend void divide (quadratic_number_standard         & c,
			    const quadratic_number_standard & a,
			    const bigrational      & n);

	friend void square (quadratic_number_standard & c,
			    const quadratic_number_standard & a);

	friend void power  (quadratic_number_standard & c,
			    const quadratic_number_standard & a,
			    unsigned long e);

	friend void power  (quadratic_number_standard & c,
			    const quadratic_number_standard & a,
			    const bigint & e);

	friend quadratic_number_standard
	operator - (const quadratic_number_standard & a);

	friend
	quadratic_number_standard operator + (const quadratic_number_standard & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator - (const quadratic_number_standard & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator * (const quadratic_number_standard & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator * (const bigint & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator * (const quadratic_number_standard & a,
					      const bigint & b);
	friend
	quadratic_number_standard operator / (const quadratic_number_standard & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator / (const bigrational & a,
					      const quadratic_number_standard & b);
	friend
	quadratic_number_standard operator / (const quadratic_number_standard & a,
					      const bigrational & b);
	quadratic_number_standard & operator += (const quadratic_number_standard & x);

	quadratic_number_standard & operator -= (const quadratic_number_standard & x);

	quadratic_number_standard & operator *= (const quadratic_number_standard & x);

	quadratic_number_standard & operator /= (const quadratic_number_standard & x);

	//
	//  special purpose functions
	//

	void absolute_value();
	int get_sign() const;

	void norm (bigrational & n) const;
	bigrational norm () const;

	void trace (bigrational & t) const;
	bigrational trace () const;

	friend bigrational norm (const quadratic_number_standard & y);
	friend bigrational trace(const quadratic_number_standard & y);


	void conjugate();
	friend void conjugate (quadratic_number_standard & x, const quadratic_number_standard & y);

	long b_value () const;
	long b_value (const bigint & sq_Delta) const;

	void get_relative_approximation (xbigfloat & l, long k) const;
	void get_absolute_Ln_approximation (xbigfloat & l, long k) const;
	xbigfloat get_absolute_Ln_approximation (long k) const;

	void get_absolute_Ln_approximation_old (xbigfloat & l, long k) const;
	void get_absolute_Ln_approximation_with_check (xbigfloat & l, long k) const;
	void get_absolute_Ln_approximation_with_check_old (xbigfloat & l, long k) const;

	void get_Ln_estimate(xbigfloat & l, xbigfloat & u) const;

	bool check_Ln_correctness (lidia_size_t trials = 5) const;

	//
	//  input / output
	//

	friend std::istream & operator >> (std::istream & in, quadratic_number_standard &);

	friend std::ostream & operator << (std::ostream & out, const quadratic_number_standard &);

	friend int string_to_quadratic_number_standard (const char* s, quadratic_number_standard & x);
};

// friend functions

void swap (quadratic_number_standard & x, quadratic_number_standard & y);

    //
    // comparisons
    //

bool operator == (const quadratic_number_standard & x,
		  const quadratic_number_standard & y);
bool operator == (const bigint & n,
		  const quadratic_number_standard & x);
bool operator == (const quadratic_number_standard & x,
		  const bigint & n);
bool operator == (const bigrational & n,
		  const quadratic_number_standard & x);
bool operator == (const quadratic_number_standard & x,
		  const bigrational & n);
bool operator != (const quadratic_number_standard & x,
		  const quadratic_number_standard & y);
bool operator != (const bigint & n,
		  const quadratic_number_standard & x);
bool operator != (const quadratic_number_standard & x,
		  const bigint & n);
bool operator != (const bigrational & n,
		  const quadratic_number_standard & x);
bool operator != (const quadratic_number_standard & x,
		  const bigrational & n);

    //
    //  arithmetic
    //

void negate   (quadratic_number_standard & c,
	       const quadratic_number_standard & a);
void add      (quadratic_number_standard & c,
	       const quadratic_number_standard & a,
	       const quadratic_number_standard & b);
void subtract (quadratic_number_standard & c,
	       const quadratic_number_standard & a,
	       const quadratic_number_standard & b);
void multiply (quadratic_number_standard & c,
	       const quadratic_number_standard & a,
	       const quadratic_number_standard & b);
void multiply (quadratic_number_standard & c,
	       const bigint            & n,
	       const quadratic_number_standard & a);
void multiply (quadratic_number_standard & c,
	       const quadratic_number_standard & a,
	       const bigint            & n);
void multiply (quadratic_number_standard & c,
	       const bigrational      & n,
	       const quadratic_number_standard & a);
void multiply (quadratic_number_standard & c,
	       const quadratic_number_standard & a,
	       const bigrational      & n);
void inverse (quadratic_number_standard & x,
	      const quadratic_number_standard & y);
void divide (quadratic_number_standard & c,
	     const quadratic_number_standard & a,
	     const quadratic_number_standard & b);
void divide (quadratic_number_standard         & c,
	     const bigrational      & n,
	     const quadratic_number_standard & a);
void divide (quadratic_number_standard         & c,
	     const quadratic_number_standard & a,
	     const bigrational      & n);
void square (quadratic_number_standard & c,
	     const quadratic_number_standard & a);
void power  (quadratic_number_standard & c,
	     const quadratic_number_standard & a,
	     unsigned long e);
void power  (quadratic_number_standard & c,
	     const quadratic_number_standard & a,
	     const bigint & e);

quadratic_number_standard operator - (const quadratic_number_standard & a);
quadratic_number_standard operator + (const quadratic_number_standard & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator - (const quadratic_number_standard & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator * (const quadratic_number_standard & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator * (const bigint & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator * (const quadratic_number_standard & a,
				      const bigint & b);
quadratic_number_standard operator / (const quadratic_number_standard & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator / (const bigrational & a,
				      const quadratic_number_standard & b);
quadratic_number_standard operator / (const quadratic_number_standard & a,
				      const bigrational & b);

    //
    //  special purpose functions
    //

bigrational norm (const quadratic_number_standard & y);
bigrational trace(const quadratic_number_standard & y);

void conjugate (quadratic_number_standard & x,
		const quadratic_number_standard & y);

	//
	//  input / output
	//

std::istream & operator >> (std::istream & in, quadratic_number_standard &);
std::ostream & operator << (std::ostream & out,
			    const quadratic_number_standard &);
int string_to_quadratic_number_standard (const char* s,
					 quadratic_number_standard & x);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
