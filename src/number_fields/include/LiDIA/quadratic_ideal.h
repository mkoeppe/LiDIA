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
//	Author	: Michael Jacobson (MJJ), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QUADRATIC_IDEAL_H_GUARD_
#define LIDIA_QUADRATIC_IDEAL_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
#include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
#include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
#include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
#include	"LiDIA/xbigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
#include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class qo_node;
class quadratic_order;
class quadratic_form;
class qi_class;
class qi_class_real;

class quadratic_number_standard;
class quadratic_number_power_product;



//
// Class: quadratic_ideal
//
// This class represents a normalized fractional ideal of a quadratic order,
// i.e., the module q [ aZ + (b+sqrt(Delta))/2 Z ] where a and b are
// integers and q is rational. The quadratic order to which an instance of
// this class belongs is stored in a list of current quadratic orders, and
// a pointer to its list element is stored with each instance of this
// class. The form (a, b, c) is always normalized.
//

class quadratic_ideal
{
protected:

	bigint a; // coefficient in the representation of the ideal
	bigint b; // coefficient in the representation of the ideal
	bigrational q; // coefficient in the representation of the ideal
	qo_node *QO; // pointer to a quadratic order list entry

	bool reduced_flag; // true, implies that ideal is reduced

	void normalize();
	void normalize_imag();
	void normalize_real();

	bool is_normalized_imag() const;
	bool is_normalized_real() const;

	// < MM >
	void reduce_imag(quadratic_number_standard & g, bool comp_red_number = false);
	void reduce_real(quadratic_number_standard & g, bool comp_red_number = false);

	void rho_imag(quadratic_number_standard & g);
	void rho_real(quadratic_number_standard & g);

	void inverse_rho_imag(quadratic_number_standard & g);
	void inverse_rho_real(quadratic_number_standard & g);
	// < /MM >


	friend void multiply_imag(quadratic_ideal & C, const quadratic_ideal & A,
				  const quadratic_ideal & B);
	friend void multiply_real(quadratic_ideal & C, const quadratic_ideal & A,
				  const quadratic_ideal & B);

	friend void square_imag(quadratic_ideal & C, const quadratic_ideal & A);
	friend void square_real(quadratic_ideal & C, const quadratic_ideal & A);



public:

	//
	// constructors and destructor
	//

	quadratic_ideal();
	quadratic_ideal(const bigint & a2, const bigint & b2,
			const bigrational & q2);
	quadratic_ideal(const long a2, const long b2,
			const bigrational & q2);
	quadratic_ideal(const bigint & a2, const bigint & b2,
			const bigrational & q2, const quadratic_order & QO2);
	quadratic_ideal(const long a2, const long b2,
			const bigrational & q2, const quadratic_order & QO2);
	quadratic_ideal(const quadratic_form & qf);
	quadratic_ideal(const qi_class & A);
	quadratic_ideal(const qi_class_real & A);
	quadratic_ideal(const quadratic_ideal & A);

	// < MM >
	quadratic_ideal(const quadratic_order & QO2);
	// < /MM >

	~quadratic_ideal();

	//
	// assignments
	//

	void assign_zero();
	void assign_one();
	void assign_principal(const bigint & x, const bigint & y);
	bool assign(const bigint & a2, const bigint & b2, const bigrational & q2);
	bool assign(const long a2, const long b2, const bigrational & q2);

	void assign_zero(const quadratic_order & QO2);
	void assign(const quadratic_order & QO2);
	void assign_one(const quadratic_order & QO2);
	void assign_principal(const bigint & x, const bigint & y,
			      const quadratic_order & QO2);
	bool assign(const bigint & a2, const bigint & b2,
		    const bigrational & q2, const quadratic_order & QO2);
	bool assign(const long a2, const long b2,
		    const bigrational & q2, const quadratic_order & QO2);

	void assign(const quadratic_form & qf);
	void assign(const qi_class & B);
	void assign(const qi_class_real & B);
	void assign(const quadratic_ideal & B);

	bool assign_order(const quadratic_order & QO2);
	// < MM >
	bool set_order(const quadratic_order & QO2);
	// < /MM >

	quadratic_ideal & operator = (const quadratic_ideal & A);

	// < MM >
private:
	void assign_module_of (const bigrational & nq,
			       quadratic_number_standard q1,
			       quadratic_number_standard q2);

public:
	bool assign (const bigrational & nq,
		     const quadratic_number_standard & q1,
		     const quadratic_number_standard & q2);
	// < /MM >


	//
	// access functions
	//

	bigint get_a() const;
	bigint get_b() const;
	bigrational get_q() const;
	bigint get_c() const;
	const quadratic_order & which_order() const;
	friend const quadratic_order & which_order(const quadratic_ideal & A);
	bigint discriminant() const;

	// < MM >
	const quadratic_order & get_order() const;
	// < /MM >



	//
	// arithmetic operations
	//

#ifndef HEADBANGER
	friend void add(quadratic_ideal & C, const quadratic_ideal & A,
			const quadratic_ideal & B);
	friend void intersect(quadratic_ideal & C, const quadratic_ideal & A,
			      const quadratic_ideal & B);
#endif
	friend void multiply(quadratic_ideal & C, const quadratic_ideal & A,
			     const quadratic_ideal & B);
	friend void divide(quadratic_ideal & C, const quadratic_ideal & A,
			   const quadratic_ideal & B);
	void invert();
	friend void inverse(quadratic_ideal & A, const quadratic_ideal & B);
	friend quadratic_ideal inverse(const quadratic_ideal & B);
	void conjugate();
	friend void get_conjugate(quadratic_ideal & A, const quadratic_ideal & B);
	friend quadratic_ideal get_conjugate(const quadratic_ideal & B);
	friend void square(quadratic_ideal & C, const quadratic_ideal & A);
	friend void power(quadratic_ideal & C, const quadratic_ideal & A,
			  const bigint & i);
	friend void power(quadratic_ideal & C, const quadratic_ideal & A,
			  const long i);

	friend quadratic_ideal operator - (const quadratic_ideal & A);

#ifndef HEADBANGER
	friend quadratic_ideal operator + (const quadratic_ideal & A,
					   const quadratic_ideal & B);
	friend quadratic_ideal operator & (const quadratic_ideal & A,
					   const quadratic_ideal & B);
#endif
	friend quadratic_ideal operator * (const quadratic_ideal & A,
					   const quadratic_ideal & B);
	friend quadratic_ideal operator / (const quadratic_ideal & A,
					   const quadratic_ideal & B);

#ifndef HEADBANGER
	quadratic_ideal & operator += (const quadratic_ideal & A);
	quadratic_ideal & operator &= (const quadratic_ideal & A);
#endif
	quadratic_ideal & operator *= (const quadratic_ideal & A);
	quadratic_ideal & operator /= (const quadratic_ideal & A);


	// < MM >
	friend void multiply (quadratic_ideal & J,
			      const quadratic_ideal & I,
			      const quadratic_number_standard & g);

	friend void divide (quadratic_ideal & J,
			    const quadratic_ideal & I,
			    const quadratic_number_standard & g);
	// < /MM >

	//
	// comparisons
	//

	bool is_zero() const;
	bool is_one() const;
	bool is_equal(const quadratic_ideal & B) const;
	bool is_subset(const quadratic_ideal & B) const;
	bool is_proper_subset(const quadratic_ideal & B) const;

	friend bool operator == (const quadratic_ideal &A, const quadratic_ideal &B);
	friend bool operator != (const quadratic_ideal &A, const quadratic_ideal &B);

	friend bool operator <= (const quadratic_ideal &A, const quadratic_ideal &B);
	friend bool operator < (const quadratic_ideal &A, const quadratic_ideal &B);
	friend bool operator >= (const quadratic_ideal &A, const quadratic_ideal &B);
	friend bool operator > (const quadratic_ideal &A, const quadratic_ideal &B);

	friend bool operator ! (const quadratic_ideal & A);



	//
	// basic functions
	//

	friend void swap(quadratic_ideal & A, quadratic_ideal & B);


	//
	// high level functions
	//

	bigrational smallest_rational() const;
	bool is_integral() const;
	bigint conductor() const;
	bool is_invertible() const;
	bigrational norm() const;
	void ring_of_multipliers(quadratic_order & rm);
	friend bool generate_prime_ideal(quadratic_ideal & A, const bigint & p);
	friend bool generate_prime_ideal(quadratic_ideal & A, const bigint & p,
					 const quadratic_order & QO2);
	// < MM >
	bigint denominator() const;
	void multiply_by_denominator();
	// < /MM >

	//
	// reduction
	//

	bool is_reduced() const;
	void reduce();
	void reduce(quadratic_number_standard & g);

	void rho();
	void rho(quadratic_number_standard & g);
	friend void apply_rho(quadratic_ideal & A, const quadratic_ideal & B);
	friend quadratic_ideal apply_rho(const quadratic_ideal & B);

	void inverse_rho();
	void inverse_rho(quadratic_number_standard & g);
	friend void apply_inverse_rho(quadratic_ideal & A, const quadratic_ideal & B);
	friend quadratic_ideal apply_inverse_rho(const quadratic_ideal & B);

	// for debugging purposes; should always return true;
	bool is_normalized () const;


	// < MM >
	void local_close(quadratic_number_standard & alpha,
			 xbigfloat & a,
			 xbigfloat t,
			 long k);

	void order_close(quadratic_number_power_product & alpha,
			 xbigfloat & a,
			 xbigfloat t,
			 long k);

	void close      (quadratic_number_power_product & alpha,
			 xbigfloat & a,
			 xbigfloat t,
			 long k);
	// < /MM >


	//
	// equivalence and principality testing
	//

	bool is_equivalent(const quadratic_ideal & B) const;
	bool is_principal() const;



	//
	// order, DL, subgroup
	//

	bigint order_in_CL() const;
	bool DL(quadratic_ideal & G, bigint & x) const;
	friend base_vector< bigint > subgroup(base_vector< quadratic_ideal > & G);
	bigfloat regulator();
	bigint class_number();
	base_vector< bigint > class_group();


	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, quadratic_ideal & A);
	friend std::ostream & operator << (std::ostream & out, const quadratic_ideal & A);
};

// friend functions

void multiply_imag(quadratic_ideal & C, const quadratic_ideal & A,
		   const quadratic_ideal & B);
void multiply_real(quadratic_ideal & C, const quadratic_ideal & A,
		   const quadratic_ideal & B);

void square_imag(quadratic_ideal & C, const quadratic_ideal & A);
void square_real(quadratic_ideal & C, const quadratic_ideal & A);

    //
    // access functions
    //

const quadratic_order & which_order(const quadratic_ideal & A);

    //
    // arithmetic operations
    //

#ifndef HEADBANGER
void add(quadratic_ideal & C, const quadratic_ideal & A,
	 const quadratic_ideal & B);
void intersect(quadratic_ideal & C, const quadratic_ideal & A,
	       const quadratic_ideal & B);
#endif
void multiply(quadratic_ideal & C, const quadratic_ideal & A,
	      const quadratic_ideal & B);
void divide(quadratic_ideal & C, const quadratic_ideal & A,
	    const quadratic_ideal & B);
void inverse(quadratic_ideal & A, const quadratic_ideal & B);
quadratic_ideal inverse(const quadratic_ideal & B);
void get_conjugate(quadratic_ideal & A, const quadratic_ideal & B);
quadratic_ideal get_conjugate(const quadratic_ideal & B);
void square(quadratic_ideal & C, const quadratic_ideal & A);
void power(quadratic_ideal & C, const quadratic_ideal & A,
	   const bigint & i);
void power(quadratic_ideal & C, const quadratic_ideal & A,
	   const long i);

quadratic_ideal operator - (const quadratic_ideal & A);

#ifndef HEADBANGER
quadratic_ideal operator + (const quadratic_ideal & A,
			    const quadratic_ideal & B);
quadratic_ideal operator & (const quadratic_ideal & A,
			    const quadratic_ideal & B);
#endif
quadratic_ideal operator * (const quadratic_ideal & A,
			    const quadratic_ideal & B);
quadratic_ideal operator / (const quadratic_ideal & A,
			    const quadratic_ideal & B);

    // < MM >
void multiply (quadratic_ideal & J,
	       const quadratic_ideal & I,
	       const quadratic_number_standard & g);
void divide (quadratic_ideal & J,
	     const quadratic_ideal & I,
	     const quadratic_number_standard & g);
    // < /MM >

    //
    // comparisons
    //

bool operator == (const quadratic_ideal &A, const quadratic_ideal &B);
bool operator != (const quadratic_ideal &A, const quadratic_ideal &B);

bool operator <= (const quadratic_ideal &A, const quadratic_ideal &B);
bool operator < (const quadratic_ideal &A, const quadratic_ideal &B);
bool operator >= (const quadratic_ideal &A, const quadratic_ideal &B);
bool operator > (const quadratic_ideal &A, const quadratic_ideal &B);

bool operator ! (const quadratic_ideal & A);

    //
    // basic functions
    //

void swap(quadratic_ideal & A, quadratic_ideal & B);

    //
    // high level functions
    //

bool generate_prime_ideal(quadratic_ideal & A, const bigint & p);
bool generate_prime_ideal(quadratic_ideal & A, const bigint & p,
			  const quadratic_order & QO2);

    //
    // reduction
    //

void apply_rho(quadratic_ideal & A, const quadratic_ideal & B);
quadratic_ideal apply_rho(const quadratic_ideal & B);
void apply_inverse_rho(quadratic_ideal & A, const quadratic_ideal & B);
quadratic_ideal apply_inverse_rho(const quadratic_ideal & B);

    //
    // order, DL, subgroup
    //

base_vector< bigint > subgroup(base_vector< quadratic_ideal > & G);

    //
    // input/output
    //

std::istream & operator >> (std::istream & in, quadratic_ideal & A);
std::ostream & operator << (std::ostream & out, const quadratic_ideal & A);


// key function for hash tables
bigint quadratic_ideal_key(const quadratic_ideal & G);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_IDEAL_H_GUARD_
