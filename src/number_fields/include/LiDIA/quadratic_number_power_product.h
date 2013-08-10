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


#ifndef LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_H_GUARD_
#define LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_RESIDUE_CLASS_LIST_H_GUARD_
# include	"LiDIA/base/residue_class_list.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_BASIS_H_GUARD_
# include	"LiDIA/quadratic_number_power_product_basis.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quadratic_order;
class quadratic_ideal;



class quadratic_number_power_product
{
private:

	typedef bigint exp_type;

	static xbigfloat                             xbigfloat_dummy;
	static quadratic_number_power_product_basis  basis_dummy;
	static residue_class_list< quadratic_number_power_product_basis > basis_list;

	residue_class< quadratic_number_power_product_basis > *basis;
	math_vector< exp_type > exp;

	static int info;
	static int default_debug_verification_value;
	int do_debug_verification;

	quadratic_number_power_product_basis* verify_preconditions(char * s) const;




#if 0
	// For evaluation modulo only.
	//
private:
	class quadratic_number_mod
	{
	private:
		bigmod a, b, d;
		static bigint Delta;

	public:
		quadratic_number_mod();
		~quadratic_number_mod();
		quadratic_number_mod & operator = (const quadratic_number_mod &);
		static void set_discriminant();
		void assign_one();
		void multiply (const quadratic_number_mod &,
			       const quadratic_number_mod &);
		void power (const quadratic_number_mod &,
			    const bigint &);
	};
#endif

	//
	//  constructor
	//
public:

	quadratic_number_power_product();

	//
	//
	//
	~quadratic_number_power_product();


	//
	//  Set the basis to \MATH{b}.
	//
	void set_basis (const quadratic_number_power_product_basis & b);

	//
	//  Set the basis to the basis of \MATH{x}.
	//
	void set_basis (const quadratic_number_power_product & x);

	//
	//  Set the exponents to \MATH{e}.
	//
	void set_exponents (const base_vector< exp_type > & e);

	//
	//  Resetting. Frees memory.
	//
	void reset();


	//
	//  Assign \MATH{x}.
	//
	quadratic_number_power_product & operator = (
		const quadratic_number_power_product & x);

	//
	//  Assign \MATH{x}.
	//
	quadratic_number_power_product & operator = (
		const quadratic_number_standard & x);

	//
	//  Assign \MATH{x}.
	//
	void assign (const quadratic_number_power_product & x);


	//
	//  Assign \MATH{g} with exponent \MATH{1}.
	//
	void assign (const quadratic_number_standard & g);

	//
	//  Assign \MATH{1} of order \MATH{O}.
	//
	void assign_one(const quadratic_order & O);

	//
	//  Let \MATH{O} be the order of \MATH{I}. It is assumed that
	//  \MATH{I} is a principal ideal in \MATH{O} with \MATH{a O = I}
	//  and \MATH{|t - Ln a| < 2^{-4}}, \MATH{sign(a) = s}.
	//  Then \MATH{a} is uniquely determined by \MATH{I, t, s}.
	//  The function \MATH{assign} initializes the power product by a
	//  compact representation of \MATH{a}. If the above conditions are
	//  not fulfilled, the result is undefined.
	//
	void assign (const quadratic_ideal & I,
		     const xbigfloat & t,
		     int s);

	//
	//  Swap the power products.
	//
	friend void swap (
		quadratic_number_power_product & x,
		quadratic_number_power_product & y);


	//
	//  Access
	//

	const quadratic_number_power_product_basis & get_basis () const;

	const base_vector< exp_type > & get_exponents () const;

	const quadratic_order & get_order() const;

	bool is_initialized() const;

private:

	bool is_one_simple(const quadratic_number_power_product_basis*) const;


public:

	int get_sign () const;

	//
	//  Arithmetic
	//

	void negate ();

	void multiply (const quadratic_number_power_product & x,
		       const quadratic_number_standard & q);

	void multiply (const quadratic_number_standard & x,
		       const quadratic_number_power_product & q);

	void multiply (const quadratic_number_power_product & x,
		       const quadratic_number_power_product & y);

	void invert (const quadratic_number_power_product & x);
	void invert_base_elements ();

	void divide (const quadratic_number_power_product & x,
		     const quadratic_number_standard & y);

	void divide (const quadratic_number_standard & x,
		     const quadratic_number_power_product & y);

	void divide (const quadratic_number_power_product & x,
		     const quadratic_number_power_product & y);

	void square (const quadratic_number_power_product & x);

	void power (const quadratic_number_power_product & x,
		    const bigint & e);

	quadratic_number_power_product & operator *= (const quadratic_number_standard & x);

	//
	//  Comparison
	//

	friend bool operator == (const quadratic_number_power_product & x,
				 const quadratic_number_power_product & y);

	friend bool operator == (const quadratic_number_power_product & x,
				 const quadratic_number_standard & q);

	friend bool operator == (const quadratic_number_standard & q,
				 const quadratic_number_power_product & x);

	friend bool operator != (const quadratic_number_power_product & x,
				 const quadratic_number_power_product & y);

	friend bool operator != (const quadratic_number_power_product & x,
				 const quadratic_number_standard & q);

	friend bool operator != (const quadratic_number_standard & q,
				 const quadratic_number_power_product & x);

	//
	//  Conjugate
	//

	void conjugate();

	//
	//   Norm
	//

	void norm_modulo(bigint & num,
			 bigint & den,
			 const bigint & m) const;


	//
	//  Evaluates the power product.
	//

	quadratic_number_standard evaluate() const;
	//	quadratic_number_standard evaluate_modulo(const bigint & m) const;


	//
	//  Let \MATH{a = *this}. The function computes an absolute \MATH{k}
	//  approximation \MATH{l} to the \MATH{Ln a}.
	//
	void get_absolute_Ln_approximation(xbigfloat & l, long k) const;
	void absolute_Ln_approximation(xbigfloat & l, long k) const;

	//
	//  Returns an absolute \MATH{k} approximation to the \MATH{Ln} of
	//  the power product.
	//
	xbigfloat get_absolute_Ln_approximation(long k) const;
	xbigfloat absolute_Ln_approximation(long k) const;

	//
	//  Let \MATH{a = *this}. The function computes a relative \MATH{k}
	//  approximation \MATH{l} to \MATH{Ln a}. It requires that \MATH{|Ln a| > 2^m}.
	//
	void get_relative_Ln_approximation(xbigfloat & l, long k, long m) const;
	void relative_Ln_approximation(xbigfloat & l, long k, long m) const;

	//
	//  Let \MATH{a = *this}. The function returns a relative \MATH{k}
	//  approximation to \MATH{Ln a}. It requires that
	//  \MATH{|Ln a| > 2^m}.
	//
	xbigfloat get_relative_Ln_approximation(long k, long m) const;
	xbigfloat relative_Ln_approximation(long k, long m) const;


	//
	//  Let \MATH{a = *this}. The function \MATH{compact\_representation}
	//  computes a compact representation of \MATH{a}. It
	//  assumes that \MATH{I} is the principal ideal generated
	//  by \MATH{a}}.
	//
	void compact_representation(const quadratic_ideal & I);

	//
	// Let \MATH{p = *this.} Let \MATH{O} be a quadratic order,
	// \MATH{rho} its fundamental unit, and let \MATH{R = |Ln(rho)|}
	// denote the regulator of \MATH{O}. If \MATH{|b-b(Ln(rho))|\LEQ
	// 1} and \MATH{p} is a non-rational number in the field of
	// fractions of \MATH{O}, then this functions computes \MATH{z},
	// such that either \MATH{-R/2 < Ln p - z R \LEQ R} or \MATH{0\LEQ
	// Ln p - z R < R} and also an absolute \MATH{k} approximation to
	// \MATH{a = Ln(p / (z \CDOT rho))}. If the conditions are not
	// fulfilled, the result is undefined.
	//
	// If \MATH{k\GEQ 4}, then \MATH{l} and the ideal \MATH{I = pO} can
	// be used to compute a compact representation of \MATH{a} using
	// the \TT{assign} function of the class. If the caller is only
	// interested in \MATH{z}, he should choose \MATH{k = -b+6}.}
	//

	void get_short_principal_ideal_generator (

		xbigfloat & l,
		bigint & z,
		const quadratic_number_power_product & rho,
		long b,
		long k) const;

	//
	//   generating unit
	//

	bool is_rational (long m) const;

	void generating_unit (
		const quadratic_number_power_product & q1,
		const quadratic_number_power_product & q2,
		long l);

	void generating_unit (
		bigint & x, bigint & y,
		bigint & M1, bigint & M2,
		const quadratic_number_power_product & q1,
		const quadratic_number_power_product & q2,
		long l);

	void generating_unit (
		base_vector< quadratic_number_power_product > & q,
		long l,
		int strategy = 0);

	void generating_unit (
		const matrix< bigint > & M,
		const base_vector< quadratic_number_standard > & q,
		long l,
		int strategy = 0);

	void generating_unit (const char * units_input_file);

	static void remove_rationals_from_quasi_units (
		base_vector< quadratic_number_power_product > & q,
		long l);


private:

	static void rgcd (bigint & x, bigint & y,
			  bigint & M1, bigint & M2,
			  const quadratic_number_power_product & q1,
			  const quadratic_number_power_product & q2,
			  long l);

	static void cfrac_expansion_bounded_den (bigint & conv_num,
						 bigint & conv_den,
						 bigint num,
						 bigint den,
						 const bigint & S);

	static long estimate_pp_accuracy (const matrix< bigint > & M,
					  long m,
					  xbigfloat c);

	//
	//  compact_representation
	//

public:

	void compact_representation_of_unit();


	//
	//  unit tests
	//

	bool could_be_unit_norm_test (int trials = 5);
	bool could_be_quasi_unit_close_test ();
	int  could_be_unit_refinement_test ();


	//
	//  Fundamental unit computation.
	//

	bool fundamental_unit_if_regulator_leq (const xbigfloat & x,
						const quadratic_order & O);

	void fundamental_unit (const quadratic_number_power_product & beta,
			       const bigint & h);


	//
	//   input / output
	//

	friend std::istream & operator >> (std::istream & in,
					   quadratic_number_power_product & b);

protected:

	void read (std::istream & in);

	friend std::ostream & operator << (std::ostream & out,
					   const quadratic_number_power_product & b);

	void write (std::ostream & out) const;

};

// friend functions

void swap (quadratic_number_power_product & x,
	   quadratic_number_power_product & y);

    //
    //  Comparison
    //

bool operator == (const quadratic_number_power_product & x,
		  const quadratic_number_power_product & y);
bool operator == (const quadratic_number_power_product & x,
		  const quadratic_number_standard & q);
bool operator == (const quadratic_number_standard & q,
		  const quadratic_number_power_product & x);
bool operator != (const quadratic_number_power_product & x,
		  const quadratic_number_power_product & y);
bool operator != (const quadratic_number_power_product & x,
		  const quadratic_number_standard & q);
bool operator != (const quadratic_number_standard & q,
		  const quadratic_number_power_product & x);

    //
    //   input / output
    //

std::istream & operator >> (std::istream & in,
			    quadratic_number_power_product & b);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_H_GUARD_
