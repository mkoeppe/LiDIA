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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QI_CLASS_REAL_H_GUARD_
#define LIDIA_QI_CLASS_REAL_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_QI_CLASS_H_GUARD_
# include	"LiDIA/qi_class.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quadratic_form;
class quadratic_ideal;



//
// Class: qi_class_real
//
// This class represents an element in the class group of a real quadratic
//    order, i.e., a reduced, primitive, invertible ideal
//    aZ + (b+sqrt(Delta))/2 Z  where Delta is the discriminant of the real
//    quadratic order.  This class is derived from qi_class and shares all
//    its properties.  In addition, a distance is associated with each
//    instance of this class, so that computations in the infrastructure are
//    possible.
//


class qi_class_real : public qi_class
{
protected:

	xbigfloat d; // distance associated with the representative

	static int info; // whether to print timings or not
	static int do_verify; // whether to verify computations


	void normalize_real();
	void reduce();


public:

	//
	// constructors and destructor
	//

	qi_class_real();
	qi_class_real(const bigint & a2, const bigint & b2, const xbigfloat & dist);
	qi_class_real(const long a2, const long b2, const xbigfloat & dist);
	qi_class_real(const bigint & a2, const bigint & b2);
	qi_class_real(const long a2, const long b2);
	qi_class_real(const quadratic_form & qf);
	qi_class_real(const quadratic_ideal & A);
	qi_class_real(const qi_class & A);
	qi_class_real(const qi_class_real & A);

	~qi_class_real();



	//
	// initialization
	//

	static void verbose(int state);
	static void verification(int level);



	//
	// assignment
	//

	void assign_zero();
	void assign_one();
	void assign_principal(const bigint & x, const bigint & y,
			      const xbigfloat & dist);
	bool assign(const bigint & a2, const bigint & b2, const xbigfloat & dist);
	bool assign(const long a2, const long b2, const xbigfloat & dist);
	void assign(const quadratic_form & qf, const xbigfloat & dist);
	void assign(const quadratic_ideal & B, const xbigfloat & dist);
	void assign(const qi_class & B, const xbigfloat & dist);

	void assign_principal(const bigint & x, const bigint & y);
	bool assign(const bigint & a2, const bigint & b2);
	bool assign(const long a2, const long b2);
	void assign(const quadratic_form & qf);
	void assign(const quadratic_ideal & B);
	void assign(const qi_class & B);

	void assign(const qi_class_real & B);

	qi_class_real & operator = (const qi_class_real & A);



	//
	// access functions
	//

	bigfloat get_distance() const;
	xbigfloat get_distance_x() const;



	//
	// arithmetic operations
	//


	friend void multiply(qi_class_real & C, const qi_class_real & A,
			     const qi_class_real & B);
	friend void multiply_real(qi_class_real & C, const qi_class_real & A,
				  const qi_class_real & B);

	void invert();
	friend void inverse(qi_class_real & A, const qi_class_real & B);
	friend qi_class_real inverse(const qi_class_real & B);

	friend void divide(qi_class_real & C, const qi_class_real & A,
			   const qi_class_real & B);

	friend void square(qi_class_real & C, const qi_class_real & A);
	friend void square_real(qi_class_real & C, const qi_class_real & A);

	friend void power(qi_class_real & C, const qi_class_real & A,
			  const bigint & i);
	friend void power(qi_class_real & C, const qi_class_real & A, const long i);
	friend void power_real(qi_class_real & C, const qi_class_real & A,
			       const bigint & i);
	friend void power_real(qi_class_real & C, const qi_class_real & A,
			       const long i);

	friend qi_class_real nearest(const qi_class_real & S, const bigfloat & E);

	friend qi_class_real operator - (const qi_class_real & A);

	friend qi_class_real operator * (const qi_class_real & A,
					 const qi_class_real & B);
	friend qi_class_real operator / (const qi_class_real & A,
					 const qi_class_real & B);

	qi_class_real & operator *= (const qi_class_real & A);
	qi_class_real & operator /= (const qi_class_real & A);



	//
	// basic functions
	//

	friend void swap(qi_class_real & A, qi_class_real & B);



	//
	// high_level functions
	//

	friend bool generate_prime_ideal(qi_class_real & A, const bigint & p);



	//
	// reduction operators
	//

	void rho();
	friend void apply_rho(qi_class_real & A, const qi_class_real & B);
	friend qi_class_real apply_rho(const qi_class_real & B);

	void inverse_rho();
	friend void apply_inverse_rho(qi_class_real & A, const qi_class_real & B);
	friend qi_class_real apply_inverse_rho(const qi_class_real & B);

	bigfloat convert_distance(const qi_class_real & A,
				  const bigfloat & dist) const;

#ifndef HEADBANGER
	void adjust_pos(const bigfloat & dist);
	void adjust_neg(const bigfloat & dist);
	void adjust_abs(const bigfloat & dist);

	bigfloat make_list(bigint & oa, const bigfloat & dist,
			   hash_table< qi_class_real > & HT);
	int make_list(bigint & oa, const bigfloat & dist, const qi_class_real & G,
		      const qi_class_real & U, hash_table< qi_class_real > & HT);
	void make_list(const bigfloat & dist, long x, hash_table< ideal_node > & HT);
	void make_list(const bigfloat & dist, long x,
		       indexed_hash_table< ideal_node > & HT);
#endif



	//
	// equivalence and principality testing
	//

	bool is_principal() const;
	bool is_principal_buch() const;
	bool is_principal_subexp() const;

	bool is_principal(bigfloat & dist) const;
	bool is_principal_buch(bigfloat & dist) const;
	bool is_principal_subexp(bigfloat & dist) const;

	bool verify_principal(const bigfloat & x, const bool is_prin) const;


	bool is_equivalent(const qi_class_real & B) const;
	bool is_equivalent_buch(const qi_class_real & B) const;
	bool is_equivalent_subexp(const qi_class_real & B) const;

	bool is_equivalent(const qi_class_real & B, bigfloat & dist) const;
	bool is_equivalent_buch(const qi_class_real & B, bigfloat & dist) const;
	bool is_equivalent_subexp(const qi_class_real & B, bigfloat & dist) const;

	bool verify_equivalent(const qi_class_real & B, const bigfloat & x,
			       const bool is_equiv) const;

	bool DL(const qi_class_real & G, bigint & x, bigfloat & dist) const;
	bool verify_DL(const qi_class_real & G, const bigint & x,
		       const bigfloat & dist, const bool is_DL) const;

	friend base_vector< bigint > subgroup(base_vector< qi_class_real > & G);

	// subgroup computation with no knowledge of class number
	friend base_vector< bigint > subgroup_BJT(base_vector< qi_class_real > & G,
						  base_vector< long > & v);

#ifndef HEADBANGER
	// subgroup computation when class number is known
	friend base_vector< bigint > subgroup_h(base_vector< qi_class_real > & G);
#endif


	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, qi_class_real & A);
	friend std::ostream & operator << (std::ostream & out, const qi_class_real & A);
};



// key function for hash tables
bigint qi_class_real_key(const qi_class_real & G);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifndef LIDIA_QUADRATIC_ORDER_H_GUARD_
# include	"LiDIA/quadratic_order.h"
#endif
#ifndef LIDIA_QUADRATIC_IDEAL_H_GUARD_
# include	"LiDIA/quadratic_ideal.h"
#endif
#ifndef LIDIA_QUADRATIC_FORM_H_GUARD_
# include	"LiDIA/quadratic_form.h"
#endif



#endif	// LIDIA_QI_CLASS_REAL_H_GUARD_
