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


#ifndef LIDIA_QUADRATIC_FORM_H_GUARD_
#define LIDIA_QUADRATIC_FORM_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_PAIR_H_GUARD_
# include	"LiDIA/pair.h"
#endif
#ifndef LIDIA_MATRIX_GL2Z_H_GUARD_
# include	"LiDIA/matrix_GL2Z.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class qi_class;
class qi_class_real;
class quadratic_ideal;
class quadratic_order;



//
// Class: quadratic_form
//
// This class represents a binary quadratic form, i.e., the polynomial in
//    x and y given by a x^2 + b xy + c y^2.  If the form is not irregular,
//    then a pointer to an element in the list of quadratic orders
//    corresponging to the quadratic order of the same discriminant as the form
//    is stored with the form.
//

class quadratic_form
{
protected:

	bigint a; // coefficient of x^2
	bigint b; // coefficient of xy
	bigint c; // coefficient of y^2
	bigint Delta; // discriminant of the form
	bigint rootD; // sqrt(Delta)
	bigint PEA_L; // nucomp and nudulp constant

	bigint norm_number();
	bigint norm_number_pos_def();
	bigint norm_number_indef();

	bool is_reduced_pos_def() const;
	bool is_reduced_indef() const;
	bool is_reduced_irregular() const;

	void almost_reduce_irregular(matrix_GL2Z & U);
	void reduce_irregular();
	void reduce_irregular(matrix_GL2Z & U);
	void reduce_pos_def();
	void reduce_pos_def(matrix_GL2Z & U);
	void reduce_indef();
	void reduce_indef(matrix_GL2Z & U);

	bool prop_equivalent_pos_def(const quadratic_form & g) const;
	bool prop_equivalent_neg_def(const quadratic_form & g) const;
	bool prop_equivalent_indef(const quadratic_form & g) const;
	bool prop_equivalent_irregular(const quadratic_form & g) const;

	bool prop_equivalent_pos_def(const quadratic_form & g, matrix_GL2Z & U) const;
	bool prop_equivalent_neg_def(const quadratic_form & g, matrix_GL2Z & U) const;
	bool prop_equivalent_indef(const quadratic_form & g, matrix_GL2Z & U) const;
	bool prop_equivalent_irregular(const quadratic_form &g, matrix_GL2Z &U) const;

	bool comp_reps_irregular(sort_vector< pair < bigint, bigint > > & Reps,
				 const bigint & N);



public:

	//
	// constructors and destructor
	//

	quadratic_form();
	quadratic_form(const bigint & a2, const bigint & b2, const bigint & c2);
	quadratic_form(const long a2, const long b2, const long c2);
	quadratic_form(const qi_class & A);
	quadratic_form(const qi_class_real & A);
	quadratic_form(const quadratic_ideal & A);
	quadratic_form(const quadratic_form & f);

	~quadratic_form();


	//
	// assignments
	//

	void assign_zero();
	void assign_one();
	void assign_one(const bigint & newDelta);
	void assign(const bigint & a2, const bigint & b2, const bigint & c2);
	void assign(const long a2, const long b2, const long c2);

	void assign(const qi_class & A);
	void assign(const qi_class_real & A);
	void assign(const quadratic_ideal & A);
	void assign(const quadratic_form & g);

	quadratic_form & operator = (const quadratic_form & g);


	//
	// access functions
	//

	bigint get_a() const;
	bigint get_b() const;
	bigint get_c() const;
	bigint discriminant() const;
	quadratic_order & which_order() const;
	friend quadratic_order & which_order(const quadratic_form & f);



	//
	// arithmetic operations
	//

	friend void compose_regular(quadratic_form & f, const quadratic_form & g1,
				    const quadratic_form & g2);
	friend void compose(quadratic_form & f, const quadratic_form & g1,
			    const quadratic_form & g2);
	friend void compose_reduce(quadratic_form & f, const quadratic_form & g1,
				   const quadratic_form & g2);
	friend void nucomp(quadratic_form & f, const quadratic_form & g1,
			   const quadratic_form & g2);
	void conjugate();
	friend void get_conjugate(quadratic_form & f, const quadratic_form & g);
	friend quadratic_form get_conjugate(const quadratic_form & f);
	friend void divide(quadratic_form & f, const quadratic_form & g1,
			   const quadratic_form & g2);
	friend void divide_reduce(quadratic_form & f, const quadratic_form & g1,
				  const quadratic_form & g2);
	friend void square(quadratic_form & f, const quadratic_form & g);
	friend void square_reduce(quadratic_form & f, const quadratic_form & g);
	friend void nudupl(quadratic_form & f, const quadratic_form & g);
	friend void power(quadratic_form & f, const quadratic_form & g,
			  const bigint & i);
	friend void power(quadratic_form & f, const quadratic_form & g, const long i);
	friend void power_reduce(quadratic_form & f, const quadratic_form & g,
				 const bigint & i);
	friend void power_reduce(quadratic_form & f, const quadratic_form & g,
				 const long i);
	friend void nupower(quadratic_form & f, const quadratic_form & g,
			    const bigint & i);
	friend void nupower(quadratic_form & f, const quadratic_form & g,
			    const long i);

	friend quadratic_form operator - (const quadratic_form & f);

	friend quadratic_form operator * (const quadratic_form & f,
					  const quadratic_form & g);
	friend quadratic_form operator / (const quadratic_form & f,
					  const quadratic_form & g);

	quadratic_form & operator *= (const quadratic_form & f);
	quadratic_form & operator /= (const quadratic_form & f);



	//
	// comparisons
	//

	bool is_zero() const;
	bool is_one() const;
	bool is_equal(const quadratic_form & g) const;
	int compare(const quadratic_form & g) const;
	int abs_compare(const quadratic_form & g) const;


	friend bool operator == (const quadratic_form & f, const quadratic_form & g);
	friend bool operator != (const quadratic_form & f, const quadratic_form & g);

	friend bool operator <= (const quadratic_form & f, const quadratic_form & g);
	friend bool operator < (const quadratic_form & f, const quadratic_form & g);
	friend bool operator >= (const quadratic_form & f, const quadratic_form & g);
	friend bool operator > (const quadratic_form & f, const quadratic_form & g);

	friend bool operator ! (const quadratic_form & f);



	//
	// basic functions
	//

	friend void swap(quadratic_form & f, quadratic_form & g);



	//
	// high level functions
	//

	int definiteness() const;
	bool is_pos_definite() const;
	bool is_pos_semidefinite() const;
	bool is_indefinite() const;
	bool is_neg_definite() const;
	bool is_neg_semidefinite() const;
	bool is_regular() const;

	bigint content() const;
	bool is_primitive() const;

	bigint eval(const bigint & x , const bigint & y) const;
	bigint operator () (const bigint & x, const bigint & y);

	void transform(const matrix_GL2Z & U);

	friend bool generate_prime_form(quadratic_form & f, const bigint & p);
	friend bool generate_prime_form(quadratic_form & f, const bigint & p,
					const bigint & newDelta);



	//
	// normalization and reduction
	//

	bool is_normal() const;
	void normalize();
	void normalize(matrix_GL2Z & U);
	void normalize_regular();

	bool is_reduced() const;
	void reduce();
	void reduce(matrix_GL2Z & U);
	void reduce_regular();

	void rho();
	void rho(matrix_GL2Z & U);
	void inverse_rho();
	void inverse_rho(matrix_GL2Z & U);

	void rho_indef();



	//
	// equivalence and principality testing
	//

	bool is_equivalent(const quadratic_form & g) const;
	bool is_prop_equivalent(const quadratic_form & g) const;
	bool is_principal() const;

	bool is_equivalent(const quadratic_form & g, matrix_GL2Z & U) const;
	bool is_prop_equivalent(const quadratic_form & g, matrix_GL2Z & U) const;
	bool is_principal(matrix_GL2Z & U) const;



	//
	// order, DL, subgroup
	//

	bigint order_in_CL() const;
	bool DL(quadratic_form & g, bigint & x) const;
	friend base_vector< bigint > subgroup(base_vector< quadratic_form > & G);
	bigfloat regulator();
	bigint class_number();
	base_vector< bigint > class_group();

	void fundamental_automorphism(matrix_GL2Z &);



	//
	// Representations
	//

#ifndef HEADBANGER
	bool representations(sort_vector< pair < bigint, bigint > > & Reps,
			     const bigint & N);
#endif



	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, quadratic_form & f);
	friend std::ostream & operator << (std::ostream & out, const quadratic_form & f);
};

// friend functions
 
quadratic_order & which_order(const quadratic_form & f);

    //
    // arithmetic operations
    //

void compose_regular(quadratic_form & f, const quadratic_form & g1,
		     const quadratic_form & g2);
void compose(quadratic_form & f, const quadratic_form & g1,
	     const quadratic_form & g2);
void compose_reduce(quadratic_form & f, const quadratic_form & g1,
		    const quadratic_form & g2);
void nucomp(quadratic_form & f, const quadratic_form & g1,
	    const quadratic_form & g2);
void get_conjugate(quadratic_form & f, const quadratic_form & g);
quadratic_form get_conjugate(const quadratic_form & f);
void divide(quadratic_form & f, const quadratic_form & g1,
	    const quadratic_form & g2);
void divide_reduce(quadratic_form & f, const quadratic_form & g1,
		   const quadratic_form & g2);
void square(quadratic_form & f, const quadratic_form & g);
void square_reduce(quadratic_form & f, const quadratic_form & g);
void nudupl(quadratic_form & f, const quadratic_form & g);
void power(quadratic_form & f, const quadratic_form & g,
	   const bigint & i);
void power(quadratic_form & f, const quadratic_form & g, const long i);
void power_reduce(quadratic_form & f, const quadratic_form & g,
		  const bigint & i);
void power_reduce(quadratic_form & f, const quadratic_form & g,
		  const long i);
void nupower(quadratic_form & f, const quadratic_form & g,
	     const bigint & i);
void nupower(quadratic_form & f, const quadratic_form & g,
	     const long i);

quadratic_form operator - (const quadratic_form & f);
quadratic_form operator * (const quadratic_form & f,
			   const quadratic_form & g);
quadratic_form operator / (const quadratic_form & f,
			   const quadratic_form & g);

    //
    // comparisons
    //

bool operator == (const quadratic_form & f, const quadratic_form & g);
bool operator != (const quadratic_form & f, const quadratic_form & g);
bool operator <= (const quadratic_form & f, const quadratic_form & g);
bool operator < (const quadratic_form & f, const quadratic_form & g);
bool operator >= (const quadratic_form & f, const quadratic_form & g);
bool operator > (const quadratic_form & f, const quadratic_form & g);
bool operator ! (const quadratic_form & f);


    //
    // basic functions
    //

void swap(quadratic_form & f, quadratic_form & g);


    //
    // high level functions
    //

bool generate_prime_form(quadratic_form & f, const bigint & p);
bool generate_prime_form(quadratic_form & f, const bigint & p,
			 const bigint & newDelta);

    //
    // order, DL, subgroup
    //

base_vector< bigint > subgroup(base_vector< quadratic_form > & G);

    //
    // input/output
    //

std::istream & operator >> (std::istream & in, quadratic_form & f);
std::ostream & operator << (std::ostream & out, const quadratic_form & f);



// key function for hash tables
bigint quadratic_form_key(const quadratic_form & G);


// Computes the reduced representatives of the class group
base_vector< quadratic_form >
compute_class_group(const bigint & Delta);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_FORM_H_GUARD_
