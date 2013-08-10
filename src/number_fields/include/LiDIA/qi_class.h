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


#ifndef LIDIA_QI_CLASS_H_GUARD_
#define LIDIA_QI_CLASS_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class qo_node;
class qi_class_real;
class quadratic_ideal;
class quadratic_form;
class quadratic_order;
class ideal_node;
class quadratic_number_standard;
template< class T > class indexed_hash_table;
template< class T > class hash_table;
class rational_factorization;



//
// Class: qi_class
//
// This class represents an element in the class group of a quadratic order,
//    i.e., a reduced, primitive, invertible ideal aZ + (b+sqrt(Delta))/2 Z
//    where Delta is the discriminant of the quadratic order.  If the
//    quadratic order is imaginary, then this representation is unique.
//    Otherwise, the element is represented by a finite set of reduced
//    representatives, which can be cycled through by applying the reduction
//    operators rho and rho^-1.
// All instances of this class are assumed to belong to the same quadratic
//    order. This quadratic order is stored in the value current_order, which
//    is a pointer to an element in the list of current quadratic orders.
// A qi_class is represented by the bigint coefficients a and b.  The
//    following information is stored for use in some computations:
//
//     Delta = discriminant of the current quadratic order
//     rootD = floor(sqrt(|Delta))
//     rd = sqrt(Delta)
//     PEA_L = floor((|Delta|/4)^1/4) (used for nucomp and nudupl)
//     ordmult = largest multiple of an element's order found by order_shanks
//     omfact = factorization of ordmult
//

#define OBJT_RB 9		// for order, use BJT if D < 10^12
#define OSHANKS_RB 10		// for order, use Shanks if D < 10^17
#define DLBJT_RB 9		// for DL, use BJT if D < 10^0

#define OBJT_IB 12		// for order, use BJT if |D| < 10^12
#define OSHANKS_IB 12		// for order, use Shanks if |D| < 10^12
#define DLBJT_IB 13		// for DL, use BJT if |D| < 10^13

class qi_class
{
protected:

	bigint a; // coefficient in ideal representation
	bigint b; // coefficient in ideal representation

	static bigint Delta; // discriminant of the current order
	static bigint rootD; // floor(sqrt(|Delta|))
	static xbigfloat rd; // sqrt(|Delta|)
	static long xprec; // xbigfloat precision to be used

	static qo_node *current_order; // current quadratic order list entry

	static bigint ordmult; // largest exponent from order_shanks
	static rational_factorization omfact; // factorization of ordmult

	static int info; // whether to print timings or not
	static int do_verify; // whether to verify computations

	bigint conductor() const;
	bool is_invertible() const;

	void normalize();
	void normalize_imag();
	void normalize_real();

	void reduce();
	void reduce_imag();
	void reduce_real();
	void reduce_real(quadratic_number_standard & alpha);

	// order of element using BJT algorithm
	bigint oBJT_imag(long v) const;
	bigint oBJT_imag(long v, bigint & GO, bigint & NS) const;
	bigint oBJT_real() const;

	// order of element given multiple of its order
	bigint omult_imag(const bigint & h,
			  const rational_factorization & hfact) const;
	bigint omult_real(const bigint & h,
			  const rational_factorization & hfact) const;
	bigint omult_real_sub(const bigint & h,
			      const rational_factorization & hfact) const;

	// order of element using Shanks algorithm
	bigint os_imag(long v, const bigint & hstar) const;
	bigint os_real(long v, const bigint & hstar) const;

	// DL using BJT algorithm
	bool DLBJT_imag(const qi_class & G, bigint & x, long v) const;
	bool DLBJT_real(const qi_class & G, bigint & x) const;

	// DL when class number is known
	bool DLh_imag(const qi_class & G, bigint & x) const;
	bool DLh_real(const qi_class & G, bigint & x) const;

	// DL using subexponential algorithm
	bool DLsubexp(const qi_class & G, bigint & x) const;

	// structure of subgroup using BJT algorithm
	friend base_vector< bigint > subBJT_imag(base_vector< qi_class > & G,
						 base_vector< long > & v);
	friend base_vector< bigint > subBJT_real(base_vector< qi_class > & G);

	// structure of subgroup when class number is known
	friend base_vector< bigint > subh_imag(base_vector< qi_class > & G);
	friend base_vector< bigint > subh_real(base_vector< qi_class > & G);


public:

	static bigint PEA_L; // floor((Delta/4)^1/4)

	//
	// constructors and destructor
	//

	qi_class();
	qi_class(const bigint & a2, const bigint & b2);
	qi_class(const long a2, const long b2);
	qi_class(const quadratic_form & qf);
	qi_class(const quadratic_ideal & A);
	qi_class(const qi_class_real & A);
	qi_class(const qi_class & A);

	~qi_class();



	//
	// initialization
	//

	static void set_current_order(quadratic_order & QO);
	static void verbose(int state);
	static void verification(int level);



	//
	// assignment
	//

	void assign_zero();
	void assign_one();
	void assign_principal(const bigint & x, const bigint & y);
	bool assign(const bigint & a2, const bigint & b2);
	bool assign(const long a2, const long b2);

	void assign(const quadratic_form & qf);
	void assign(const quadratic_ideal & B);
	void assign(const qi_class_real & B);
	void assign(const qi_class & B);

	qi_class & operator = (const qi_class & A);



	//
	// access functions
	//

	static quadratic_order * get_current_order_ptr();
	static quadratic_order & get_current_order();
	static bigint discriminant();
	static bigint get_rootD();
	static xbigfloat get_rd();

	bigint get_a() const;
	bigint get_b() const;
	bigint get_c() const;



	//
	// arithmetic operations
	//

	friend void multiply(qi_class & C, const qi_class & A, const qi_class & B);
	friend void multiply_imag(qi_class &C, const qi_class &A, const qi_class &B);
	friend void nucomp(qi_class & C, const qi_class & A, const qi_class & B);
	friend void multiply_real(qi_class &C, const qi_class &A, const qi_class &B);
	friend void multiply_real(qi_class &C, const qi_class &A, const qi_class &B,
				  quadratic_number_standard &Cq, const quadratic_number_standard &Aq,
				  const quadratic_number_standard &Bq);

	void invert();
	friend void inverse(qi_class & A, const qi_class & B);
	friend qi_class inverse(const qi_class & B);

	friend void divide(qi_class & C, const qi_class & A, const qi_class & B);

	friend void square(qi_class & C, const qi_class & A);
	friend void square_imag(qi_class & C, const qi_class & A);
	friend void nudupl(qi_class & C, const qi_class & A);
	friend void square_real(qi_class & C, const qi_class & A);

	friend void power(qi_class & C, const qi_class & A, const bigint & i);
	friend void power(qi_class & C, const qi_class & A, const long i);
	friend void power_imag(qi_class & C, const qi_class & A, const bigint & i);
	friend void power_imag(qi_class & C, const qi_class & A, const long i);
	friend void nupower(qi_class & C, const qi_class & A, const bigint & i);
	friend void nupower(qi_class & C, const qi_class & A, const long i);
	friend void power_real(qi_class & C, const qi_class & A, const bigint & i);
	friend void power_real(qi_class & C, const qi_class & A, const long i);

	friend qi_class operator - (const qi_class & A);

	friend qi_class operator * (const qi_class & A, const qi_class & B);
	friend qi_class operator / (const qi_class & A, const qi_class & B);

	qi_class & operator *= (const qi_class & A);
	qi_class & operator /= (const qi_class & A);



	//
	// comparisons
	//

	bool is_zero() const;
	bool is_one() const;
	bool is_equal(const qi_class & B) const;

	friend bool operator == (const qi_class & A, const qi_class & B);
	friend bool operator != (const qi_class & A, const qi_class & B);

	friend bool operator ! (const qi_class & A);



	//
	// basic functions
	//

	friend void swap(qi_class & A, qi_class & B);



	//
	// high level functions
	//

	friend bool generate_prime_ideal(qi_class & A, const bigint & p);



	//
	// reduction operators
	//

	void rho();
	void rho(bigint & oa, quadratic_number_standard & alpha);
	friend void apply_rho(qi_class & A, const qi_class & B);
	friend qi_class apply_rho(const qi_class & B);

	void inverse_rho();
	friend void apply_inverse_rho(qi_class & A, const qi_class & B);
	friend qi_class apply_inverse_rho(const qi_class & B);

#ifndef HEADBANGER
	void make_cycle(long x, hash_table< ideal_node > & HT);
	void make_cycle(long x, indexed_hash_table< ideal_node > & HT);
#endif



	//
	// equivalence and principality testing
	//

	bool is_equivalent(const qi_class & B) const;
	bool is_principal() const;



	//
	// order of element in class group
	//

	bigint order_in_CL() const;

	// order computation with no knowledge of class number
	bigint order_BJT(long v) const;

	// order compution when class number is known
	bigint order_h() const;

	// order compution when multiple of order is known
	bigint order_mult(const bigint & h,
			  const rational_factorization & hfact) const;

	// order computation using estimate of L(1), class number not known
	bigint order_shanks() const;

	// subexponential order computation
	bigint order_subexp() const;

	// verification
	bool verify_order(const bigint & x) const;


	//
	// DL computation
	//

	bool DL(const qi_class & G, bigint & x) const;

	// DL computation with no knowledge of class number
	bool DL_BJT(const qi_class & G, bigint & x) const;
	bool DL_BJT(const qi_class & G, bigint & x, long v) const;

#ifndef HEADBANGER
	// DL computation when class number is known
	bool DL_h(const qi_class & G, bigint & x) const;
#endif

	// subexponential DL computation
	bool DL_subexp(const qi_class & G, bigint & x) const;

	// verification
	bool verify_DL(const qi_class & G, const bigint & x, const bool is_DL) const;


	//
	// subgroup computation
	//

	friend base_vector< bigint > subgroup(base_vector< qi_class > & G);

	// subgroup computation with no knowledge of class number
	friend base_vector< bigint > subgroup_BJT(base_vector< qi_class > & G,
						  base_vector< long > & v);

#ifndef HEADBANGER
	// subgroup computation when class number is known
	friend base_vector< bigint > subgroup_h(base_vector< qi_class > & G);
#endif



	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, qi_class & A);
	friend std::ostream & operator << (std::ostream & out, const qi_class & A);
};


// friend functions

base_vector< bigint > subBJT_imag(base_vector< qi_class > & G,
				  base_vector< long > & v);
base_vector< bigint > subBJT_real(base_vector< qi_class > & G);

    // structure of subgroup when class number is known
base_vector< bigint > subh_imag(base_vector< qi_class > & G);
base_vector< bigint > subh_real(base_vector< qi_class > & G);

    //
    // arithmetic operations
    //

void multiply(qi_class & C, const qi_class & A, const qi_class & B);
void multiply_imag(qi_class &C, const qi_class &A, const qi_class &B);
void nucomp(qi_class & C, const qi_class & A, const qi_class & B);
void multiply_real(qi_class &C, const qi_class &A, const qi_class &B);
void multiply_real(qi_class &C, const qi_class &A, const qi_class &B,
		   quadratic_number_standard &Cq,
		   const quadratic_number_standard &Aq,
		   const quadratic_number_standard &Bq);

void inverse(qi_class & A, const qi_class & B);
qi_class inverse(const qi_class & B);

void divide(qi_class & C, const qi_class & A, const qi_class & B);

void square(qi_class & C, const qi_class & A);
void square_imag(qi_class & C, const qi_class & A);
void nudupl(qi_class & C, const qi_class & A);
void square_real(qi_class & C, const qi_class & A);

void power(qi_class & C, const qi_class & A, const bigint & i);
void power(qi_class & C, const qi_class & A, const long i);
void power_imag(qi_class & C, const qi_class & A, const bigint & i);
void power_imag(qi_class & C, const qi_class & A, const long i);
void nupower(qi_class & C, const qi_class & A, const bigint & i);
void nupower(qi_class & C, const qi_class & A, const long i);
void power_real(qi_class & C, const qi_class & A, const bigint & i);
void power_real(qi_class & C, const qi_class & A, const long i);

qi_class operator - (const qi_class & A);
qi_class operator * (const qi_class & A, const qi_class & B);
qi_class operator / (const qi_class & A, const qi_class & B);

    //
    // comparisons
    //

bool operator == (const qi_class & A, const qi_class & B);
bool operator != (const qi_class & A, const qi_class & B);
bool operator ! (const qi_class & A);

    //
    // basic functions
    //

void swap(qi_class & A, qi_class & B);

    //
    // high level functions
    //

bool generate_prime_ideal(qi_class & A, const bigint & p);

    //
    // reduction operators
    //

void apply_rho(qi_class & A, const qi_class & B);
qi_class apply_rho(const qi_class & B);

void apply_inverse_rho(qi_class & A, const qi_class & B);
qi_class apply_inverse_rho(const qi_class & B);

    //
    // subgroup computation
    //

base_vector< bigint > subgroup(base_vector< qi_class > & G);

    // subgroup computation with no knowledge of class number
base_vector< bigint > subgroup_BJT(base_vector< qi_class > & G,
				   base_vector< long > & v);

#ifndef HEADBANGER
    // subgroup computation when class number is known
base_vector< bigint > subgroup_h(base_vector< qi_class > & G);
#endif



    //
    // input/output
    //

std::istream & operator >> (std::istream & in, qi_class & A);
std::ostream & operator << (std::ostream & out, const qi_class & A);



// key function for hash tables
bigint qi_class_key(const qi_class & G);






#ifndef HEADBANGER

//
// Class: ideal_node
//
// This class is simply an instance of a qi_class together with an integer
//    index.  It is used in baby-step giant-step type algorithms as an
//    element in a hash table.
//

class ideal_node
{
private:

	qi_class A;
	long index;

public:

	ideal_node() {}
	ideal_node(const qi_class & G, const long idx)
	{
		A.assign(G);
		index = idx;
	}
	ideal_node(const qi_class & G, const bigint & idx)
	{
		A.assign(G);
		idx.longify(index);
	}
	~ideal_node() {}

	const qi_class& get_A() const
	{
		return A;
	}
	long get_index() const
	{
		return index;
	}
	void assign(const qi_class & G, const long idx)
	{
		A.assign(G);
		index = idx;
	}
	void assign(const qi_class & G, const bigint & idx)
	{
		A.assign(G);
		idx.longify(index);
	}
	void assign(const ideal_node & G)
	{
		A.assign(G.A);
		index = G.index;
	}
	ideal_node& operator = (const ideal_node & B)
	{
		A.assign(B.A);
		index = B.index;
		return *this;
	}

	void assign_zero()
	{
		A.assign_zero();
		index = 0;
	}
	bool is_zero()
	{
		return A.is_zero();
	}
	friend void swap(ideal_node & A, ideal_node & B);
	friend std::istream & operator >> (std::istream & in, ideal_node & A);
	friend std::ostream & operator << (std::ostream & out, const ideal_node & A);
};


// friend functions

void swap(ideal_node & A, ideal_node & B);
std::istream & operator >> (std::istream & in, ideal_node & A);
std::ostream & operator << (std::ostream & out, const ideal_node & A);


inline bool operator == (const ideal_node & N1, const ideal_node & N2)
{
	return (N1.get_A() == N2.get_A());
}



inline bool operator != (const ideal_node & N1, const ideal_node & N2)
{
	return (N1.get_A() != N2.get_A());
}



inline bool operator < (const ideal_node & N1, const ideal_node & N2)
{
	return (N1.get_A().get_a() < N2.get_A().get_a());
}



inline bool operator <= (const ideal_node & N1, const ideal_node & N2)
{
	return ((N1 < N2) || (N1 == N2));
}



inline bool operator > (const ideal_node & N1, const ideal_node & N2)
{
	return (N1.get_A().get_a() > N2.get_A().get_a());
}



inline bool operator >= (const ideal_node & N1, const ideal_node & N2)
{
	return ((N1 > N2) || (N1 == N2));
}



// key function for hash tables
bigint ideal_node_key(const ideal_node & G);

#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifndef LIDIA_QUADRATIC_ORDER_H_GUARD_
# include	"LiDIA/quadratic_order.h"
#endif
#ifndef LIDIA_QI_CLASS_REAL_H_GUARD_
# include	"LiDIA/qi_class_real.h"
#endif
#ifndef LIDIA_QUADRATIC_IDEAL_H_GUARD_
# include	"LiDIA/quadratic_ideal.h"
#endif
#ifndef LIDIA_QUADRATIC_FORM_H_GUARD_
# include	"LiDIA/quadratic_form.h"
#endif



#endif	// LIDIA_QI_CLASS_H_GUARD_
