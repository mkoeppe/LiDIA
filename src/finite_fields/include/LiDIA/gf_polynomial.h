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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_GF_POLYNOMIAL_H_GUARD_
#define LIDIA_GF_POLYNOMIAL_H_GUARD_

// This include file supports polynomial operations over the
// class 'gf_element'
// It is based on Stefan Neis' implementation of polynomials over
// arbitrary types.
// Arithmetic (esp. modular arithmetic) is very similar to
// arithmetic in the class Fp_polynomial.
//					Thomas Pfahler

#include "LiDIA/LiDIA.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

// The following forward declarations are necessary for the instantiation
// in fact_gf_polynomial.cc of base_factor< polynomial< gf_element > > to 
// succeed.
class gf_element;

template <typename T>
class polynomial;

template <>
class polynomial< gf_element >;

void swap(LiDIA::polynomial< gf_element > &a, polynomial< gf_element > &b);

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifndef LIDIA_FIELD_POLYNOMIAL_H_GUARD_
# include	"LiDIA/field_polynomial.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
//	SPECIAL INSTATIATION
//

template <>
class polynomial< gf_element >;
class gf_poly_modulus;
class gf_poly_multiplier;

bool checked_min_poly(polynomial< gf_element > & h,
		      const polynomial< gf_element > & g, lidia_size_t r,
		      const gf_poly_modulus& F);

template<>
gf_element base_polynomial< gf_element >::
operator() (gf_element const& value) const;

template<>
class polynomial< gf_element >
{
private:

	static const galois_field& common_field(const galois_field&,
						const galois_field&);
	// returns a common superfield of a and b
	// in this version, an error is raised if a != b

	void check_coefficients();
	// looks for an initialized K field among all pol[i].get_field()
	// (if none found, set K = ffield)
	// initialize all uninitialized elements with K
	// sets ffield to K

	friend void base_polynomial< gf_element >::set_degree(lidia_size_t);
	friend void integral(field_polynomial< gf_element > &,
			     const base_polynomial< gf_element > &);
	friend class gf_element;
	static const galois_field * FIELD;
	static void build_frame(const galois_field& K)
	{
		FIELD = &K;
	}
	static void delete_frame()
	{
		FIELD = 0;
	}
	// FIELD is needed in base_polynomial::set_degree() in order
	// to indicate in which field we are working

	field_polynomial< gf_element > pol;
	galois_field ffield;

public:

	//
	// constructors and destructor
	//

	polynomial () :
		pol(), ffield()
	{ }

	polynomial (const gf_element &x) :
		pol(x), ffield(x.get_field())
	{ }

	polynomial (const polynomial< gf_element > &p) :
		pol(p.pol), ffield(p.ffield)
	{ }

	polynomial (const galois_field &K)
	{
		this->assign_zero(K);
	}

	polynomial(const gf_element *v, lidia_size_t d) :
		pol(v, d)
	{
		check_coefficients();
	}

	polynomial(const base_vector< gf_element > &v) :
		pol(v)
	{
		check_coefficients();
	}

	~polynomial()
	{ }


	//
	// friend functions
	//

	// comparisons

	friend bool operator == (const polynomial< gf_element > &a,
				 const polynomial< gf_element > &b);

	friend bool operator != (const polynomial< gf_element > &a,
				 const polynomial< gf_element > &b);

	gf_element lead_coeff() const;

	gf_element const_term() const;


	//
	// member functions
	//

	void remove_leading_zeros();

        int set_data (const gf_element * d, lidia_size_t l)
        {
	        return pol.set_data(d, l);
	}
	
        gf_element* get_data () const
        {
	        return pol.get_data();
	}

	const galois_field &get_field() const
	{
		return ffield;
	}

	void set_field(const galois_field &K);

	lidia_size_t degree() const
	{
		return pol.degree();
	}

	void set_degree(lidia_size_t d);

	bool is_monic() const;

	bool is_zero() const
	{
		return pol.is_zero();
	}

	bool is_one() const
	{
		return pol.is_one();
	}

	bool is_x() const
	{
		return pol.is_x();
	}

	//
	// assignment
	//

	void assign(const gf_element &a);

	polynomial< gf_element > & operator = (const gf_element &a);

	void assign(const polynomial< gf_element > &a);

	polynomial< gf_element > & operator = (const polynomial< gf_element > &a)
	{
		this->assign(a);
		return *this;
	}

	void assign_zero();
	void assign_one();
	void assign_x();

	void assign_zero(const galois_field &K);
	void assign_one(const galois_field &K);
	void assign_x(const galois_field &K);

	friend void swap(polynomial< gf_element > &a, polynomial< gf_element > &b);


	//
	// arithmetic procedures
	//

	friend void negate(polynomial< gf_element > & c,
			   const polynomial< gf_element > &a);

	friend void add(polynomial< gf_element > & c,
			const polynomial< gf_element > & a,
			const polynomial< gf_element > & b);

	friend void add(polynomial< gf_element > & c,
			const polynomial< gf_element > & a, const gf_element & b);

	friend void add(polynomial< gf_element > & c,
			const gf_element & b, const polynomial< gf_element > & a);


	friend void subtract(polynomial< gf_element > & c,
			     const polynomial< gf_element > & a,
			     const polynomial< gf_element > & b);

	friend void subtract(polynomial< gf_element > & c,
			     const polynomial< gf_element > & a, const gf_element & b);

	friend void subtract(polynomial< gf_element > & c,
			     const gf_element & b, const polynomial< gf_element > & a);


	friend void multiply(polynomial< gf_element > & c,
			     const polynomial< gf_element > & a,
			     const polynomial< gf_element > & b);

	friend void multiply(polynomial< gf_element > & c,
			     const polynomial< gf_element > & a, const gf_element & b);

	friend void multiply(polynomial< gf_element > & c,
			     const gf_element & b, const polynomial< gf_element > & a);

	friend void square(polynomial< gf_element > & c,
			   const polynomial< gf_element > & a);

	friend void power(polynomial< gf_element > & c,
			  const polynomial< gf_element > & a, const bigint & b);



	//
	// operator overloading
	//

	const gf_element & operator[] (lidia_size_t i) const
	{
		return pol[i];
	}
	gf_element & operator[] (lidia_size_t i)
	{
		return pol[i];
	}

	gf_element operator() (const gf_element & value) const;

	friend
	polynomial< gf_element > operator - (const polynomial< gf_element > &a);

	friend
	polynomial< gf_element > operator + (const polynomial< gf_element > &a,
					     const polynomial< gf_element > &b);
	friend
	polynomial< gf_element > operator + (const polynomial< gf_element > & a,
					     const gf_element& b);
	friend
	polynomial< gf_element > operator + (const gf_element & b,
					     const polynomial< gf_element > & a);

	friend
	polynomial< gf_element > operator - (const polynomial< gf_element > &a,
					     const polynomial< gf_element > &b);
	friend
	polynomial< gf_element > operator - (const polynomial< gf_element > &a,
					     const gf_element &b);
	friend
	polynomial< gf_element > operator - (const gf_element &a,
					     const polynomial< gf_element > &b);

	friend
	polynomial< gf_element > operator * (const polynomial< gf_element > &a,
					     const polynomial< gf_element > &b);
	friend
	polynomial< gf_element > operator * (const polynomial< gf_element > &a,
					     const gf_element &b);
	friend
	polynomial< gf_element > operator * (const gf_element &b,
					     const polynomial< gf_element > &a);

	polynomial< gf_element > & operator += (const polynomial< gf_element > &a);
	polynomial< gf_element > & operator += (const gf_element &a);

	polynomial< gf_element > & operator -= (const polynomial< gf_element > &a);
	polynomial< gf_element > & operator -= (const gf_element &a);

	polynomial< gf_element > & operator *= (const polynomial< gf_element > &a);
	polynomial< gf_element > & operator *= (const gf_element &a);

	//
	// functions
	//
	friend void derivative(polynomial< gf_element > &,
			       const polynomial< gf_element > &);
	friend
	polynomial< gf_element > derivative(const polynomial< gf_element > & a);

	friend void integral(polynomial< gf_element > &,
			     const polynomial< gf_element > &);

	friend
	polynomial< gf_element > randomize(const galois_field &K, lidia_size_t n);


	//
	// input / output
	//

	friend std::istream & operator >> (std::istream &is, polynomial< gf_element > &a);

	std::istream & read_verbose(std::istream &is);

	friend std::ostream & operator << (std::ostream &os, const polynomial< gf_element > &a);


	//
	// Division and related stuff
	//

	friend void div_rem(polynomial< gf_element > & q, polynomial< gf_element > & r,
			    const polynomial< gf_element > & a,
			    const polynomial< gf_element > & b);

	friend
	void divide(polynomial< gf_element > & c,
		    const polynomial< gf_element > & a, const gf_element & b);

	friend
	void divide(polynomial< gf_element > & q,
		    const polynomial< gf_element > & a,
		    const polynomial< gf_element > & b);

	friend
	void remainder(polynomial< gf_element > & r,
		       const polynomial< gf_element > & a,
		       const polynomial< gf_element > & b);

	polynomial< gf_element > & operator /= (const polynomial< gf_element > &a);

	polynomial< gf_element > & operator /= (const gf_element &a);

	polynomial< gf_element > & operator %= (const polynomial< gf_element > &a);



	friend void gcd(polynomial< gf_element > &d,
			const polynomial< gf_element > &aa, const polynomial< gf_element > &bb);

	friend void
	xgcd(polynomial< gf_element > &d, polynomial< gf_element > &x,
	     polynomial< gf_element > &y, const polynomial< gf_element > &aa,
	     const polynomial< gf_element > &bb);


	//
	//  modular arithmetic without pre-conditioning
	//

	friend void multiply_mod(polynomial< gf_element > & x,
				 const polynomial< gf_element > & a, const polynomial< gf_element > & b,
				 const polynomial< gf_element > & f);

	friend void square_mod(polynomial< gf_element > & x,
			       const polynomial< gf_element > & a, const polynomial< gf_element > & f);

	friend void multiply_by_x_mod(polynomial< gf_element > & h,
				      const polynomial< gf_element > & a, const polynomial< gf_element > & f);

	friend void invert_mod(polynomial< gf_element > & x,
			       const polynomial< gf_element > & a, const polynomial< gf_element > & f);

	friend bool invert_mod_status(polynomial< gf_element > & x,
				      const polynomial< gf_element > & a, const polynomial< gf_element > & f);

	friend void power_mod(polynomial< gf_element > & h,
			      const polynomial< gf_element > & g, const bigint & e,
			      const polynomial< gf_element > & f);

	friend void power_x_mod(polynomial< gf_element > & h,
				const bigint & e, const polynomial< gf_element > & f);

	friend void power_x_plus_a_mod(polynomial< gf_element > & h,
				       const gf_element & a, const bigint & e,
				       const polynomial< gf_element > & f);

	friend void cyclic_reduce(polynomial< gf_element > & x,
				  const polynomial< gf_element > & a, lidia_size_t m);


	//
	//  'plain' arithmetic (slow)
	//
	friend void plain_power(polynomial< gf_element > & c,
				const polynomial< gf_element > & a, const bigint & b);
	friend void plain_gcd(polynomial< gf_element > &d,
			      const polynomial< gf_element > &aa, const polynomial< gf_element > &bb);
	friend void plain_multiply(polynomial< gf_element > &c,
				   const polynomial< gf_element > &a, const polynomial< gf_element > &b);
	friend void plain_square(polynomial< gf_element > & c,
				 const polynomial< gf_element > & a);
	friend void plain_div_rem(polynomial< gf_element > &q,
				  polynomial< gf_element > &r, const polynomial< gf_element > &a,
				  const polynomial< gf_element > &b);
	friend void plain_divide(polynomial< gf_element > &q,
				 const polynomial< gf_element > &a, const polynomial< gf_element > &b);
	friend void plain_remainder(polynomial< gf_element > &r,
				    const polynomial< gf_element > &a, const polynomial< gf_element > &b);

	//
	//  'fast'  arithmetic (Kronecker substitution)
	//
	friend void fast_multiply(polynomial< gf_element > & c,
				  const polynomial< gf_element > & a, const polynomial< gf_element > & b);
	friend void fast_div_rem(polynomial< gf_element > & q,
				 polynomial< gf_element > & r, const polynomial< gf_element > & a,
				 const polynomial< gf_element > & b);
	friend void invert(polynomial< gf_element > & x,
			   const polynomial< gf_element > & a, lidia_size_t m);
	friend void copy_reverse(polynomial< gf_element > & x,
				 const polynomial< gf_element > & a, lidia_size_t lo, lidia_size_t hi);

	//
	//  modular arithmetic with pre-conditioning (class gf_poly_modulus)
	//
	friend class gf_poly_modulus;
	friend void remainder(polynomial< gf_element > &c,
			      const polynomial< gf_element > &a, const gf_poly_modulus &P);
	friend void multiply(polynomial< gf_element > &c,
			     const polynomial< gf_element > &a, const polynomial< gf_element > &b,
			     const gf_poly_modulus &P);
	friend void square(polynomial< gf_element > &c,
			   const polynomial< gf_element > &a, const gf_poly_modulus &P);
	friend void invert(polynomial< gf_element > &c,
			   const polynomial< gf_element > &a, const gf_poly_modulus &P);
	friend void power(polynomial< gf_element > &c,
			  const polynomial< gf_element > &a, const bigint &e,
			  const gf_poly_modulus &P);
	friend void power_x(polynomial< gf_element > &c, const bigint &e,
			    const gf_poly_modulus &P);
	friend void power_x_plus_a(polynomial< gf_element > &c,
				   const gf_element &a, const bigint &e, const gf_poly_modulus &P);



	//
	//  stuff for factoring algorithms
	//

	friend void random_basis_elt2(polynomial< gf_element > & g,
				      const base_vector< sdigit > &D, Fp_polynomial** M,
				      const galois_field& K);
	friend bool checked_min_poly(polynomial< gf_element > & h,
				     const polynomial< gf_element > & g, lidia_size_t r,
				     const gf_poly_modulus& F);

	friend void trace_map(polynomial< gf_element > & w,
			      const polynomial< gf_element > & a, lidia_size_t d,
			      const gf_poly_modulus &f, const polynomial< gf_element > & b);
}; // end class polynomial< gf_element >



// declaration of gf_polynomial's friend functions
bool operator == (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b);
bool operator != (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b);

void swap(polynomial< gf_element > &a, polynomial< gf_element > &b);


	//
	// arithmetic procedures
	//

void negate(polynomial< gf_element > & c,
	    const polynomial< gf_element > &a);

void add(polynomial< gf_element > & c,
	 const polynomial< gf_element > & a,
	 const polynomial< gf_element > & b);

void add(polynomial< gf_element > & c,
	 const polynomial< gf_element > & a, const gf_element & b);

void add(polynomial< gf_element > & c,
	 const gf_element & b, const polynomial< gf_element > & a);


void subtract(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a,
	      const polynomial< gf_element > & b);

void subtract(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a, const gf_element & b);

void subtract(polynomial< gf_element > & c,
	      const gf_element & b, const polynomial< gf_element > & a);


void multiply(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a,
	      const polynomial< gf_element > & b);

void multiply(polynomial< gf_element > & c,
	      const polynomial< gf_element > & a, const gf_element & b);

void multiply(polynomial< gf_element > & c,
	      const gf_element & b, const polynomial< gf_element > & a);

void square(polynomial< gf_element > & c,
	    const polynomial< gf_element > & a);

void power(polynomial< gf_element > & c,
	   const polynomial< gf_element > & a, const bigint & b);


polynomial< gf_element > operator + (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b);
polynomial< gf_element > operator + (const polynomial< gf_element > & a,
				     const gf_element& b);
polynomial< gf_element > operator + (const gf_element & b,
				     const polynomial< gf_element > & a);

polynomial< gf_element > operator - (const polynomial< gf_element > &a);
polynomial< gf_element > operator - (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b);
polynomial< gf_element > operator - (const polynomial< gf_element > &a,
				     const gf_element &b);
polynomial< gf_element > operator - (const gf_element &a,
				     const polynomial< gf_element > &b);

polynomial< gf_element > operator * (const polynomial< gf_element > &a,
				     const polynomial< gf_element > &b);

polynomial< gf_element > operator * (const polynomial< gf_element > &a,
				     const gf_element &b);

polynomial< gf_element > operator * (const gf_element &b,
				     const polynomial< gf_element > &a);

	//
	// functions
	//
void derivative(polynomial< gf_element > &,
		const polynomial< gf_element > &);


void integral(polynomial< gf_element > &,
			     const polynomial< gf_element > &);


polynomial< gf_element > randomize(const galois_field &K, lidia_size_t n);


	//
	// input / output
	//

std::istream & operator >> (std::istream &is, polynomial< gf_element > &a);

std::ostream & operator << (std::ostream &os,
			    const polynomial< gf_element > &a);


	//
	// Division and related stuff
	//

void div_rem(polynomial< gf_element > & q, polynomial< gf_element > & r,
	     const polynomial< gf_element > & a,
	     const polynomial< gf_element > & b);


void divide(polynomial< gf_element > & c,
	    const polynomial< gf_element > & a, const gf_element & b);

void divide(polynomial< gf_element > & q,
	    const polynomial< gf_element > & a,
	    const polynomial< gf_element > & b);

void remainder(polynomial< gf_element > & r,
	       const polynomial< gf_element > & a,
	       const polynomial< gf_element > & b);


void gcd(polynomial< gf_element > &d,
	 const polynomial< gf_element > &aa,
	 const polynomial< gf_element > &bb);

void xgcd(polynomial< gf_element > &d, polynomial< gf_element > &x,
	  polynomial< gf_element > &y, const polynomial< gf_element > &aa,
	  const polynomial< gf_element > &bb);


	//
	//  modular arithmetic without pre-conditioning
	//

void multiply_mod(polynomial< gf_element > & x,
		  const polynomial< gf_element > & a,
		  const polynomial< gf_element > & b,
		  const polynomial< gf_element > & f);

void square_mod(polynomial< gf_element > & x,
		const polynomial< gf_element > & a,
		const polynomial< gf_element > & f);

void multiply_by_x_mod(polynomial< gf_element > & h,
		       const polynomial< gf_element > & a,
		       const polynomial< gf_element > & f);

void invert_mod(polynomial< gf_element > & x,
		const polynomial< gf_element > & a,
		const polynomial< gf_element > & f);

bool invert_mod_status(polynomial< gf_element > & x,
		       const polynomial< gf_element > & a,
		       const polynomial< gf_element > & f);

void power_mod(polynomial< gf_element > & h,
	       const polynomial< gf_element > & g, const bigint & e,
	       const polynomial< gf_element > & f);

void power_x_mod(polynomial< gf_element > & h,
		 const bigint & e, const polynomial< gf_element > & f);

void power_x_plus_a_mod(polynomial< gf_element > & h,
			const gf_element & a, const bigint & e,
			const polynomial< gf_element > & f);
    
void cyclic_reduce(polynomial< gf_element > & x,
		   const polynomial< gf_element > & a, lidia_size_t m);


	//
	//  'plain' arithmetic (slow)
	//
void plain_power(polynomial< gf_element > & c,
		 const polynomial< gf_element > & a, const bigint & b);
void plain_gcd(polynomial< gf_element > &d,
	       const polynomial< gf_element > &aa,
	       const polynomial< gf_element > &bb);
void plain_multiply(polynomial< gf_element > &c,
		    const polynomial< gf_element > &a,
		    const polynomial< gf_element > &b);
void plain_square(polynomial< gf_element > & c,
		  const polynomial< gf_element > & a);
void plain_div_rem(polynomial< gf_element > &q,
		   polynomial< gf_element > &r,
		   const polynomial< gf_element > &a,
		   const polynomial< gf_element > &b);
void plain_divide(polynomial< gf_element > &q,
		  const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b);
void plain_remainder(polynomial< gf_element > &r,
		     const polynomial< gf_element > &a,
		     const polynomial< gf_element > &b);

	//
	//  'fast'  arithmetic (Kronecker substitution)
	//
void fast_multiply(polynomial< gf_element > & c,
		   const polynomial< gf_element > & a,
		   const polynomial< gf_element > & b);
void fast_div_rem(polynomial< gf_element > & q,
		  polynomial< gf_element > & r,
		  const polynomial< gf_element > & a,
		  const polynomial< gf_element > & b);
void invert(polynomial< gf_element > & x,
	    const polynomial< gf_element > & a, lidia_size_t m);
void copy_reverse(polynomial< gf_element > & x,
		  const polynomial< gf_element > & a,
		  lidia_size_t lo, lidia_size_t hi);

	//
	//  modular arithmetic with pre-conditioning (class gf_poly_modulus)
	//
void remainder(polynomial< gf_element > &c,
	       const polynomial< gf_element > &a, const gf_poly_modulus &P);
void multiply(polynomial< gf_element > &c,
	      const polynomial< gf_element > &a,
	      const polynomial< gf_element > &b,
	      const gf_poly_modulus &P);
void square(polynomial< gf_element > &c,
	    const polynomial< gf_element > &a, const gf_poly_modulus &P);
void invert(polynomial< gf_element > &c,
	    const polynomial< gf_element > &a, const gf_poly_modulus &P);
void power(polynomial< gf_element > &c,
	   const polynomial< gf_element > &a, const bigint &e,
	   const gf_poly_modulus &P);
void power_x(polynomial< gf_element > &c, const bigint &e,
	     const gf_poly_modulus &P);
void power_x_plus_a(polynomial< gf_element > &c,
		    const gf_element &a, const bigint &e,
		    const gf_poly_modulus &P);



	//
	//  stuff for factoring algorithms
	//

void random_basis_elt2(polynomial< gf_element > & g,
		       const base_vector< sdigit > &D, Fp_polynomial** M,
		       const galois_field& K);
bool checked_min_poly(polynomial< gf_element > & h,
		      const polynomial< gf_element > & g, lidia_size_t r,
		      const gf_poly_modulus& F);

void trace_map(polynomial< gf_element > & w,
	       const polynomial< gf_element > & a, lidia_size_t d,
	       const gf_poly_modulus &f, const polynomial< gf_element > & b);



typedef polynomial< gf_element > gf_polynomial;

//for class single_factor < gf_polynomial >
bool operator< (const polynomial < gf_element > &a,
		const polynomial< gf_element > &b);
bool operator <= (const polynomial< gf_element > &a,
		  const polynomial< gf_element > &b);



base_vector< gf_element > find_roots(const gf_polynomial &f);
gf_element find_root( const gf_polynomial & f );


class gf_poly_modulus
{
	gf_polynomial F_plain;
	Fp_polynomial F, F_recip;
	bool use_Fp; // flag indicating whether F, F_recip should be used

public:
	gf_poly_modulus() { }
	gf_poly_modulus(const gf_poly_modulus &p);
	gf_poly_modulus(const gf_polynomial &f);

	void build(const gf_polynomial &f);

	const gf_polynomial& modulus() const
	{
		return F_plain;
	}

	void rem21(gf_polynomial &c, const gf_polynomial &a) const;
	friend void remainder(gf_polynomial &c, const gf_polynomial &a,
			      const gf_poly_modulus &P);
	friend void multiply(gf_polynomial &c, const gf_polynomial &a,
			     const gf_polynomial &b, const gf_poly_modulus &P);
	friend void square(gf_polynomial &c, const gf_polynomial &a,
			   const gf_poly_modulus &P);
	friend void invert(gf_polynomial &c, const gf_polynomial &a,
			   const gf_poly_modulus &P);
	friend void power(gf_polynomial &c, const gf_polynomial &a,
			  const bigint &e, const gf_poly_modulus &P);
	friend void power_x(gf_polynomial &c, const bigint &e,
			    const gf_poly_modulus &P);
	friend void power_x_plus_a(gf_polynomial &c, const gf_element &a,
				   const bigint &e, const gf_poly_modulus &P);
};



#define LIDIA_CLASS_GF_POLYNOMIAL



// must be inline
inline std::istream &
operator >> (std::istream &is, polynomial< gf_element > &a)
{
	if (a.ffield.degree() == 0)
		lidia_error_handler("polynomial< gf_element >",
				    "operator >>::polynomial must be assigned to a field before any input");
	polynomial< gf_element >::build_frame(a.ffield);
	is >> a.pol;
	polynomial< gf_element >::delete_frame();
	a.remove_leading_zeros();
	a.check_coefficients();
	return is;
}



// must be inline
inline std::ostream &
operator << (std::ostream &os, const polynomial< gf_element > &a)
{
	polynomial< gf_element >::build_frame(a.ffield);
	a.pol.print_verbose(os, 'y');
	polynomial< gf_element >::delete_frame();
	return os;
}


void to_Kronecker(Fp_polynomial &g, const gf_polynomial &f, lidia_size_t lo, lidia_size_t hi);
void from_Kronecker(gf_polynomial &f, const Fp_polynomial &g, lidia_size_t lo, lidia_size_t hi);



class gf_poly_argument
{
	gf_polynomial *vec;
	lidia_size_t len;

	void inner_prod(gf_polynomial &x, const gf_polynomial &g,
			lidia_size_t lo, lidia_size_t hi) const;

public:
	gf_poly_argument();
	gf_poly_argument(const gf_poly_argument &x);
	~gf_poly_argument();

	void build(const gf_polynomial &h, const gf_poly_modulus &F,
		   lidia_size_t m);
	//computes and stores h, h^2, ..., h^m mod f

	void compose(gf_polynomial &x, const gf_polynomial &g,
		     const gf_poly_modulus &F) const;
	//x = g(h) mod F
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#include	"LiDIA/specialization/gf_polynomial.special"



#endif	// LIDIA_GF_POLYNOMIAL_H_GUARD_
