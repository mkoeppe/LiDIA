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


#ifndef LIDIA_SF_FP_POLYNOMIAL_H_GUARD_
#define LIDIA_SF_FP_POLYNOMIAL_H_GUARD_


#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//MM, Added include of single_factor.h, because general class must be
//    defined prior to an explicit class instance.
//    Self initiating compilers, e.g. HP UX 10.20 CC, eventually include
//    sf_polynomial.h before single_factor.h.



//****************************************************************************
//
//     USE THIS FILE FOR SPECIAL INSTANTIATIONS OF CLASS SINGLE_FACTOR
//
//****************************************************************************
//
//  0) Backup the files
//  		single_factor.h and
//  		single_factor.cc
//
//  1) Replace 'Fp_polynomial' by the name of the class for which you
//     want to instantiate the class single_factor
//
//  2) If necessary, declare new members in this class
//     If you have declared new members, you have to adjust the following
//     functions:
//     - all constructors, destructor
//     - operator =
//     - friend void swap(single_factor< * > &a, single_factor< * > &b)
//     - operator << and >>
//     Otherwise, you may have trouble because your 'new members' might be
//     uninitialized.
//
//  3) Implement your own primality test in the function
//     single_factor< * >::is_prime_factor()
//
//  4) If you want to get rid of 'units' (i.e. invertible elements of your
//     data type), rewrite the function single_factor< * >::extract_unit().
//     This function is used in the class factorization< * > in order to
//     collect such units separately from 'normal' (i.e. non-invertible)
//     factors.
//
//  5) Implement your own standard factorization algorithm in the function
//     single_factor< * >::factor(). Remember that the return type must be
//     factorization< * > .
//
//  6) Implement your own factorization algorithms, e.g.
//     factorization< * > single_factor< * >::split_into_two_factors()
//     factorization< * > single_factor< * >::complete_factorization()
//     factorization< * > single_factor< * >::trial_division()
//     ... and so on ...
//
//
//  For examples, see the instantiations of
//  	single_factor< bigint >
//  	single_factor< Fp_polynomial >
//  	single_factor< ideal >
//
//****************************************************************************



template< class T > class base_factor;
//template < class T > class single_factor;
template< class T > class factorization;


//****************************************************************
//                   class single_factor< Fp_polynomial >
//****************************************************************

template<>
class single_factor< Fp_polynomial > : public base_factor< Fp_polynomial >
{
	friend class factorization< Fp_polynomial >;


private:

	//additional information should be defined here

	bool know_modulus;
	//flag indicating whether the field over which the polynomial
	//is defined is known or not
	//if it is not known, (*this) was not initialized yet
	//and we assume 'rep' to equal '1'

	static bool verbose_flag;


public:

	inline static bool verbose()
	{
		return verbose_flag;
	}
	inline static void set_verbose_mode(int i)
	{
		verbose_flag = (i != 0);
	}


public:

	//
	// constructors, destructor
	//

	single_factor();
	//default value must be '1' (the neutral element of multiplication)
	//^^^^^^^^^^^^^^^^^^^^^^^^^
	//we have a problem here : we must assign '1' to a polynomial without
	//knowing the field the polynomial is defined over
	//therefore we use the additional flag 'know_modulus'

	single_factor(const single_factor< Fp_polynomial > &);
	//copy constructor
	single_factor(const Fp_polynomial&);
	//conversion for type Fp_polynomial
	~single_factor();
	//destructor


	//
	// swap
	//

public :

	friend void swap(single_factor< Fp_polynomial > & a, single_factor< Fp_polynomial > & b);


	//
	// assignment
	//
	single_factor< Fp_polynomial > & operator = (const single_factor< Fp_polynomial > & x);
	single_factor< Fp_polynomial > & operator = (const Fp_polynomial & x);
#ifndef HEADBANGER
	void assign(const single_factor< Fp_polynomial > & x);
	void assign(const Fp_polynomial & x);
#endif	// HEADBANGER

	//
	// I/O
	//
	friend std::istream & operator >> (std::istream &in, single_factor< Fp_polynomial > &f);

	//friend std::ostream & operator << (std::ostream &out, const single_factor < Fp_polynomial > &f)
	// defined for base class 'base_factor'


	//
	// queries
	//

	bool is_one() const
	{
		return (know_modulus ? rep.is_one() : 1);
	}

	bool is_prime_factor() const
	{
		return (prime_flag() == prime);
	}

	bool is_prime_factor(int test);
	//test == 0 ->no explicit primality test (only a flag is checked)
	//test != 0 ->explicit prime-test if prime_state() == unknown
	//see manual

	Fp_polynomial extract_unit();
	friend lidia_size_t ord_divide(const single_factor< Fp_polynomial > &a,
				       single_factor< Fp_polynomial > &b);



	//
	// factorization algorithms
	//

	// "do not implement function bodies in .h-file"
	// "otherwise, it will confuse the CC-compiler"

	factorization< Fp_polynomial > factor() const;
	// standard factorization algorithm, used in function
	// factorization < Fp_polynomial >::factor_all_components
	friend factorization< Fp_polynomial >
	factor(const single_factor< Fp_polynomial > & f);


	friend factorization< Fp_polynomial > square_free_decomp(const Fp_polynomial& f);
	factorization< Fp_polynomial > square_free_decomp() const;


	friend factorization< Fp_polynomial > sf_berlekamp(const Fp_polynomial& f);
	factorization< Fp_polynomial > sf_berlekamp() const;
	// Assumes f is square-free, monic and f(0) != 0

	friend factorization< Fp_polynomial > berlekamp(const Fp_polynomial& f);
	factorization< Fp_polynomial > berlekamp() const;
	// f must be non-zero


	friend factorization< Fp_polynomial > ddf(const Fp_polynomial& f);
	factorization< Fp_polynomial > ddf() const;
	// f must be square-free and monic

	friend factorization< Fp_polynomial > edf(const Fp_polynomial &f,
						  lidia_size_t d);
	factorization< Fp_polynomial > edf(lidia_size_t d) const;
	// f must be the product of distinc monic irred. polynomials of degree d

	friend factorization< Fp_polynomial > sf_can_zass(const Fp_polynomial& f);
	factorization< Fp_polynomial > sf_can_zass() const;
	// Assumes f is square-free and monic

	friend factorization< Fp_polynomial > can_zass(const Fp_polynomial& f);
	factorization< Fp_polynomial > can_zass() const;
	// f must be non-zero


	//
	// misc
	//

	//we need new comparisons because of the flag 'know_modulus'
	friend bool operator == (const single_factor< Fp_polynomial > & a,
				 const single_factor< Fp_polynomial > & b);
	friend bool operator != (const single_factor< Fp_polynomial > & a,
				 const single_factor< Fp_polynomial > & b);

	// need special versions of the following because of the flag 'know_modulus':
#ifndef HEADBANGER
	friend void multiply(single_factor< Fp_polynomial > &c,
			     const single_factor< Fp_polynomial > &a,
			     const single_factor< Fp_polynomial > &b);
	friend void divide(single_factor< Fp_polynomial > &c,
			   const single_factor< Fp_polynomial > &a,
			   const single_factor< Fp_polynomial > &b);
#endif	// HEADBANGER
	friend void gcd(single_factor< Fp_polynomial > &c,
			const single_factor< Fp_polynomial > &a,
			const single_factor< Fp_polynomial > &b);
	friend single_factor< Fp_polynomial > operator * (const single_factor< Fp_polynomial > &a,
							  const single_factor< Fp_polynomial > &b);
	friend single_factor< Fp_polynomial > operator / (const single_factor< Fp_polynomial > &a,
							  const single_factor< Fp_polynomial > &b);
};

// friend functions

void swap(single_factor< Fp_polynomial > & a,
	  single_factor< Fp_polynomial > & b);

    //
    // I/O
    //
std::istream & operator >> (std::istream &in,
			    single_factor< Fp_polynomial > &f);

    //
    // queries
    //

lidia_size_t ord_divide(const single_factor< Fp_polynomial > &a,
			single_factor< Fp_polynomial > &b);

    //
    // factorization algorithms
    //

factorization< Fp_polynomial >
factor(const single_factor< Fp_polynomial > & f);

factorization< Fp_polynomial > square_free_decomp(const Fp_polynomial& f);

factorization< Fp_polynomial > sf_berlekamp(const Fp_polynomial& f);
    // Assumes f is square-free, monic and f(0) != 0
factorization< Fp_polynomial > berlekamp(const Fp_polynomial& f);
    // f must be non-zero

factorization< Fp_polynomial > ddf(const Fp_polynomial& f);
    // f must be square-free and monic

factorization< Fp_polynomial > edf(const Fp_polynomial &f,
				   lidia_size_t d);
    // f must be the product of distinc monic irred. polynomials of degree d

factorization< Fp_polynomial > sf_can_zass(const Fp_polynomial& f);
    // Assumes f is square-free and monic

factorization< Fp_polynomial > can_zass(const Fp_polynomial& f);
    // f must be non-zero

    //
    // misc
    //

    //we need new comparisons because of the flag 'know_modulus'
bool operator == (const single_factor< Fp_polynomial > & a,
		  const single_factor< Fp_polynomial > & b);
inline bool operator != (const single_factor< Fp_polynomial > & a,
		  const single_factor< Fp_polynomial > & b)
{
    return !(a == b);
}

// need special versions of the following because of the flag 'know_modulus':
#ifndef HEADBANGER
void multiply(single_factor< Fp_polynomial > &c,
	      const single_factor< Fp_polynomial > &a,
	      const single_factor< Fp_polynomial > &b);
void divide(single_factor< Fp_polynomial > &c,
	    const single_factor< Fp_polynomial > &a,
	    const single_factor< Fp_polynomial > &b);
#endif	// HEADBANGER
void gcd(single_factor< Fp_polynomial > &c,
	 const single_factor< Fp_polynomial > &a,
	 const single_factor< Fp_polynomial > &b);

single_factor< Fp_polynomial >
operator * (const single_factor< Fp_polynomial > &a,
	    const single_factor< Fp_polynomial > &b);

single_factor< Fp_polynomial >
operator / (const single_factor< Fp_polynomial > &a,
	    const single_factor< Fp_polynomial > &b);




inline std::istream &
operator >> (std::istream &in, single_factor< Fp_polynomial > &f)
{
	in >> f.rep;
	f.set_prime_flag(decomposable_object::unknown);
	f.know_modulus = true;
	return in;
}



inline single_factor< Fp_polynomial >
operator * (const single_factor< Fp_polynomial > &a,
	    const single_factor< Fp_polynomial > &b)
{
	single_factor< Fp_polynomial > tmp;

	multiply(tmp, a, b);
	return tmp;
}



inline single_factor< Fp_polynomial >
operator / (const single_factor< Fp_polynomial > &a,
	    const single_factor< Fp_polynomial > &b)
{
	single_factor< Fp_polynomial > tmp;

	divide(tmp, a, b);
	return tmp;
}



void swap(single_factor< Fp_polynomial > & a, single_factor< Fp_polynomial > & b);


void factor_quadratic_pol(factorization< Fp_polynomial > & F, const Fp_polynomial &f);
void sf_berlekamp_work(factorization< Fp_polynomial > &factors,
		       const Fp_polynomial &x_to_the_p, const Fp_poly_modulus &F);
void sf_can_zass_work(factorization< Fp_polynomial > &factors,
		      const Fp_polynomial &x_to_the_p, const Fp_poly_modulus &F);


void find_irred_factors(factorization< Fp_polynomial > &factors,
			const Fp_polynomial & f, const Fp_polynomial & g,
			const base_vector< bigint > &roots);

void append_irred_factor(factorization< Fp_polynomial > &F,
			 const Fp_polynomial &f, lidia_size_t i = 1);

factorization< Fp_polynomial > factor_binomial(const Fp_polynomial &f,
					       int FACTOR_P_MINUS_ONE = 0);
void factor(factorization< Fp_polynomial > &F, const Fp_polynomial &f);

void factor_generic(factorization< Fp_polynomial > &factors,
		    const Fp_polynomial &ff,
		    void (*do_work)(factorization< Fp_polynomial > &, const Fp_polynomial &));



// ddf  
void old_ddf(factorization< Fp_polynomial > & F, const Fp_polynomial& f, const Fp_polynomial& h);
void ddf(factorization< Fp_polynomial > &F, const Fp_polynomial& f, const Fp_polynomial& h);
// Performs distinct-degree factorization.
// h = X^p mod f
// the exponents of 'F' are the degrees of the irred. factors of 'f' !!!


// edf  
void edf(factorization< Fp_polynomial > &factors, const Fp_poly_modulus& F, const Fp_polynomial& b, lidia_size_t d);
// Performs equal-degree factorization.
// b = X^p (mod f = F.polynomial()).
// d = degree of irreducible factors of f
// Space for the trace-map computation can be controlled via ComposeBound.

void edf1(Fp_polynomial& factor, const Fp_poly_modulus& F, const Fp_polynomial& b,
          lidia_size_t d);
// like edf, but returns just a single factor

void root_edf(factorization< Fp_polynomial > & factors, const Fp_polynomial& f);
// edf for d == 1

void mystery_edf(Fp_polynomial & factor, const Fp_poly_modulus & F,
		 const Fp_polynomial & b);
// assumes b = X^p mod f and p > n
// this a probabilistic algorithm that computes a monic factor of f
// It never returns 1, and with very high probability (at least for
// large p) it actually returns an irreducible factor.




//in these functions, the work is done:
void square_free_decomp(factorization< Fp_polynomial > &u, const Fp_polynomial& f);
void sf_berlekamp(factorization< Fp_polynomial > &F, const Fp_polynomial& f);
void berlekamp(factorization< Fp_polynomial > &F, const Fp_polynomial& f);
void can_zass(factorization< Fp_polynomial > &F, const Fp_polynomial& f);

extern int OLD_DDF_GCD_BLOCKING_FACTOR;
//see old_DDF.cc
extern int DDF_GCD_BLOCKING_FACTOR;
//see DDF.cc



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SF_FP_POLYNOMIAL_H_GUARD_
