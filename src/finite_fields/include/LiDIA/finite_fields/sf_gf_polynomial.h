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


#ifndef LIDIA_SF_GF_POLYNOMIAL_H_GUARD_
#define LIDIA_SF_GF_POLYNOMIAL_H_GUARD_



#ifndef LIDIA_GF_POLYNOMIAL_H_GUARD_
# include	"LiDIA/gf_polynomial.h"
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


//***************************************************************************
//
//    USE THIS FILE FOR SPECIAL INSTANTIATIONS OF CLASS SINGLE_FACTOR
//
//***************************************************************************

//  0) Backup the files
//  		single_factor.h and
//  		single_factor.c
//
//  1) Replace 'gf_polynomial' by the name of the class for which you
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
//  	single_factor< gf_polynomial >
//  	single_factor< ideal >
//
//***************************************************************************



template< class T > class base_factor;
//template < class T > class single_factor;
template< class T > class factorization;


//****************************************************************
//                 class single_factor< gf_polynomial >
//****************************************************************

template<>
class single_factor< gf_polynomial > : public base_factor< gf_polynomial >
{
	friend class factorization< gf_polynomial >;


private:
	//additional information should be defined here

	bool know_field;
	//flag indicating whether the field over which the polynomial
	//is defined is known or not
	//if it is not known, (*this) was not initialized yet
	//and we assume 'rep' to equal '1'

	static bool verbose_flag;


public:

	static bool verbose()
	{
		return verbose_flag;
	}

	static void set_verbose_mode(int i)
	{
		verbose_flag = i;
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
	//therefore we use the additional flag 'know_field'

	single_factor(const single_factor< gf_polynomial > &);
	//copy constructor
	single_factor(const gf_polynomial&);
	//conversion for type gf_polynomial
	~single_factor();
	//destructor


	//
	// swap
	//
public :
	friend void swap(single_factor< gf_polynomial > & a, single_factor< gf_polynomial > & b);



	//
	// assignment
	//
	single_factor< gf_polynomial > & operator = (const single_factor< gf_polynomial > & x);
	single_factor< gf_polynomial > & operator = (const gf_polynomial & x);
#ifndef HEADBANGER
	void assign(const single_factor< gf_polynomial > & x);
	void assign(const gf_polynomial & x);
#endif		// HEADBANGER

	//
	// I/O
	//
	friend std::istream & operator >> (std::istream &in, single_factor< gf_polynomial > &f);

	//friend std::ostream & operator << (std::ostream &out, const single_factor < gf_polynomial > &f)
	// defined for base class 'base_factor'


	//
	// queries
	//

	bool is_one() const
	{
		return (know_field ? rep.is_one() : 1);
	}

	bool is_prime_factor() const
	{
		return (prime_flag() == prime);
	}

	bool is_prime_factor(int test);
	//test == 0 ->no explicit primality test (only a flag is checked)
	//test != 0 ->explicit prime-test if prime_state() == unknown
	//see manual

	gf_polynomial extract_unit();
	friend lidia_size_t ord_divide(
		const single_factor< gf_polynomial > &a,
		single_factor< gf_polynomial > &b);



	//
	// factorization algorithms
	//

	// "do not implement function body in .h-file"
	// "otherwise, it will confuse the CC-compiler"

	factorization< gf_polynomial > factor() const;
	// standard factorization algorithm, used in function
	// factorization < gf_polynomial >::factor_all_components
	friend factorization< gf_polynomial >
	factor(const single_factor< gf_polynomial > & f);

	friend factorization< gf_polynomial > square_free_decomp(const gf_polynomial& f);
	factorization< gf_polynomial > square_free_decomp() const;

	friend factorization< gf_polynomial > sf_berlekamp(const gf_polynomial& f);
	factorization< gf_polynomial > sf_berlekamp() const;
	// Assumes f is square-free, monic and f(0) != 0

	friend factorization< gf_polynomial > berlekamp(const gf_polynomial& f);
	factorization< gf_polynomial > berlekamp() const;
	// f must be non-zero


	friend factorization< gf_polynomial > ddf(const gf_polynomial& f);
	factorization< gf_polynomial > ddf() const;
	// f must be square-free and monic

	friend factorization< gf_polynomial > edf(const gf_polynomial& f, lidia_size_t d);
	factorization< gf_polynomial > edf(lidia_size_t d) const;
	// f must be square-free and monic


	friend factorization< gf_polynomial > sf_can_zass(const gf_polynomial& f);
	factorization< gf_polynomial > sf_can_zass() const;
	// Assumes f is square-free and monic

	friend factorization< gf_polynomial > can_zass(const gf_polynomial& f);
	factorization< gf_polynomial > can_zass() const;
	// f must be non-zero



	//
	// misc
	//

	//we need new comparisons because of the flag 'know_field'
	friend bool operator == (const single_factor< gf_polynomial > & a,
				 const single_factor< gf_polynomial > & b);
	friend bool operator != (const single_factor< gf_polynomial > & a,
				 const single_factor< gf_polynomial > & b)
	{
		return !(a == b);
	}

	// need special versions of the following because of the flag 'know_field':
#ifndef HEADBANGER
	friend void multiply(single_factor< gf_polynomial > &c,
			     const single_factor< gf_polynomial > &a,
			     const single_factor< gf_polynomial > &b);
	friend void divide(single_factor< gf_polynomial > &c,
			   const single_factor< gf_polynomial > &a,
			   const single_factor< gf_polynomial > &b);
#endif	//HEADBANGER
	friend void gcd(single_factor< gf_polynomial > &c,
			const single_factor< gf_polynomial > &a,
			const single_factor< gf_polynomial > &b);
};



inline std::istream &
operator >> (std::istream &in, single_factor< gf_polynomial > &f)
{
	in >> f.rep;
	f.set_prime_flag(decomposable_object::unknown);
	f.know_field = true;
	return in;
}



inline single_factor< gf_polynomial >
operator * (const single_factor< gf_polynomial > &a,
	    const single_factor< gf_polynomial > &b)
{
	single_factor< gf_polynomial > tmp;

	multiply(tmp, a, b);
	return tmp;
}



inline single_factor< gf_polynomial >
operator / (const single_factor< gf_polynomial > &a,
	    const single_factor< gf_polynomial > &b)
{
	single_factor< gf_polynomial > tmp;

	divide(tmp, a, b);
	return tmp;
}



void swap(single_factor< gf_polynomial > & a, single_factor< gf_polynomial > & b);


void find_irred_factors(factorization< gf_polynomial > &factors,
			const gf_polynomial & f, const gf_polynomial & g,
			const base_vector< gf_element > &roots);

void append_irred_factor(factorization< gf_polynomial > &F, const gf_polynomial &f, lidia_size_t e = 1);

void factor(factorization< gf_polynomial > &F, const gf_polynomial &f);


void ddf(factorization< gf_polynomial > & F, const gf_polynomial& f, const gf_polynomial& h);
// Performs distinct-degree factorization.
// h = X^p mod f
// the exponents of 'F' are the degrees of the irred. factors of 'f' !!!


void edf(factorization< gf_polynomial > &factors, const Fp_poly_modulus& F, const gf_polynomial& b, lidia_size_t d);
// Performs equal-degree factorization.
// b = X^p mod f.
// d = degree of irreducible factors of f


void square_free_decomp(factorization< gf_polynomial > &u, const gf_polynomial& f);
void berlekamp(factorization< gf_polynomial > &F, const gf_polynomial& f);
void can_zass(factorization< gf_polynomial > &F, const gf_polynomial& f);



bool det_irred_test(const gf_polynomial &f);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SF_GF_POLYNOMIAL_H_GUARD_
