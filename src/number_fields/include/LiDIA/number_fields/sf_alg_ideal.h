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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SF_ALG_IDEAL_H_GUARD_
#define LIDIA_SF_ALG_IDEAL_H_GUARD_


#ifndef LIDIA_ALG_NUMBER_H_GUARD_
# include	"LiDIA/alg_number.h"
#endif
#ifndef LIDIA_PRIME_IDEAL_H_GUARD_
# include	"LiDIA/prime_ideal.h"
#endif
#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class rational_factorization;

//***************************************************************************
//
//   USE THIS FILE FOR SPECIAL INSTANTIATIONS OF CLASS SINGLE_FACTOR
//
//***************************************************************************

//  0) Backup the files
//                  single_factor.h and
//                  single_factor.c
//
//  1) Replace 'alg_ideal' by the name of the class for which you
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
//  For examples, see the instantiations of
//          single_factor< bigint >
//          single_factor< Fp_polynomial >
//          single_factor< alg_ideal >
//
//***************************************************************************

template< class T > class factorization;

//****************************************************************
//                 class single_factor< alg_ideal >
//****************************************************************



typedef single_factor< alg_ideal > sf_alg_ideal;



template<>
class single_factor< alg_ideal > : public decomposable_object
{
	friend class factorization< alg_ideal >;


private:
	//additional information should be defined here

	// additional information needed for prime ideals.
	// for "normal ideals" this is always zero.
	prime_ideal rep_p;
	alg_ideal rep_a;

	static bool verbose_flag;

public:
	static bool verbose()
	{
		return verbose_flag;
	}
	static void set_verbose_mode(bool b)
	{
		verbose_flag = b;
	}

public:

	//
	// constructors, destructor
	//

	single_factor();
	// default value must be '1' (the neutral element of multiplication)
	// ^^^^^^^^^^^^^^^^^^^^^^^^^

	single_factor(const single_factor< alg_ideal > &);
	//copy constructor
	single_factor(const alg_ideal&);
	//conversion for type alg_ideal
	single_factor(const prime_ideal &);
	//conversion for type prime_ideal
	~single_factor() {}
	//destructor

	//
	// queries
	//

	bool is_prime_factor() const
	{
		return (prime_flag() == decomposable_object::prime);
	}

	bool is_prime_factor(int test);
	// test == 0 ->no explicit primality test (only a flag is checked)
	// test != 0 ->explicit prime-test if prime_state() == unknown
	//see manual

	alg_ideal extract_unit();

	bool is_one() const
	{
		return (!is_prime_factor() && rep_a.is_one());
	}

	//
	// access function
	//

	const alg_ideal& base() const
	{
		if (is_prime_factor()) {
			lidia_error_handler("single_factor< alg_ideal >",
					    "component was prime, but "
					    "adressed as non-prime");
		}

		return rep_a;
	}

	alg_ideal& base ()
	{
		if (is_prime_factor()) {
			lidia_error_handler("single_factor< alg_ideal >",
					    "component was prime, but "
					    "adressed as non-prime");
		}

		return rep_a;
	}

	const prime_ideal&  prime_base() const
	{
    	        if (!is_prime_factor())  {
		         lidia_error_handler("single_factor< alg_ideal >",
					     "component was not prime, but "
					     "adressed as prime");
		}

	    return rep_p;
	}

	prime_ideal&  prime_base()
	{
  	        if (!is_prime_factor()) {
			lidia_error_handler("single_factor< alg_ideal >",
					    "component was not prime, but "
					    "adressed as prime");
		}
		
		return rep_p;
	}

	//
	// swap
	//

public:
	friend void swap(single_factor< alg_ideal > & a, single_factor< alg_ideal > & b);

	//
	// assignment
	//

	single_factor< alg_ideal > & operator = (const single_factor< alg_ideal > & x);
	single_factor< alg_ideal > & operator = (const alg_ideal & x);
	single_factor< alg_ideal > & operator = (const prime_ideal & x);
#ifndef HEADBANGER
	void assign(const single_factor< alg_ideal > & x);
	void assign(const alg_ideal & x);
	void assign(const prime_ideal & x);
#endif	// HEADBANGER

	//
	// arithmetic operations
	//

#ifndef HEADBANGER
	friend void multiply(single_factor< alg_ideal > & c,
			     const single_factor< alg_ideal > & a,
			     const single_factor< alg_ideal > & b);

	friend void divide(single_factor< alg_ideal > & c,
			   const single_factor< alg_ideal > & a,
			   const single_factor< alg_ideal > & b);
	// ********************************************************************
	// *WARNING : the software cannot check if 'b' is a divisor of 'a' !!!*
	// ********************************************************************
#endif

	friend lidia_size_t ord_divide(const single_factor< alg_ideal > &a,
				       single_factor< alg_ideal > &b);



	//
	// factorization algorithms
	//

	factorization< alg_ideal > factor(int upper_bound = 34) const;
	// standard factorization algorithm, used in function
	// factorization < alg_ideal >::factor_all_components
	void factor(factorization< alg_ideal > &, int upper_bound = 34) const;
	friend void factor(factorization< alg_ideal > &, const alg_ideal &,
			   int upper_bound);
	friend factorization< alg_ideal > factor(const alg_ideal &,
						 int upper_bound);

	friend void decompose_prime(prime_ideal * &factor, lidia_size_t & num,
				    const bigint & p, const order & O);
	factorization< alg_ideal > finish(rational_factorization &) const;
	friend factorization< alg_ideal > finish(const alg_ideal &,
						 rational_factorization &);
	void finish(factorization< alg_ideal > &,
		    rational_factorization &) const;
	friend void finish(factorization< alg_ideal > &, const alg_ideal &,
			   rational_factorization &);

	factorization< alg_ideal > trialdiv(unsigned int upper_bound = 1000000,
					    unsigned int lower_bound = 1) const;
	void trialdiv(factorization< alg_ideal > &,
		      unsigned int upper_bound = 1000000,
		      unsigned int lower_bound = 1) const;
	friend factorization< alg_ideal > trialdiv(const alg_ideal &,
						   unsigned int upper_bound,
						   unsigned int lower_bound );
	friend void trialdiv(factorization< alg_ideal > &, const alg_ideal &,
			     unsigned int upper_bound,
			     unsigned int lower_bound);

	factorization< alg_ideal > ecm(int upper_bound = 34,
				       int lower_bound = 6, int step = 3) const;
	friend factorization< alg_ideal > ecm(const alg_ideal &,
					      int upper_bound,
					      int lower_bound, int step);
	void ecm(factorization< alg_ideal > &,
		 int upper_bound = 34, int lower_bound = 6, int step = 3) const;
	friend void ecm(factorization< alg_ideal > &, const alg_ideal &,
			int upper_bound, int lower_bound, int step);

	factorization< alg_ideal > mpqs() const;
	friend factorization< alg_ideal > mpqs(const alg_ideal &);
	void mpqs(factorization< alg_ideal > &) const;
	friend void mpqs(factorization< alg_ideal > &, const alg_ideal &);

	//
	// misc
	//

	//we need new comparisons because of the flag 'know_modulus'
	friend bool operator == (const single_factor< alg_ideal > & a,
				 const single_factor< alg_ideal > & b);
	friend bool operator != (const single_factor< alg_ideal > & a,
				 const single_factor< alg_ideal > & b);

	friend void gcd(single_factor< alg_ideal > &c,
			const single_factor< alg_ideal > &a,
			const single_factor< alg_ideal > &b);
	friend single_factor< alg_ideal > operator /(const single_factor< alg_ideal > &a,
						     const single_factor< alg_ideal > &b);
	//
	// a specialization of the output routine which produces nicer output
	//
	friend std::ostream & operator << (std::ostream &out,
					   const single_factor< alg_ideal > &f);
	friend std::istream & operator >> (std::istream &in,
					   single_factor< alg_ideal > &f);
};


// friend functions of class single_factor < alg_ideal >

void swap(single_factor< alg_ideal > & a, single_factor< alg_ideal > & b);

    //
    // arithmetic operations
    //

#ifndef HEADBANGER
void multiply(single_factor< alg_ideal > & c,
	      const single_factor< alg_ideal > & a,
	      const single_factor< alg_ideal > & b);

void divide(single_factor< alg_ideal > & c,
	    const single_factor< alg_ideal > & a,
	    const single_factor< alg_ideal > & b);
    // ********************************************************************
    // *WARNING : the software cannot check if 'b' is a divisor of 'a' !!!*
    // ********************************************************************
#endif

lidia_size_t ord_divide(const single_factor< alg_ideal > &a,
			single_factor< alg_ideal > &b);
    
    //
    // factorization algorithms
    //

void factor(factorization< alg_ideal > &, const alg_ideal &,
			   int upper_bound = 34);
factorization< alg_ideal > factor(const alg_ideal &,
				  int upper_bound = 34);

void decompose_prime(prime_ideal * &factor, lidia_size_t & num,
		     const bigint & p, const order & O);
factorization< alg_ideal > finish(const alg_ideal &,
				  rational_factorization &);
void finish(factorization< alg_ideal > &, const alg_ideal &,
	    rational_factorization &);

factorization< alg_ideal > trialdiv(const alg_ideal &,
				    unsigned int upper_bound = 1000000,
				    unsigned int lower_bound = 1);
void trialdiv(factorization< alg_ideal > &, const alg_ideal &,
	      unsigned int upper_bound = 1000000,
	      unsigned int lower_bound = 1);

factorization< alg_ideal > ecm(const alg_ideal &,
			       int upper_bound = 34,
			       int lower_bound = 6, int step = 3);
void ecm(factorization< alg_ideal > &, const alg_ideal &,
	 int upper_bound = 34, int lower_bound = 6, int step = 3);

factorization< alg_ideal > mpqs(const alg_ideal &);
void mpqs(factorization< alg_ideal > &, const alg_ideal &);

    //
    // misc
    //

    //we need new comparisons because of the flag 'know_modulus'
bool operator == (const single_factor< alg_ideal > & a,
		  const single_factor< alg_ideal > & b);
bool operator != (const single_factor< alg_ideal > & a,
		  const single_factor< alg_ideal > & b);

void gcd(single_factor< alg_ideal > &c,
	 const single_factor< alg_ideal > &a,
	 const single_factor< alg_ideal > &b);
single_factor< alg_ideal > operator /(const single_factor< alg_ideal > &a,
				      const single_factor< alg_ideal > &b);
    //
    // a specialization of the output routine which produces nicer output
    //
std::ostream & operator << (std::ostream &out,
			    const single_factor< alg_ideal > &f);
std::istream & operator >> (std::istream &in,
			    single_factor< alg_ideal > &f);



// Name:	factor_p
// Input:	prime (prime number to be decomposed)
//	        &Order (the Order in which the decomposition takes place
//		&number (number of prime ideals occuring in the decomposition
//	        &ideals (ideals[0], ideals[2], ..., ideals[number-1] are these
//		prime ideals)
//
// factor_p is a procedur to decompose an integer prime in a p-maximal
// order of a number field

// Note: The prime ideals are allocated by this procedure, but you have to free
// them yourself!

void factor_p(const bigint & prime, const order & Order,
              lidia_size_t &number, prime_ideal* & ideals);


void gcd(single_factor< alg_ideal > &c,
	 const single_factor< alg_ideal > &a,
	 const single_factor< alg_ideal > &b);


//
// inline functions
//

inline bool
operator != (const single_factor< alg_ideal > & a,
	     const single_factor< alg_ideal > & b)
{
	return !(a == b);
}



inline single_factor< alg_ideal >
operator * (const single_factor< alg_ideal > &a,
	    const single_factor< alg_ideal > &b)
{
	single_factor< alg_ideal > tmp;


	multiply(tmp, a, b);
	return tmp;
}



inline single_factor< alg_ideal >
operator / (const single_factor< alg_ideal > &a,
	    const single_factor< alg_ideal > &b)
{
	single_factor< alg_ideal > tmp;


	divide(tmp, a, b);
	return tmp;
}



inline std::ostream &
operator << (std::ostream &out, const single_factor< alg_ideal > &f)
{
	if (f.prime_flag() == decomposable_object::prime)
		out << f.rep_p;
	else
		out << f.rep_a;
	return out;
}



inline std::istream &
operator >> (std::istream &in, single_factor< alg_ideal > &f)
{
	in >> f.rep_a;
	f.set_prime_flag(decomposable_object::unknown);
	return in;
}



//
// comparisons
//

inline bool
operator< (const single_factor < alg_ideal > & a,
	   const single_factor< alg_ideal > & b)
{
	return false;
}



inline bool
operator > (const single_factor< alg_ideal > & a,
	    const single_factor< alg_ideal > & b)
{
	return false;
}



inline bool
operator <= (const single_factor< alg_ideal > & a,
	     const single_factor< alg_ideal > & b)
{
	return false;
}



inline bool
operator >= (const single_factor< alg_ideal > & a,
	     const single_factor< alg_ideal > & b)
{
	return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SF_ALG_IDEAL_H_GUARD_
