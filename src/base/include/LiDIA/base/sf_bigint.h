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
//	$Id: sf_bigint.h,v 2.12 2006/03/06 12:08:36 lidiaadm Exp $
//
//	Author	: Emre Binisik (EB),  Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SF_BIGINT_H_GUARD_
#define LIDIA_SF_BIGINT_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BASE_FACTOR_H_GUARD_
# include	"LiDIA/base/base_factor.h"
#endif
#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif

#include        <vector>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class base_factor;
template< class T > class single_factor;
template< class T > class factorization;
class ecm_primes;


// the following information is used for storing already tested
// algorithms (we expect algorithms use size sqrt(dec_length))
// size of ECM is stored from Bit 16 - 24

enum fact_info {
	TD = 1,
	PM1 = 2,
	PP1 = 4
};



class bigmod;


typedef single_factor< bigint > sf_bigint;


template<>
class single_factor< bigint > : public base_factor< bigint >
{
	friend class factorization< bigint >;

private:
	static bool info;
	unsigned long algo_info;


public:

	static bool verbose()
	{
		return info;
	}
	static void set_verbose_flag(int i)
	{
		info = (i != 0);
	}


	single_factor(); //default value must be '1'
	single_factor(const single_factor< bigint > &);
	single_factor(const bigint&);
	~single_factor();


	void swap(single_factor< bigint > &b);

	single_factor< bigint > & operator = (const single_factor< bigint > & x);
	single_factor< bigint > & operator = (const bigint & x);
	void assign(const single_factor< bigint > & x);
	void assign(const bigint & x);

	bool is_one() const
	{
		return rep.is_one();
	}

	bool is_prime_factor() const
	{
		return (prime_flag() == prime);
	}

	bool is_prime_factor(int test);
//test == 0 ->no explicit primality test (only a flag is checked)
//test != 0 ->explicit prime-test if prime_state() == unknown

	bigint extract_unit();

	friend lidia_size_t ord_divide(const single_factor< bigint > & a,
				       const single_factor< bigint > & b);

	friend void gcd (single_factor< bigint > & r, const single_factor< bigint > &a,
			 const single_factor< bigint > & b);

public:
	factorization< bigint > TrialDiv(unsigned int upper_bound = 1000000,
					 unsigned int lower_bound = 1);
	friend factorization< bigint > TrialDiv(const bigint & N,
						const unsigned int upper_bound = 1000000,
						const unsigned int lower_bound = 1);

	factorization< bigint > PollardRho(int size = 7);
	friend factorization< bigint > PollardRho(const bigint& x, int size = 7);

	factorization< bigint > PollardPminus1(int size = 9);
	friend factorization< bigint > PollardPminus1(const bigint& x, int size = 9);

	factorization< bigint > WilliamsPplus1(int size = 9);
	friend factorization< bigint > WilliamsPplus1(const bigint& x, int size = 9);

	factorization< bigint > Fermat();
	friend factorization< bigint > Fermat(const bigint& x);

	factorization< bigint > ECM(int upper_bound = 34, int lower_bound = 6,
				    int step = 3, bool jump_to_QS = false);
	friend factorization< bigint > ECM(const bigint & x,
					   int upper_bound = 34,
					   int lower_bound = 6,
					   int step = 3);

	factorization< bigint > MPQS();
	friend factorization< bigint > MPQS(const bigint & x);

	factorization< bigint > factor(int size = 34) const;
	friend factorization< bigint > sf_factor(const bigint & x,
						 int size = 34);
	factorization< bigint > completely_factor() const;
	friend factorization< bigint > completely_factor(const bigint & x);


private:

	void TrialDiv(factorization< bigint > & f,
		      unsigned int lower_bound,
		      unsigned int upper_bound,
		      ecm_primes& prim);

	void PollardRho(factorization< bigint > & f,
			unsigned int B);

	void Fermat(factorization< bigint > & f,
		    unsigned int B);

	void PollardPminus1(factorization< bigint > & f,
			    unsigned int B1,
			    unsigned int B3,
			    const bigint & a,
			    ecm_primes& prim);

	bool PollardPminus1_VSC(factorization< bigint > & f,
				const bigmod & Q,
				unsigned int B1,
				unsigned int B3,
				ecm_primes& prim);

	bool WilliamsPplus1_VSC(factorization< bigint > & f,
				const bigmod & Q,
				unsigned int B1,
				unsigned int B3,
				ecm_primes& prim);

	void WilliamsPplus1(factorization< bigint > & f,
			    unsigned int B1,
			    unsigned int B3,
			    const bigint & a,
			    ecm_primes& prim);

private:        // private members need for ECM

	static const unsigned int ecm_params[29][3];
	unsigned int ecm_read_max(int stell) const;
	int ecm_job_planing(int strat[30], int buf[30][5]);

	void ECM(factorization< bigint > & tmp,
		 int jobs_avail,
		 int job_buffer[30][5],
		 ecm_primes & prim, bool jump_to_QS = false);

private: // private members for MPQS

        struct qs_param_record {
          unsigned int digits;
          double T;
          unsigned int M;
          unsigned int size_FB;
          unsigned int P_ONCE;
          unsigned int P_TOTAL;
          unsigned int smallstart;
        };

	static qs_param_record const qs_params[66];
	qs_param_record const& qs_read_par(int stellen);

	void MPQS(factorization< bigint > & tmp,
		  ecm_primes & prim);

	double zeitqs(unsigned int i, bool t = false);
	bool check_gcd(factorization< bigint > &, const bigmod &, const bigint&);
};



inline void
swap(single_factor< bigint > & a, single_factor< bigint > & b)
{
	a.swap(b);
}


// Added by G.A. - friend functions should be declared also outside the class
factorization< bigint > TrialDiv(const bigint & N,
				 const unsigned int upper_bound,
				 const unsigned int lower_bound);

factorization< bigint > PollardRho(const bigint& x, int size);

factorization< bigint > PollardPminus1(const bigint& x, int size);

factorization< bigint > WilliamsPplus1(const bigint& x, int size);

factorization< bigint > Fermat(const bigint& x);

factorization< bigint > ECM(const bigint & x,
			    int upper_bound,
			    int lower_bound,
			    int step);

factorization< bigint > MPQS(const bigint & x);

factorization< bigint > sf_factor(const bigint & x,
				  int size);

factorization< bigint > completely_factor(const bigint & x);

factorization< bigint > PollardRho(const bigint &, int);

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SF_BIGINT_H_GUARD_
