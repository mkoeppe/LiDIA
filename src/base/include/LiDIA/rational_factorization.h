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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
#define LIDIA_RATIONAL_FACTORIZATION_H_GUARD_



#ifndef LIDIA_COMPARATOR_H_GUARD_
# include	"LiDIA/comparator.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_ECM_PRIMES_H_GUARD_
# include	"LiDIA/base/ecm_primes.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class ecm_primes;
class rational_factorization;

class rf_single_factor
{
	friend class rational_factorization;
	friend std::istream & operator >> (std::istream &, rational_factorization &);

private:
	bigint single_base;
	unsigned int factor_state;
	int single_exponent;

	static unsigned int state(unsigned int x)
	{
		return (x & 3);
	}

public:
	rf_single_factor() {};
	~rf_single_factor() {};
	rf_single_factor & operator = (const rf_single_factor &);

	friend void swap(rf_single_factor & a, rf_single_factor & b);

	static const unsigned int dont_know;
	static const unsigned int prime;
	static const unsigned int not_prime;

private:

	friend bool operator < (const rf_single_factor & a, const rf_single_factor & b);
	friend bool operator <= (const rf_single_factor & a, const rf_single_factor & b);
	friend bool operator == (const rf_single_factor & a, const rf_single_factor & b);

	friend std::istream & operator >> (std::istream &, rf_single_factor &);
	friend std::ostream & operator << (std::ostream &, const rf_single_factor &);
};

void swap(rf_single_factor & a, rf_single_factor & b);


class rational_factorization
{

private:
	sort_vector< rf_single_factor > factors;
	int isign;
	static int info;
	unsigned int decomp_state;

	void compose();
	void sort()
	{
		factors.sort();
	}

	void trialdiv(lidia_size_t index, unsigned int lower_bound,
		      unsigned int upper_bound, ecm_primes & prim);

	void ecm(lidia_size_t index, int jobs_avail, int job_buffer[30][5],
		 ecm_primes & prim);

	void ecm(lidia_size_t index, int rest_number_of_factors, int jobs_avail,
		 int job_buffer[30][5], ecm_primes & prim);

	static const unsigned int ecm_params[29][3];
	unsigned int ecm_read_max(int stell);
	int ecm_job_planing(int strat[30], int buf[30][5]);

	static const float qs_params[66][7];
	rational_factorization & mpqs_impl(lidia_size_t index,
					   ecm_primes & prim);
	rational_factorization & mpqs(lidia_size_t index, ecm_primes & prim);

	void qs_read_par(unsigned int stellen, double &T, unsigned int &M,
			 unsigned int &groesse_FB, unsigned int &P_ONCE,
			 unsigned int &POLY, unsigned int &P_TOTAL,
			 unsigned int &smallstart);
	int create_FB(unsigned int groesse, const bigint &kN, int **FB,
		      ecm_primes &prim);
	bool qs_build_factors(const bigint & N, const bigint & kN,
			      unsigned int index, int *FB);
	int compute_multiplier(const bigint & N, int bis, ecm_primes &prim);
	static void compose(sort_vector< rf_single_factor > & v);
	static void refine2(sort_vector< rf_single_factor > & v, rf_single_factor & sf,
			    const rf_single_factor   & a, const rf_single_factor & b);

	double zeitqs(unsigned int i, bool t = false);

public:

	rational_factorization();
	rational_factorization(int n, int exp = 1);
	rational_factorization(unsigned int n, int exp = 1);
	rational_factorization(long n, int exp = 1);
	rational_factorization(unsigned long n, int exp = 1);
	rational_factorization(const bigint & n, int exp = 1);
	rational_factorization(const bigrational & n, int exp = 1);
	rational_factorization(const rational_factorization &, int exp = 1);
	~rational_factorization();

	rational_factorization & operator = (const rational_factorization &);

#ifndef HEADBANGER
	void assign (long n, int exp = 1);
	void assign (const bigint & n, int exp = 1);
	void assign (const bigrational & n, int exp = 1);
	void assign (const rational_factorization & f, int exp = 1);
#endif


	void convert(base_vector< bigint > & basis, base_vector< int > & exponent);
	bigint base(lidia_size_t index) const;
	int exponent(lidia_size_t index) const;
	void set_exponent(lidia_size_t index, int expo);

	int sign() const
	{
		return isign;
	}

	void verbose(int a)
	{
		info = a;
	}

	lidia_size_t no_of_comp() const
	{
		return factors.size();
	}

	bool is_prime_factor(lidia_size_t index);
	bool is_prime_factorization();

	//
	// input / output
	//

	void pretty_print(std::ostream &);

	friend std::istream & operator >> (std::istream &, rational_factorization &);
	friend std::ostream & operator << (std::ostream &, const rational_factorization &);

	//
	//  functions for computing with rational_factorizations
	//

	void invert();
	void square();
	friend void multiply (rational_factorization &,
			      const rational_factorization &,
			      const rational_factorization &);
	friend void divide (rational_factorization &,
			    const rational_factorization &,
			    const rational_factorization &);

	void refine();
	bool refine(const bigint &);
	bool refine_comp(lidia_size_t index, const bigint &);

	friend bool operator == (const rational_factorization &,
				 const rational_factorization &);

	//
	// factoring functions
	//


	void trialdiv_comp(lidia_size_t index,
			   unsigned int upper_bound = 1000000,
			   unsigned int lower_bound = 1);
	void trialdiv(unsigned int upper_bound = 1000000,
		      unsigned int lower_bound = 1);

	void ecm_comp(lidia_size_t index, int upper_bound = 34,
		      int lower_bound = 6, int step = 3);
	void ecm(int upper_bound = 34, int lower_bound = 6, int step = 3);

	rational_factorization & mpqs_comp(lidia_size_t index);

	void factor_comp(lidia_size_t index, int upper_bound = 34);
	void factor(int upper_bound = 34);
	void factor_completely()
	{
		while (!is_prime_factorization())
			factor();
	}

	friend rational_factorization trialdiv(const bigint & N,
					       unsigned int upper_bound,
					       unsigned int lower_bound);
	friend rational_factorization ecm(const bigint & N, int upper_bound,
					  int lower_bound, int step);
	friend rational_factorization mpqs(const bigint & N);
	friend rational_factorization factor(const bigint & N);
	friend rational_factorization factor_completely(const bigint & N);

	friend bigrational evaluate(const rational_factorization & r);
	friend bigint evaluate_to_bigint(const rational_factorization & r);

};

std::istream & operator >> (std::istream &, rational_factorization &);
std::ostream & operator << (std::ostream &, const rational_factorization &);

inline
void invert(rational_factorization & f,
	    const rational_factorization & a)
{
    f.assign(a);
    f.invert();
}

inline
void square(rational_factorization & f,
	    const rational_factorization & a)
{
    f.assign(a);
    f.square();
}

inline
rational_factorization operator *(const rational_factorization & a,
				  const rational_factorization & b)
{
    rational_factorization c;
    
    multiply(c, a, b);
    return c;
}

inline
rational_factorization operator / (const rational_factorization & a,
				   const rational_factorization & b)
{
    rational_factorization c;
    
    divide(c, a, b);
    return c;
}

bool operator == (const rational_factorization &,
		  const rational_factorization &);

inline
bool operator != (const rational_factorization & a,
		  const rational_factorization & b)
{
    return !(a == b);
}

rational_factorization trialdiv(const bigint & N,
				unsigned int upper_bound = 1000000,
				unsigned int lower_bound = 1);
rational_factorization ecm(const bigint & N, int upper_bound = 34,
			   int lower_bound = 6, int step = 3);
rational_factorization mpqs(const bigint & N);
rational_factorization factor(const bigint & N);
rational_factorization factor_completely(const bigint & N);

bigrational evaluate(const rational_factorization & r);
bigint evaluate_to_bigint(const rational_factorization & r);




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
