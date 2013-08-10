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
//	Author	: Frank J. Lehmann (FL), Thomas Pfahler (TPF)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_CRT_TABLE_H_GUARD_
#define LIDIA_CRT_TABLE_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// FIXME: not portable
#define SDIGIT_MAX 0x7fffffff



//***************************************************************
//
//				  class crt_table
//
//***************************************************************

class crt_table
{
	friend class crt;

public :

	static const char MODULAR;
	static const char INTEGER;


private :

	base_vector< sdigit > primes;
	lidia_size_t num_primes;

	char mode;
	bigint p;

	lidia_size_t reference_counter;

	base_vector< bigint > coeff_mod_p;
	base_vector< double > x;

	bigint M, M_over_2;
	bigint minus_M_mod_p; // = - (product of all primes) mod p

	sdigit (*NEXT_PRIME) (const sdigit & q);

	static sdigit next_prime (const sdigit & udq);

	void gen_primes (bigint  B , const sdigit & bound);
	void set_primes (const base_vector< sdigit > & Prime_vec);

	void build_tables ();

private :

	void SINGLE_REMAINDER(sdigit &s, const bigint &b, lidia_size_t index)
	{
		remainder(s, b, primes[index]);
	}

	void SINGLE_REMAINDER(sdigit &s, const bigmod &b, lidia_size_t index)
	{
		remainder(s, b.mantissa(), primes[index]);
	}

	void MULTI_REMAINDER(sdigit *s, const bigint *b, lidia_size_t bsize, lidia_size_t index)
	{
		for (lidia_size_t i = 0; i < bsize; i++) {
			remainder(s[i], b[i], primes[index]);
		}
	}

	void MULTI_REMAINDER(sdigit *s, const bigmod *b, lidia_size_t bsize, lidia_size_t index)
	{
		for (lidia_size_t i = 0; i < bsize; i++) {
			remainder(s[i], b[i].mantissa(), primes[index]);
		}
	}

	void NORMALIZE (bigint & big , bigint & p_over_2 , const bigint & pr)
	{
		if (big.is_positive()) {
			if (big >= p_over_2)
				subtract (big, big, pr);
		}
		else {
			p_over_2.negate();

			if (big <= p_over_2)
				add (big, big, pr);

			p_over_2.negate();
		}
	}

public:

	//
	//  *****  initialization function  *****
	//

	void init (const bigint &B, const bigint & P);
	void init (const bigint &B, const sdigit & prime_bound, const bigint & P);
	void init (const base_vector< sdigit > & Primes, const bigint & P);

	void init (const bigint &B)
	{
		init (B, bigint(0));
	}
	void init (const bigint &B, const sdigit & prime_bound)
	{
		init (B, prime_bound, bigint(0));
	}
	void init (const base_vector< sdigit > & Primes)
	{
		init (Primes, bigint(0));
	}


	void clear ();

	//
	// *****  constructors for INTEGER and MODULAR mode respectively  *****
	//

	crt_table ();
	crt_table (const bigint &B, const bigint & P);
	crt_table (const bigint &B, const sdigit & prime_bound, const bigint & P);
	crt_table (const base_vector< sdigit > & Primes, const bigint & P);

	crt_table (const bigint &B);
	crt_table (const bigint &B, const sdigit & prime_bound);
	crt_table (const base_vector< sdigit > & Primes);


	~crt_table();

	//
	// *****  some tools  *****
	//

	int number_of_primes() const
	{
		return num_primes;
	}

	sdigit get_prime (lidia_size_t index) const;

	lidia_size_t  how_much_to_use    (bigint B) const;
	void set_prime_generator (sdigit (*np) (const sdigit & q) = NULL);

	void info () const;

};
// class crt_table



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_CRT_TABLE_H_GUARD_
