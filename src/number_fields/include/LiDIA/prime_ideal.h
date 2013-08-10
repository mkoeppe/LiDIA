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


#ifndef LIDIA_PRIME_IDEAL_H_GUARD_
#define LIDIA_PRIME_IDEAL_H_GUARD_


#ifndef LIDIA_ALG_NUMBER_H_GUARD_
# include	"LiDIA/alg_number.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// A prime ideal with respect to the base which zbasis and gen2 use.
// A prime ideal always knows its two element presentation.


class prime_ideal
{

	bigint gen1;
	alg_number gen2;
	lidia_size_t e, f; // ramifaction, inertia
	mutable alg_number valu; // A value helpful for computing valuations

	void compute_valu() const;

public:
	prime_ideal();
	prime_ideal(const bigint&, const alg_number&,
		    lidia_size_t ram = 0, lidia_size_t inert = 0);
	prime_ideal(const bigint&, const nf_base * O = nf_base::current_base);
	~prime_ideal() {};

	// Convert to ideal
	operator alg_ideal() const
	{
		return alg_ideal(gen1, gen2);
	}

	// Access Functions
	const bigint & first_generator() const
	{
		return gen1;
	}
	const alg_number & second_generator() const
	{
		return gen2;
	}
	lidia_size_t ramification_index() const
	{
		return e;
	}
	lidia_size_t degree_of_inertia() const
	{
		return f;
	}
	const bigint & base_prime() const
	{
		return gen1;
	}

	// swap
	friend void swap(prime_ideal &a, prime_ideal &b);

	// Higher Level Functions
	friend bigint norm(const prime_ideal & p);
	friend bigint exponent(const prime_ideal & p);
	friend long ord(const prime_ideal &, const alg_ideal &); // ord_p(m);
	friend long ord(const prime_ideal &, const alg_number &);
	// ord of principal ideal
	friend void power(alg_ideal &, const prime_ideal &, const bigint &);

	// Comparision -- needed since cast is ambiguous
	bool operator == (const prime_ideal &p) const
	{
		return (gen1 == p.gen1 && gen2 == p.gen2);
	}

	bool operator != (const prime_ideal &p) const
	{
		return (gen1 != p.gen1 || gen2 != p.gen2);
	}

	// NOT YET IMPLEMENTED FUNCTIONS:
	bool is_principal() const;
	alg_number generator() const;
	// return generator if ideal is principal else error

	// In-/Output:
	friend std::ostream& operator << (std::ostream &, const prime_ideal &);
	friend std::istream& operator >> (std::istream &, prime_ideal &);
};



// Higher Level Functions
inline bigint norm(const prime_ideal & p)
{
	bigint c;

	power(c, p.gen1, p.f);
	return c;
} // = gen^f



inline bigint exponent(const prime_ideal & p)
{
	return p.gen1;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_PRIME_IDEAL_H_GUARD_
