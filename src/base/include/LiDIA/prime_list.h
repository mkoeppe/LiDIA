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
//	Author	: Anja Jantke (AJ), Patrick Theobald (PT),
//                Dirk Schramm (DS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PRIME_LIST_H_GUARD_
#define LIDIA_PRIME_LIST_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// DATA TYPES
//

#define PRIME_LIST_NUMBER       unsigned long
// used to store primes, may be signed or unsigned

#define PRIME_LIST_COUNTER      long
// used during calculation of primes, may be signed or unsigned
// the greatest prime that can be calculated is about
// (MAX(PRIME_LIST_COUNTER) / 2) ^ 2
// where MAX(TYPE) is the greatest number that TYPE can store

#define PRIME_LIST_FLOAT_NUMBER double
// used during calculation of primes
// the mantissa should have at least as many bits as the greatest
// prime number

#define PRIME_LIST_DIFF         unsigned char
const int prime_list_max_diff = 2 * 255; // maximum difference
// PRIME_LIST_DIFF is used to hold the difference between
// neighbored primes, may be signed or unsigned
// the maximum difference is MAX(PRIME_LIST_DIFF) * 2
// this limits the prime numbers that can be stored in the list
//
// MAX(PRIME_LIST_DIFF) | max. diff. | max. prime         | conjecture*
// ---------------------+------------+--------------------+----------------
//        127           |      254   | > 387000000 > 2^28 | ~ 8.3 * 10^6
//        255 (default) |      510   | > 3 * 10^11 > 2^38 | ~ 6.4 * 10^9
//      32767           |    65534   |                    | ~ 1.5 * 10^111
//      65535           |   131070   |                    | ~ 1.7 * 10^157
//
// * conjecture (Shanks): first separation of width d appears
//   at approximate e ^ sqrt(d)

#define PRIME_LIST_SIEVE        char
// used to store boolean value in standard sieve-arrays

#define PRIME_LIST_BIT_SIEVE    unsigned long
const int prime_list_bit_sieve_size = SIZEOF_LONG * 8; // bit-size of PRIME_LIST_BIT_SIEVE
// used for bit-sieve-arrays
// should be unsigned and have the same size as processor registers


const int prime_list_block_size = 1024;
// number of differences per prime-block
// smaller blocks mean faster random access (get_prime, get_index,
// is_element) but need more memory all in all and vice versa


class prime_list
{

//
// INTERNAL LIST MANAGEMENT
//

private:

	class prime_block
	{
	public:
		PRIME_LIST_NUMBER first_prime;
		PRIME_LIST_DIFF   diff[prime_list_block_size];
		prime_block       *prev_block;
		prime_block       *next_block;
	};

	PRIME_LIST_NUMBER lower_bound;
	PRIME_LIST_NUMBER upper_bound;
	lidia_size_t      number_of_primes;
	PRIME_LIST_NUMBER first_prime;
	PRIME_LIST_NUMBER last_prime;

	lidia_size_t      number_of_blocks;
	prime_block       *first_block;
	prime_block       *last_block;
	prime_block       **block_list;
	lidia_size_t      first_diff_index;
	lidia_size_t      last_diff_index;

	mutable lidia_size_t      current_index;
	mutable PRIME_LIST_NUMBER current_prime;
	mutable prime_block       *current_block;
	mutable lidia_size_t      current_diff_index;


	PRIME_LIST_NUMBER block_first_prime(prime_block *block) const
	{
		return (block == first_block) ? first_prime : block->first_prime;
	}

	PRIME_LIST_NUMBER block_last_prime(prime_block *block) const
	{
		return (block == last_block) ? last_prime : block->next_block->first_prime;
	}

	lidia_size_t block_first_diff_index(prime_block *block) const
	{
		return (block == first_block) ? first_diff_index : 0;
	}

	lidia_size_t block_last_diff_index(prime_block *block) const
	{
		return (block == last_block) ? last_diff_index : prime_list_block_size - 1;
	}

	lidia_size_t block_size(prime_block *block) const
	{
		return block_last_diff_index(block) - block_first_diff_index(block) + 1;
	}

	void init_prime_list();

	void create_prime_list(PRIME_LIST_NUMBER prime);

	void add_next_prime(PRIME_LIST_NUMBER prime);

	void check_and_add_next_prime(PRIME_LIST_NUMBER prime);

	void add_prev_prime(PRIME_LIST_NUMBER prime);

	void check_and_add_prev_prime(PRIME_LIST_NUMBER prime);

	void create_block_list();

	lidia_size_t find_prime(PRIME_LIST_NUMBER prime, bool set_current_prime = false) const;

	void release_prime_list();

//
// CONSTRUCTORS
//

public:

	prime_list();
	prime_list(PRIME_LIST_NUMBER init_upper_bound, char mode = '6');
	prime_list(PRIME_LIST_NUMBER init_lower_bound, PRIME_LIST_NUMBER init_upper_bound, char mode = '6');
	prime_list(const char *filename, lidia_size_t max_number_of_primes = 0);
	prime_list(const prime_list& A);

//
// DESTRUCTOR
//

	~prime_list();

//
// GENERATION
//

private:

	void sieve(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending, char mode);
	void sieve_e(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending);
	void sieve_6k(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending);
	void sieve_ebit(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending);
	void sieve_6kbit(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending);
	void sieve_int(PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending);

//
// BASIC FUNCTIONS
//

public:

	PRIME_LIST_NUMBER get_lower_bound() const
	{
		return lower_bound;
	}

	PRIME_LIST_NUMBER get_upper_bound() const
	{
		return upper_bound;
	}

	void set_lower_bound(PRIME_LIST_NUMBER new_lower_bound, char mode = '6');

	void set_upper_bound(PRIME_LIST_NUMBER new_upper_bound, char mode = '6');

	void resize(PRIME_LIST_NUMBER new_lower_bound, PRIME_LIST_NUMBER new_upper_bound, char mode = '6');

	bool is_empty() const
	{
		return number_of_primes == 0;
	}

	lidia_size_t get_number_of_primes() const
	{
		return number_of_primes;
	}

	bool is_element(PRIME_LIST_NUMBER prime) const;

	lidia_size_t get_index(PRIME_LIST_NUMBER prime) const;

//
// ACCESS FUNCTIONS
//

public:

	PRIME_LIST_NUMBER get_current_prime() const
	{
		return current_prime;
	}

	lidia_size_t get_current_index() const
	{
		return current_index;
	}

	PRIME_LIST_NUMBER get_first_prime() const;

	PRIME_LIST_NUMBER get_first_prime(PRIME_LIST_NUMBER prime) const;

	PRIME_LIST_NUMBER get_next_prime() const;

	PRIME_LIST_NUMBER get_last_prime() const;

	PRIME_LIST_NUMBER get_last_prime(PRIME_LIST_NUMBER prime) const;

	PRIME_LIST_NUMBER get_prev_prime() const;

	PRIME_LIST_NUMBER get_prime(lidia_size_t index) const;

	PRIME_LIST_NUMBER operator [] (lidia_size_t index) const
	{
		return get_prime(index);
	}

//
// INPUT/OUTPUT
//

public:

	void load_from_file(const char *filename, lidia_size_t max_number_of_primes = 0);

	void save_to_file(const char *filename, bool descending = false) const;

//
// ASSIGNMENT
//

public:

	void assign(const prime_list& A);

	prime_list& operator = (const prime_list& A)
	{
		assign(A);
		return *this;
	}

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_PRIME_LIST_H_GUARD_
