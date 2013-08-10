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
//	Author	: Thorsten Rottschaefer (TR)
//                Victor Shoup (VS) and Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FFT_PRIME_H_GUARD_
#define LIDIA_FFT_PRIME_H_GUARD_


#ifndef LIDIA_BIT_REVERSE_TABLE_H_GUARD_
# include	"LiDIA/finite_fields/bit_reverse_table.h"
#endif
#ifndef LIDIA_UDIGIT_MOD_H_GUARD_
# include	"LiDIA/udigit_mod.h"
#endif
#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



typedef udigit fft_prime_t;

class fft_prime
{
private:

	//
	// member variables
	//

        static bit_reverse_table rev_table;

	fft_prime_t  p; // prime number that defines the prime field
	fft_prime_t* RootTable; // RootTable[j] = w^{2^{current_degree-j}},
	// where w is a primitive 2^current_degree root of unity
	// in GF(p)
	fft_prime_t* RootInvTable; // RootInvTable[j] = 1/RootTable[j] mod p
	fft_prime_t* TwoInvTable; // TwoInvTable[j] = 1/2^j mod p

	lidia_size_t current_degree; // 2^{current_degree} is the largest possible convolution size, if the current computed tables are used.

        lidia_size_t max_degree; // 2^{max_degree} is the largest possible convolution size for p.

        static fft_prime_t* A;
        static lidia_size_t size_A;

        static fft_prime_t* B;
        static lidia_size_t size_B;

        static fft_prime_t* C;
        static lidia_size_t size_C;

	//
	// private member functions
	//

	bool build_tables (lidia_size_t l);
	// computes RootTable, RootInvTable, and TwoInvTable
	// for a convolution size less or equal to 2^{l}.
	// If l <= current_degree, do nothing.
        // If the pointers are zero, initialize them.
        // Otherwise free them and reinitialize.
	// Set current_degree = l.
	// Return true, if the tables could be computed
        // (i.e, 2^l primitve root of unity exists);
	// Return false, otherwise.

	void FFT(fft_prime_t* A, const fft_prime_t* a, lidia_size_t k, fft_prime_t q, const fft_prime_t* root);
        // Performs a 2^k-point convolution mod q.
	// Assumes, that a has size 2^k
	// and root[i] is a primitve 2^{i}-th root of unity, 1 <= i <= k.

	bool evaluate(fft_prime_t *o, int type, const void* a, lidia_size_t d, lidia_size_t k);
	// Performs a 2^k (point, value) representation mod p.
        // Assumes, that "a" represents a degree d polynomial modulo p.
        //
        // if p = 0, lidia_error_handler is called
	// if k <= current_degree ok = true; else ok = build_tables (k).
	// if (ok) fills "a" with zeroes to get size 2^k, if necessary.
	// if (ok) perform the convolution by calling FFT
	// return ok.


	bool interpolate2(fft_prime_t*, const fft_prime_t*, lidia_size_t,
			  lidia_size_t, lidia_size_t);
   	// special version for Fp_polynomial

	bool interpolate(fft_prime_t* o, const fft_prime_t* a, lidia_size_t k);
        // Performs a 2^k coefficient representation mod p.
	// Assumes, that "a" represents a 2^k (point, value) representation
	// modulo p.
	//
	// if p = 0, lidia_error_handler is called
	// if k <= current_degree ok = true; else ok = build_tables (k).
	// if (ok) fills "a" with zeroes to get size 2^k, if necessary.
	// if (ok) perform the convolution by calling FFT
	// return ok.
	//
        // This function assumes, that the user exactly knows what he is doing;
        // i.e., he must not change the tables of the fft_prime by evaluating
        // another polynomial with a larger current_degree (other root of unity),
        // before calling the interpolate routine for the first polynomial.

public:
	void pointwise_multiply(fft_prime_t* x, const fft_prime_t* a,
				const fft_prime_t* b, lidia_size_t k) const;
   	// x[i] = a[i]*b[i] mod p; i = 0 .. (2^k-1)

	void pointwise_add(fft_prime_t* x, const fft_prime_t* a,
			   const fft_prime_t* b, lidia_size_t k) const;
	// x[i] = a[i]+b[i] mod p; i = 0 .. (2^k-1)

	void pointwise_subtract(fft_prime_t* x, const fft_prime_t* a,
				const fft_prime_t* b, lidia_size_t k) const;
	// x[i] = a[i]-b[i] mod p; i = 0 .. (2^k-1)

public:

	//
	// public functions
	//

	fft_prime ();
        // initializes all pointers with NULL
        // and all other variables with zero

	fft_prime(const fft_prime & a);
   	// copy constructor

	~fft_prime ();
        // deletes all pointers,
        // if they are not equal to NULL

	fft_prime & operator = (const fft_prime & a);
	void assign(const fft_prime & a);
   	// assignment

	void set_prime (fft_prime_t q);
   	// assigns p = q, nothing else

	inline fft_prime_t get_prime() const
	{
		return p;
	}

	lidia_size_t get_max_degree ();
        // returns actual max_degree, if p > 0

	friend bool multiply_fft (void* o, int type,
				  const void* const a,
				  lidia_size_t da,
				  const void* const b,
				  lidia_size_t db,
				  fft_prime & fftprime);
        // a and b represent polynomials of degree da and db respectively,
        // modulo fftprime.p.
        // Degree < 0 represents the zero polynomial.
        //
	// This computes the product polynomial into x using the FFT, i.e.
        // evaluate + pointwise_multiply + interpolate.
        //
        // Adequate space for x (da+db+1) must be allocated by the caller.
        // Squarings are automatically detected and optimized.
        // Output polynomial may alias the input polynomial.
        //
        // if type is 0 o, a and b have to be of type fft_prime_t
        // if type is 1 o, a and b have to be of type udigit_mod

	friend class base_fft_rep;
	friend class modular_fft_rep;
	friend class fft_rep; // these call evaluate & interpolate
	friend class fft_data; // calls build_tables
};

bool multiply_fft (void* o, int type, const void* const a,
		   lidia_size_t da, const void* const b,
		   lidia_size_t db, fft_prime & fftprime);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FFT_PRIME_H_GUARD_
