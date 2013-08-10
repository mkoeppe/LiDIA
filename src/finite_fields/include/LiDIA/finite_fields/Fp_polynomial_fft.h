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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FP_POLYNOMIAL_FFT_H_GUARD_
#define LIDIA_FP_POLYNOMIAL_FFT_H_GUARD_


#ifndef LIDIA_FFT_PRIME_H_GUARD_
# include	"LiDIA/fft_prime.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



lidia_size_t next_power_of_two(lidia_size_t m);
lidia_size_t square_root(lidia_size_t n);


//--------------------------------------------------------------
// NOTE : for implementations see file Fp_polynomial/fft_reps.c
//--------------------------------------------------------------


class Fp_polynomial;
class crt_table;
class crt;
class fft_prime;
class poly_mod_rep;

class fft_data
{
	friend class base_fft_rep;
	friend class modular_fft_rep;
	friend class fft_rep;

public:
	struct crt_item
	{
		bigint p;
		crt_table * ct;
		crt_item * next;
	};
	struct item
	{
		fft_prime * prime;
		lidia_size_t num_primes;
		lidia_size_t max_degree;

		crt_item * crt_list;

		int ref_counter;
		item * next;

		~item();
	};

private:
	static item* head;

public:
	item * node; //node != 0  implies
	const crt_item * crt_node; //	crt_node != 0 and
	lidia_size_t k; //	k > 0  (or at least k >= 0)

public:
	fft_data();
	fft_data(const fft_data &);
	fft_data(const bigint &p, lidia_size_t l);
	~fft_data();

	void release(); //decr. ref_counter, may delete object
	void add_ref() const
	{
		if (node)
			node->ref_counter++;
	}

	fft_data & operator = (const fft_data &);
	void init(const bigint &p, lidia_size_t l);

	lidia_size_t number_of_primes() const;

	friend bool operator == (const fft_data& a, const fft_data& b);
};



inline bool operator == (const fft_data& a, const fft_data& b)
{
	return (a.node == b.node) && (a.crt_node == b.crt_node);
}



inline bool operator != (const fft_data& a, const fft_data& b)
{
	return !(a == b);
}



class base_fft_rep
{
protected:
	static fft_prime_t * stat_vec; //vector of size (1 << stat_k)
	static lidia_size_t  stat_k; //stat_k >= max_k for all (mod_)fft_rep

	fft_data fd;
	int k;
	crt *c;

	base_fft_rep(); //disable
	~base_fft_rep();

	void init_(const fft_data &);

	void to_mod_rep(fft_prime_t *x, const Fp_polynomial &f,
			lidia_size_t lo, lidia_size_t hi, lidia_size_t index);
	void combine(fft_prime_t *x, lidia_size_t l, lidia_size_t index);
	void get_res(bigint *x, lidia_size_t l);


	friend void add_mul(modular_fft_rep &x, const modular_fft_rep &a,
			    const modular_fft_rep &b, const modular_fft_rep &c,
			    const modular_fft_rep &d, lidia_size_t index);

};

//    initialization of modular_fft_rep (and fft_rep):
//    ------------------------------------------------
//
//    modular_fft_rep();
//	does nothing
//    modular_fft_rep(const fft_data &F);
//	initializes fd AND k/max_k (calls set_size(F.node->max_degree))
//    init(const fft_data &F)
//   	only initializes fd, sets k to -1
//    set_size(lidia_size_t l)
//	only affects k/max_k;
//	fd must be set before, l must be <= fd.node->max_degree
//
//    It is an error to modify fd after it has been initialized.


class modular_fft_rep : public base_fft_rep
{
private:
	// fft_data fd; ->see base_fft_rep
	// int k; ->see base_fft_rep
	int max_k;
	fft_prime_t *vec;

	modular_fft_rep(const modular_fft_rep &); //disable
public:

	void init(const fft_data &F)
	{
		init_(F);
	}

	modular_fft_rep();
	modular_fft_rep(const fft_data &F) : max_k(-1), vec(0)
	{
		init(F);
		set_size(F.k);
	}
	//an error is raised if F is not initialized

	~modular_fft_rep();

	//const fft_data & data() const { return fd; }

	void set_size(lidia_size_t);

	// computes an n = 2^k point convolution.
	// if deg(f) >= 2^k, then x is first reduced modulo X^n-1.
	void to_modular_fft_rep(const Fp_polynomial &f, lidia_size_t index)
	{
		base_fft_rep::to_mod_rep(vec, f, 0, f.degree(), index);
	}
	void to_modular_fft_rep(const Fp_polynomial &f,
				lidia_size_t lo, lidia_size_t hi, lidia_size_t index)
	{
		base_fft_rep::to_mod_rep(vec, f, lo, hi, index);
	}

	void to_modular_fft_rep(const poly_mod_rep &a, lidia_size_t lo,
				lidia_size_t hi, lidia_size_t index);
	//see "fft_arith.cc"

	void from_modular_fft_rep(lidia_size_t lo, lidia_size_t hi,
				  lidia_size_t index);

	void get_result(Fp_polynomial &x, lidia_size_t, lidia_size_t);
	void get_result_ptr(bigint *x, lidia_size_t, lidia_size_t);

	friend void multiply(modular_fft_rep &x, const modular_fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void add(modular_fft_rep &x, const modular_fft_rep &a,
			const modular_fft_rep &b, lidia_size_t index);
	friend void subtract(modular_fft_rep &x, const modular_fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void add_mul(modular_fft_rep &x, const modular_fft_rep &a,
			    const modular_fft_rep &b, const modular_fft_rep &c,
			    const modular_fft_rep &d, lidia_size_t index);


	friend void multiply(modular_fft_rep &x, const fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void subtract(modular_fft_rep &x, const fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void reduce(fft_rep &xb, const modular_fft_rep &a,
			   lidia_size_t, lidia_size_t index);
	friend void reduce(modular_fft_rep &x, const fft_rep &a,
			   lidia_size_t, lidia_size_t index);

};



class fft_rep : public base_fft_rep
{
private:
	// fft_data fd; ->see base_fft_rep
	// int k; ->see base_fft_rep
	int max_k;
	fft_prime_t **tbl;
	int num_primes;

	fft_rep(const fft_rep &); //disable
public:
	fft_rep();
	fft_rep(const fft_data &F) : max_k(-1), tbl(0), num_primes(0)
	{
		init(F);
		set_size(F.k);
	}
	//an error is raised if F is not initialized

	~fft_rep();

	//const fft_data & data() const { return fd; }

	void init(const fft_data &);
	void set_size(lidia_size_t);

	void to_fft_rep(const Fp_polynomial &f, lidia_size_t lo, lidia_size_t hi);
	void to_fft_rep(const Fp_polynomial &f)
	{
		to_fft_rep(f, 0, f.degree());
	}

	void to_fft_rep(const poly_mod_rep &a, lidia_size_t lo, lidia_size_t hi);
	//see "fft_arith.cc"

	void from_fft_rep(Fp_polynomial &x, lidia_size_t lo, lidia_size_t hi);

	friend void multiply(fft_rep &x, const fft_rep &a, const fft_rep &b);
	friend void subtract(fft_rep &x, const fft_rep &a, const fft_rep &b);
	friend void reduce(fft_rep &x, const fft_rep &a, lidia_size_t k);


	friend void multiply(modular_fft_rep &x, const fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void subtract(modular_fft_rep &x, const fft_rep &a,
			     const modular_fft_rep &b, lidia_size_t index);
	friend void reduce(fft_rep &xb, const modular_fft_rep &a,
			   lidia_size_t, lidia_size_t index);
	friend void reduce(modular_fft_rep &x, const fft_rep &a,
			   lidia_size_t, lidia_size_t index);

	//*************************************************************************
	// three special purpose functions, used only in
	// void update_map(base_vector< bigint > & x, const base_vector< bigint > & a,
	//                 const poly_multiplier& B, const poly_modulus& F)
	//*************************************************************************
	void rev_to_fft_rep(const base_vector< bigint > &x, lidia_size_t lo,
			    lidia_size_t hi, lidia_size_t offset);
	void rev_from_fft_rep(base_vector< bigint > &x, lidia_size_t lo,
			      lidia_size_t hi);
	void add_expand(const fft_rep& a);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FP_POLYNOMIAL_FFT_H_GUARD_
