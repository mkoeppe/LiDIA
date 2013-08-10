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


#ifndef LIDIA_FP_POLYNOMIAL_UTIL_H_GUARD_
#define LIDIA_FP_POLYNOMIAL_UTIL_H_GUARD_


#ifndef HEADBANGER


#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//***************************************************************
//
//			    tools
//
//***************************************************************




class my_timer
{
	timer t;
	char *msg;
public:
	my_timer();
	~my_timer();

	void start(const char* message);
	void stop();
};





/****************************************************************

    			class poly_matrix

****************************************************************/

class poly_matrix
{
private:
	Fp_polynomial elts[2][2];

	poly_matrix();
	poly_matrix(const poly_matrix&); // disable


	void iter_half_gcd(Fp_polynomial& U, Fp_polynomial& V, lidia_size_t d_red);


	void multiply_and_kill(poly_matrix& B, poly_matrix& C);
	// (*this) = B*C, B and C are destroyed
public:

	poly_matrix(const bigint &p);

	~poly_matrix() { }

	void kill();

	poly_matrix & operator = (const poly_matrix&);

	Fp_polynomial& operator() (lidia_size_t i, lidia_size_t j)
	{
		//NO CHECK
		return elts[i][j];
	}

	const Fp_polynomial& operator() (lidia_size_t i, lidia_size_t j) const
	{
		//NO CHECK
		return elts[i][j];
	}

	void half_gcd(const Fp_polynomial& U, const Fp_polynomial& V, lidia_size_t d_red);
	// deg(U) > deg(V), 1 <= d_red <= deg(U)+1.
	//
	// This computes a 2 x 2 polynomial matrix M_out such that
	//    M_out * (U, V)^T = (U', V')^T,
	// where U', V' are consecutive polynomials in the Euclidean remainder
	// sequence of U, V, and V' is the polynomial of highest degree
	// satisfying deg(V') <= deg(U) - d_red.

	void xhalf_gcd(Fp_polynomial& U, Fp_polynomial& V, lidia_size_t d_red);
	// same as above, except that U is replaced by U', and V by V'

	void multiply(Fp_polynomial& U, Fp_polynomial& V) const;
	// (U, V)^T = (*this) * (U, V)^T

};



//****************************************************************
//
//		    class poly_argument
//
//****************************************************************

// The routine poly_argument::build (see below) which is implicitly called
// by the various compose and update_map routines builds a table
// of polynomials.  If poly_arg_bound > 0, then only up to
// about poly_arg_bound bigints mod p are allocated.  Setting this too
// low may lead to slower execution.
// If poly_arg_bound <= 0, then it is ignored, and space is allocated
// so as to maximize speed.
// Initially, poly_arg_bound = 0.

// If a single h is going to be used with many g's
// then you should build a poly_argument for h,
// and then use the compose routine below.
// poly_argument::build computes and stores h, h^2, ..., h^m mod f.
// After this pre-computation, composing a polynomial of degree
// roughly n with h takes n/m multiplies mod f, plus n^2
// scalar multiplies.
// Thus, increasing m increases the space requirement and the pre-computation
// time, but reduces the composition time.
// If poly_arg_bound > 0, a table of size less than m may be built.


//template < class T > class factorization;

class poly_argument
{
	base_vector< Fp_polynomial > H;

	//PolyargBound can only be changed if no poly_arguments
	//have been declared yet
	static lidia_size_t poly_arg_bound;
	static bool change_enable;

	static
	void compose_InnerProd(Fp_polynomial& x, const Fp_polynomial& v,
			       lidia_size_t low, lidia_size_t high,
			       const base_vector< Fp_polynomial > & H, lidia_size_t n,
			       bigint *t);

public:
	poly_argument()
	{
		change_enable = false;
		debug_handler("poly_argument", "poly_argument (void)");
	}

	~poly_argument()
	{
		debug_handler("poly_argument", "destructor");
	}

	static void set_poly_arg_bound(lidia_size_t n);
	static lidia_size_t get_poly_arg_bound()
	{
		return poly_arg_bound;
	}

	void build(const Fp_polynomial& h, const Fp_poly_modulus& F, lidia_size_t m);
	//m must be > 0, otherwise an error is raised


	void compose(Fp_polynomial& x, const Fp_polynomial& g,
		     const Fp_poly_modulus& F) const;
    	//called by the various compose-routines below


	friend void project_powers(base_vector< bigint > & x,
				   const base_vector< bigint > & a, lidia_size_t k,
				   const poly_argument& H, const Fp_poly_modulus& F);
};



//****************************************************************
//
//		    class int_factor
//		fac_vec " = " vector< int_factor >
//
//****************************************************************

class fac_vec;

class int_factor
{
public:
	lidia_size_t q; // base
	int a; // exponent
	int b;
	double len;

private:
	int_factor() {}; // allow creation only in class fac_vec
	~int_factor() {};
	int_factor & operator = (const int_factor c)
	{
		q = c.q;
		a = c.a;
		len = c.len;
		b = c.b;
		return *this;
	}

	friend class fac_vec;
};



class fac_vec
{
	int_factor *fvec;
	int num_factors;

	fac_vec(); // disable
	fac_vec(const fac_vec &);

	void compute_factorization(lidia_size_t n);
	void sort();

public:
	fac_vec(lidia_size_t n);

	~fac_vec()
	{
		debug_handler("fac_vec", "destructor");
		delete[] fvec;
	}

	void clear(int lo, int hi);
    	//fvec[i].b = 0   for   i = lo, .., hi

	const int_factor& operator[] (int n) const
	{
		debug_handler("fac_vec", "operator [] (int) const");
		return fvec[n];
	}

	int_factor& operator[] (int n)
	{
		debug_handler("fac_vec", "operator [] (int)");
		return fvec[n];
	}

	int number_of_factors() const
	{
		debug_handler("fac_vec", "number_of_factors (void)");
		return num_factors;
	}

	lidia_size_t prod(int lo, int hi) const;
    	//returns product of q^a for i = lo, .., hi
};



//****************************************************************
//
//			    fractions
//
//****************************************************************

// ************   addition, test for equality   ************

void add_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
	      const Fp_polynomial& b, const Fp_polynomial& c, const Fp_polynomial& d,
	      const Fp_poly_modulus& F);
// computes (x, y) = (a*d + b*c mod f, b*d mod f)
// ALIAS RESTRICTION: inputs may not alias outputs

void subtract_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
		   const Fp_polynomial& b, const Fp_polynomial& c, const Fp_polynomial& d,
		   const Fp_poly_modulus& F);
// computes (x, y) = (a*d - b*c mod f, b*d mod f)
// ALIAS RESTRICTION: inputs may not alias outputs


void plain_add_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
		    const Fp_polynomial& b, const Fp_polynomial& c, const Fp_polynomial& d,
		    const Fp_poly_modulus& F);

void plain_subtract_frac(Fp_polynomial& x, Fp_polynomial& y,
			 const Fp_polynomial& a, const Fp_polynomial& b, const Fp_polynomial& c,
			 const Fp_polynomial& d, const Fp_poly_modulus& F);

// same as above, but slower



bool eq_frac(const Fp_polynomial& a, const Fp_polynomial& b,
	     const Fp_polynomial& c, const Fp_polynomial& d, const Fp_poly_modulus& F);
// if (a*d - b*c = 0 mod f) then 1 else 0


bool plain_eq_frac(const Fp_polynomial& a, const Fp_polynomial& b,
		   const Fp_polynomial& c, const Fp_polynomial& d, const Fp_poly_modulus& F);
// same as above, but slower




// ************   class eq_frac_info   ************

// next come data structures and routines for testing if a fixed
// a/b is equal to one of several c/d mod f, using a fast probabilistic test
// With a very small probability, two unequal fractions may be
// claimed to be equal.

class eq_frac_info
{
	base_vector< bigint > L1, L2;
	bigint modulus;

public:

	void build(const Fp_polynomial& a, const Fp_polynomial& b,
		   const Fp_poly_modulus& F);

	bool eq_frac(const Fp_polynomial& c, const Fp_polynomial& d) const;
};



#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FP_POLYNOMIAL_UTIL_H_GUARD_
