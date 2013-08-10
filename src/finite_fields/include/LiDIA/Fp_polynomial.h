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


#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
#define LIDIA_FP_POLYNOMIAL_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
# ifndef LIDIA_COMPARATOR_H_GUARD_
# include	"LiDIA/comparator.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGINT_POLYNOMIAL_H_GUARD_
# include	"LiDIA/bigint_polynomial.h"
#endif
#ifndef LIDIA_RESIDUE_CLASS_LIST_H_GUARD_
# include	"LiDIA/base/residue_class_list.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// forward declarations required later on

class Fp_polynomial;
class Fp_pol_ref;
class base_fft_rep;
class fft_rep;
class modular_fft_rep;
class Fp_poly_modulus;
class Fp_poly_multiplier;
class poly_matrix;
class poly_mod_rep;
class eq_frac_info;
class poly_argument;



#ifndef HEADBANGER


inline void
Remainder (bigint & x, const bigint & a, const bigint & p)
{
	remainder(x, a, p);
	if (x.is_lt_zero()) {
		add(x, x, p);
	}
}



inline void
AddMod (bigint & c, const bigint & a, const bigint & b, const bigint & p)
{
	add(c, a, b);
	if (c >= p) {
		subtract(c, c, p);
	}
}



inline void
SubMod (bigint & c, const bigint & a, const bigint & b, const bigint & p)
{
	subtract(c, a, b);
	if (c.is_lt_zero()) {
		add(c, c, p);
	}
}



inline void
MulMod (bigint & c, const bigint & a, const bigint & b, const bigint & p)
{
	multiply(c, a, b);
	Remainder(c, c, p);
}



inline void
NegateMod (bigint & x, const bigint & a, const bigint & p)
{
	if (a.is_zero()) {
		x.assign_zero();
	}
	else {
		subtract(x, p, a);
	}
}



inline void
InvMod(bigint & x, const bigint & a, const bigint & p)
{
	// x = a^{-1} mod n, 0 <= x < p
	// error is raised if inverse not defined
	if (xgcd_left(x, a, p).is_one()) {
		Remainder(x, x, p);
		return;
	}
	else {
		lidia_error_handler("Fp_polynomial", "InvMod(...)::inverse does not exist");
	}
}

#endif 	// HEADBANGER




//************************************************************
//	    class Fp_pol_ref
//    tries to differentiate between the use of
//    Fp_polynomial::operator[] as an l-value and
//    as an r-value
//************************************************************

class Fp_pol_ref
{
	friend class Fp_polynomial; // MM


private:

	class Fp_polynomial &p;
	lidia_size_t	 ix;

	Fp_pol_ref(); // not needed
	Fp_pol_ref (const Fp_pol_ref & a) : p(a.p), ix(a.ix)
	{ }


public:

	Fp_pol_ref(Fp_polynomial &f, lidia_size_t i) : p(f), ix(i)
	{ }

	//l-values
	Fp_pol_ref & operator = (const bigint &a);
	void assign_zero ();
	void assign_one ();
	void assign (const bigint &a);

	//r-value
	operator bigint ();
};



#ifndef HEADBANGER

//************************************************************
//		class crossover_class
//    class for computing crossover-points with respect to
//    bitlengths of used moduli
//************************************************************

#define CROV_NUM_VALUES 10

class crossover_class
{
	int	x[CROV_NUM_VALUES];
	int	y_fftmul[CROV_NUM_VALUES];
	int	y_fftdiv[CROV_NUM_VALUES];
	int	y_fftrem[CROV_NUM_VALUES];
	int	y_inv[CROV_NUM_VALUES];
	int	y_gcd[CROV_NUM_VALUES];
	int	y_xgcd[CROV_NUM_VALUES];
	double	y2_fftmul[CROV_NUM_VALUES];
	double	y2_fftdiv[CROV_NUM_VALUES];
	double	y2_fftrem[CROV_NUM_VALUES];
	double	y2_inv[CROV_NUM_VALUES];
	double	y2_gcd[CROV_NUM_VALUES];
	double	y2_xgcd[CROV_NUM_VALUES];
	int	halfgcd;
	int	log2_newton;

	crossover_class (const crossover_class &); //disable


public :

	crossover_class ();
	void init(const int x_val[CROV_NUM_VALUES],
		  const int fftmul_val[CROV_NUM_VALUES],
		  const int fftdiv_val[CROV_NUM_VALUES],
		  const int fftrem_val[CROV_NUM_VALUES],
		  const int inv_val[CROV_NUM_VALUES],
		  const int gcd_val[CROV_NUM_VALUES],
		  const int xgcd_val[CROV_NUM_VALUES]);

	int fftmul_crossover (const bigint &modulus) const;
	int fftdiv_crossover (const bigint &modulus) const;
	int fftrem_crossover (const bigint &modulus) const;
	int inv_crossover (const bigint &modulus) const;
	int halfgcd_crossover (const bigint &modulus) const;
	int gcd_crossover (const bigint &modulus) const;
	int xgcd_crossover (const bigint &modulus) const;
	int log2_newton_crossover (const bigint &modulus) const;
};

#endif 	// HEADBANGER



//************************************************************
//
//			 Fp_polynomial
//
// The class Fp_polynomial implements polynomial arithmetic modulo p.
// Polynomials are represented as vectors of bigint with values in the
// range 0..p-1.
//
// If f is a Fp_polynomial, then f.coeff is a bigint*.
// The zero polynomial is represented as a zero length vector.
// Otherwise, f.coeff[0] is the constant-term, and f.coeff[f.c_length-1]
// is the leading coefficient, which is always non-zero.
// Use the member function remove_leading_zeros() to strip leading zeros.
//
//**************************************************************



class Fp_polynomial
{
#ifndef HEADBANGER
private:
	static bigint	ZERO; //coeff. value for zero polynomial
	static residue_class_list< bigint > *L; //list of moduli

	bigint		*coeff; //coefficient vector
	lidia_size_t	c_length;
	lidia_size_t	c_alloc;
	residue_class< bigint > *Mp; //pointer to modulus
#endif

	void set_degree (lidia_size_t n)
	{
		debug_handler_l("Fp_polynomial", "set_degree (lidia_size_t)", 2);
		set_max_degree(n);
		c_length = comparator< lidia_size_t >::max(n+1, 0);
	}

	void read_x_term (std::istream &datei, bigint &coeff, lidia_size_t &exp);


public:

	static crossover_class	crossovers; // crossover-points for fft-arithm.
	// depending on architecture


	//***************************************************************
	//
	//	Constructors, Destructors, and Basic Functions
	//
	//***************************************************************

	Fp_polynomial () :
		coeff(0),
		c_length(0),
		c_alloc(0),
		Mp(0)
	{
		debug_handler_l("Fp_polynomial", "Fp_polynomial()", 1);
	}

	Fp_polynomial (const Fp_polynomial& a);
	Fp_polynomial (const polynomial< bigint > & f, const bigint & p);

	~Fp_polynomial();

	void kill();

	void set_modulus (const Fp_polynomial &f);
	//sets modulus to f.modulus(), assigns zero polynomial
	void set_modulus (const bigint &p);
	//sets modulus to p, assigns zero polynomial

	const bigint & modulus() const;

	lidia_size_t degree () const
		// note that the zero polynomial has degree -1.
	{
		debug_handler_l("Fp_polynomial", "degree(void)", 2);
		return c_length - 1;
	}

#ifndef HEADBANGER
	void comp_modulus (const Fp_polynomial& x, const char* fctname) const;
	void set_max_degree (lidia_size_t n);
#endif 	// HEADBANGER

	//***************************************************************
	//
	//	routines for manipulating coefficients
	//
	//***************************************************************

#ifndef HEADBANGER
	void get_coefficient (bigint& a, lidia_size_t i) const;
	void set_coefficient (const bigint& a, lidia_size_t i);
	void set_coefficient (lidia_size_t i);
#endif 	// HEADBANGER

	const bigint& lead_coeff () const;
	const bigint& const_term () const;

	void remove_leading_zeros ();
	void make_monic ();

	friend class Fp_pol_ref;
	Fp_pol_ref operator[] (lidia_size_t i)
	{
		return Fp_pol_ref(*this, i);
	}
	// operator[] with context sensitivity
	// tries to differentiate between use as l-value and r-value

	const bigint & operator[] (lidia_size_t i) const;
	//always r-value


	//***************************************************************
	//
	//	assignments
	//
	//****************************************************************/

	Fp_polynomial & operator = (const Fp_polynomial& a);
	Fp_polynomial & operator = (const bigint& c);

#ifndef HEADBANGER
	void assign (const bigint & a);
	//void assign (const bigmod& a);
	//void assign (const base_vector < bigmod > & a);
	void assign (const base_vector< bigint > & a, const bigint &p);
	void assign (const Fp_polynomial& a);
#endif	// HEADBANGER

	void assign_zero ()
	{
		debug_handler_l ("Fp_polynomial", "assign_zero (void)", 2);
		if (modulus().is_zero()) {
			lidia_error_handler("Fp_polynomial", "assign_zero(void)::modulus = 0");
		}

		c_length = 0; //don't delete the coeff. vector !
	}

	void assign_one ()
	{
		debug_handler_l ("Fp_polynomial", "assign_one (void)", 2);
		if (modulus().is_zero()) {
			lidia_error_handler("Fp_polynomial", "assign_one(void)::modulus = 0");
		}

		set_degree(0);
		coeff[0].assign_one();
	}

	void assign_x ();

	void assign_zero (const bigint &p); //sets modulus to p
	void assign_one (const bigint &p); //	"
	void assign_x (const bigint &p); //	"

	void randomize (lidia_size_t n); //random polynomial of degree n


	//***************************************************************
	//
	//	comparisons
	//
	//***************************************************************

	friend bool operator == (const Fp_polynomial& a, const Fp_polynomial& b);
	friend bool operator != (const Fp_polynomial& a, const Fp_polynomial& b);

	bool is_zero () const;
	bool is_one () const;
	bool is_x () const;

	bool is_monic () const
	{
		return (lead_coeff().is_one());
	}

	bool is_binomial() const;


	//***************************************************************
	//
	//	operators
	//
	//***************************************************************

	Fp_polynomial & operator += (const Fp_polynomial &a);
	Fp_polynomial & operator += (const bigint &a);
	Fp_polynomial & operator -= (const Fp_polynomial &a);
	Fp_polynomial & operator -= (const bigint &a);
	Fp_polynomial & operator *= (const Fp_polynomial &a);
	Fp_polynomial & operator *= (const bigint &a);
	Fp_polynomial & operator /= (const Fp_polynomial &a);
	Fp_polynomial & operator /= (const bigint &a);
	Fp_polynomial & operator %= (const Fp_polynomial &a);

	bigint operator () (const bigint &value) const;


	//***************************************************************
	//
	//	procedural versions for arithmetics
	//
	//***************************************************************


#ifndef HEADBANGER
	void negate ();
	friend void negate (Fp_polynomial& x, const Fp_polynomial& a);
	friend void add (Fp_polynomial& x, const Fp_polynomial& a,
			 const Fp_polynomial& b);
	friend void add (Fp_polynomial & x, const Fp_polynomial& a,
			 const bigint& b);
        friend void add (Fp_polynomial& x, const bigint& a,
			 const Fp_polynomial& b);

	friend void subtract (Fp_polynomial& x, const Fp_polynomial& a,
			      const Fp_polynomial& b);
	friend void subtract (Fp_polynomial & x, const Fp_polynomial& a,
			      const bigint& b);
	friend void subtract (Fp_polynomial& x, const bigint& a,
			      const Fp_polynomial& b);

	friend void multiply (Fp_polynomial& x, const Fp_polynomial& a,
			      const Fp_polynomial& b);
	friend void multiply (Fp_polynomial & x, const Fp_polynomial& a,
			      const bigint& b);
	friend void multiply (Fp_polynomial & x, const bigint& a,
			      const Fp_polynomial& b);
	friend void multiply_by_scalar (Fp_polynomial &c,
					const Fp_polynomial &a,
					const bigint &b);

	// These always use "classical" arithmetic
	friend void plain_mul (Fp_polynomial& x, const Fp_polynomial& a,
			       const Fp_polynomial& b);
	friend void plain_sqr (Fp_polynomial& x, const Fp_polynomial& a);

	// These always use FFT arithmetic
	friend void fft_mul (Fp_polynomial& x, const Fp_polynomial& a,
			     const Fp_polynomial& b);
	friend void fft_sqr (Fp_polynomial& x, const Fp_polynomial& a);
#endif 	// HEADBANGER

	friend void square(Fp_polynomial& x, const Fp_polynomial& a);



#ifndef HEADBANGER
	friend void div_rem (Fp_polynomial& q, Fp_polynomial& r,
			     const Fp_polynomial& a, const Fp_polynomial& b);
	// q = a/b, r = a%b

	friend void divide (Fp_polynomial& q, const Fp_polynomial& a,
			    const Fp_polynomial& b);
	friend void divide (Fp_polynomial& q, const bigint& a,
			    const Fp_polynomial& b);
	friend void divide (Fp_polynomial& q, const Fp_polynomial& a,
			    const bigint& b);
	// q = a/b

	friend void remainder (Fp_polynomial& r, const Fp_polynomial& a,
			       const Fp_polynomial& b);
	// r = a%b
#endif 	// HEADBANGER

	friend void invert (Fp_polynomial& c, const Fp_polynomial& a,
			    lidia_size_t m);
	// computes c = a^{-1} % x^m
	// constant term must be non-zero


#ifndef HEADBANGER
	// These always use "classical" arithmetic
	friend void plain_div_rem (Fp_polynomial& q, Fp_polynomial& r,
				   const Fp_polynomial& a,
				   const Fp_polynomial& b);
	friend void plain_div (Fp_polynomial& q, const Fp_polynomial& a,
			       const Fp_polynomial& b);
	friend void plain_rem (Fp_polynomial& r, const Fp_polynomial& a,
			       const Fp_polynomial& b);

	// These always use FFT arithmetic
	friend void fft_div_rem (Fp_polynomial& q, Fp_polynomial& r,
				 const Fp_polynomial& a,
				 const Fp_polynomial& b);
	friend void fft_div (Fp_polynomial& q, const Fp_polynomial& a,
			     const Fp_polynomial& b);
	friend void fft_rem (Fp_polynomial& r, const Fp_polynomial& a,
			     const Fp_polynomial& b);

	friend void plain_inv (Fp_polynomial& x, const Fp_polynomial& a,
			       lidia_size_t m);
	// always uses "classical" algorithm
	// ALIAS RESTRICTION: input may not alias output

	friend void newton_inv (Fp_polynomial& x, const Fp_polynomial& a,
				lidia_size_t m);
	// uses a Newton Iteration with the FFT.
	// ALIAS RESTRICTION: input may not alias output
#endif	// HEADBANGER



	//***************************************************************
	//
	//	Miscellaneaous
	//
	//****************************************************************/

	friend void swap (Fp_polynomial& x, Fp_polynomial& y);
	// swap x & y (only pointers are swapped)

	friend void trunc (Fp_polynomial& c, const Fp_polynomial& a,
			   lidia_size_t m);
	// c = a % x^m

	friend void shift_right (Fp_polynomial& c, const Fp_polynomial& a,
				 lidia_size_t n);
	// c = a/x^n

	friend void shift_left (Fp_polynomial& c, const Fp_polynomial& a,
				lidia_size_t n);
	// c = a*x^n

	friend void derivative (Fp_polynomial& c, const Fp_polynomial& a);
	// c = derivative of a
        friend Fp_polynomial derivative (const Fp_polynomial& a);

	friend void add_multiple (Fp_polynomial &f,
				  const Fp_polynomial &g,
				  const bigint &s,
				  lidia_size_t n,
				  const Fp_polynomial &h);
	// f = g + s*x^n*h, n >= 0

	friend void cyclic_reduce (Fp_polynomial& x, const Fp_polynomial& a,
				   lidia_size_t m);
	// computes c = a mod x^m-1

	friend void copy_reverse (Fp_polynomial& x, const Fp_polynomial& a,
				  lidia_size_t lo, lidia_size_t hi);
	// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
	// input may not alias output

	void build_from_roots (const base_vector< bigint > & a);
	// computes the polynomial (x-a[0]) ... (x-a[n-1]),
        // where n = a.length()
	// implementation: see file "fftbuild.c"

	friend void multiply_by_x_mod (Fp_polynomial& c,
				       const Fp_polynomial& a,
				       const Fp_polynomial& f);
	// c = (a * x) mod f


	//********************************************************************
	//
	// I/O
	//
	// two I/O formats:
	//
	// 1)   [a_n ... a_1 a_0] mod p
	// represents the polynomial a_n*x^n + ... + a_1*x + a_0 mod p.
	//
	// 2)   a_n*x^n + ... + a_1*x + a_0 mod p
	// is self-explanatory.
	// The default is 2).
	// 
	// On output, all coefficients will be integers between 0 and p-1,
	// and a_n not zero (the zero polynomial is [ ]).
	// On input, the coefficients are arbitrary integers which are
	// then reduced modulo p, and leading zeros stripped. The function
	// Fp_polynomial::read(std::istream&) can handle both formats; the decision
	// is made according to the first character being '[' or not.
	//
	//*********************************************************************

	void read (std::istream & s);
	void write (std::ostream & s) const; //writes format 1)
	void pretty_read (std::istream & s); //reads format 2)
	void pretty_print (std::ostream & s = std::cout) const; //writes format 2)


	//***************************************************************
	//
	//	special functions that have to be friends
	//
	//***************************************************************

	friend class base_fft_rep;
	friend class fft_rep;
	friend class modular_fft_rep;
	friend class Fp_poly_modulus;
	friend class Fp_poly_multiplier;
	friend class poly_matrix;
	friend class poly_mod_rep;
	friend class eq_frac_info;
	friend class poly_argument;


	friend void remainder (Fp_polynomial &, const Fp_polynomial &, const Fp_poly_modulus &);
	//set_degree


}; //end class Fp_polynomial

// declaration of Fp_polynomial's friends

bool operator == (const Fp_polynomial& a, const Fp_polynomial& b);
bool operator != (const Fp_polynomial& a, const Fp_polynomial& b);

void negate (Fp_polynomial& x, const Fp_polynomial& a);
void add (Fp_polynomial& x, const Fp_polynomial& a,
	  const Fp_polynomial& b);
void add (Fp_polynomial & x, const Fp_polynomial& a,
	  const bigint& b);
void add (Fp_polynomial& x, const bigint& a,
	  const Fp_polynomial& b);

void subtract (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_polynomial& b);
void subtract (Fp_polynomial & x, const Fp_polynomial& a,
	       const bigint& b);
void subtract (Fp_polynomial& x, const bigint& a,
	       const Fp_polynomial& b);

void multiply (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_polynomial& b);
void multiply (Fp_polynomial & x, const Fp_polynomial& a,
	       const bigint& b);
void multiply (Fp_polynomial & x, const bigint& a,
	       const Fp_polynomial& b);
void multiply_by_scalar (Fp_polynomial &c,
			 const Fp_polynomial &a,
			 const bigint &b);

	// These always use "classical" arithmetic
void plain_mul (Fp_polynomial& x, const Fp_polynomial& a,
		const Fp_polynomial& b);
void plain_sqr (Fp_polynomial& x, const Fp_polynomial& a);

	// These always use FFT arithmetic
void fft_mul (Fp_polynomial& x, const Fp_polynomial& a,
	      const Fp_polynomial& b);
void fft_sqr (Fp_polynomial& x, const Fp_polynomial& a);

void square(Fp_polynomial& x, const Fp_polynomial& a);



#ifndef HEADBANGER
void div_rem (Fp_polynomial& q, Fp_polynomial& r,
	      const Fp_polynomial& a, const Fp_polynomial& b);
	// q = a/b, r = a%b

void divide (Fp_polynomial& q, const Fp_polynomial& a,
	     const Fp_polynomial& b);
void divide (Fp_polynomial& q, const bigint& a,
	     const Fp_polynomial& b);
void divide (Fp_polynomial& q, const Fp_polynomial& a,
	     const bigint& b);
	// q = a/b

void remainder (Fp_polynomial& r, const Fp_polynomial& a,
		const Fp_polynomial& b);
	// r = a%b
#endif 	// HEADBANGER

void invert (Fp_polynomial& c, const Fp_polynomial& a,
	     lidia_size_t m);
	// computes c = a^{-1} % x^m
	// constant term must be non-zero


#ifndef HEADBANGER
	// These always use "classical" arithmetic
void plain_div_rem (Fp_polynomial& q, Fp_polynomial& r,
		    const Fp_polynomial& a,
		    const Fp_polynomial& b);
void plain_div (Fp_polynomial& q, const Fp_polynomial& a,
		const Fp_polynomial& b);
void plain_rem (Fp_polynomial& r, const Fp_polynomial& a,
		const Fp_polynomial& b);

	// These always use FFT arithmetic
void fft_div_rem (Fp_polynomial& q, Fp_polynomial& r,
		  const Fp_polynomial& a,
		  const Fp_polynomial& b);
void fft_div (Fp_polynomial& q, const Fp_polynomial& a,
	      const Fp_polynomial& b);
void fft_rem (Fp_polynomial& r, const Fp_polynomial& a,
	      const Fp_polynomial& b);

void plain_inv (Fp_polynomial& x, const Fp_polynomial& a,
		lidia_size_t m);
	// always uses "classical" algorithm
	// ALIAS RESTRICTION: input may not alias output

void newton_inv (Fp_polynomial& x, const Fp_polynomial& a,
		 lidia_size_t m);
	// uses a Newton Iteration with the FFT.
	// ALIAS RESTRICTION: input may not alias output
#endif	// HEADBANGER



	//***************************************************************
	//
	//	Miscellaneaous
	//
	//****************************************************************/

void swap (Fp_polynomial& x, Fp_polynomial& y);
	// swap x & y (only pointers are swapped)

void trunc (Fp_polynomial& c, const Fp_polynomial& a,
	    lidia_size_t m);
	// c = a % x^m

void shift_right (Fp_polynomial& c, const Fp_polynomial& a,
		  lidia_size_t n);
	// c = a/x^n

void shift_left (Fp_polynomial& c, const Fp_polynomial& a,
		 lidia_size_t n);
	// c = a*x^n

void derivative (Fp_polynomial& c, const Fp_polynomial& a);
	// c = derivative of a
Fp_polynomial derivative (const Fp_polynomial& a);

void add_multiple (Fp_polynomial &f,
		   const Fp_polynomial &g,
		   const bigint &s,
		   lidia_size_t n,
		   const Fp_polynomial &h);
	// f = g + s*x^n*h, n >= 0

void cyclic_reduce (Fp_polynomial& x, const Fp_polynomial& a,
		    lidia_size_t m);
	// computes c = a mod x^m-1

void copy_reverse (Fp_polynomial& x, const Fp_polynomial& a,
		   lidia_size_t lo, lidia_size_t hi);
	// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
	// input may not alias output

void multiply_by_x_mod (Fp_polynomial& c,
			const Fp_polynomial& a,
			const Fp_polynomial& f);
	// c = (a * x) mod f


void remainder (Fp_polynomial &, const Fp_polynomial &,
		const Fp_poly_modulus &);

//***************************************************************
//
//                      operators
//
//***************************************************************

Fp_polynomial operator - (const Fp_polynomial &a);
Fp_polynomial operator + (const Fp_polynomial &a, const Fp_polynomial &b);
Fp_polynomial operator + (const bigint &a, const Fp_polynomial &b);
Fp_polynomial operator + (const Fp_polynomial &a, const bigint &b);
Fp_polynomial operator - (const Fp_polynomial &a, const Fp_polynomial &b);
Fp_polynomial operator - (const bigint &a, const Fp_polynomial &b);
Fp_polynomial operator - (const Fp_polynomial &a, const bigint &b);
Fp_polynomial operator * (const Fp_polynomial &a, const Fp_polynomial &b);
Fp_polynomial operator * (const bigint &a, const Fp_polynomial &b);
Fp_polynomial operator * (const Fp_polynomial &a, const bigint &b);
Fp_polynomial operator / (const Fp_polynomial &a, const Fp_polynomial &b);
Fp_polynomial operator / (const bigint &a, const Fp_polynomial &b);
Fp_polynomial operator / (const Fp_polynomial &a, const bigint &b);
Fp_polynomial operator % (const Fp_polynomial &a, const Fp_polynomial &b);


// for class factorization < Fp_polynomial > :
//	used in class single_factor < Fp_polynomial > :
bool operator < (const Fp_polynomial &a, const Fp_polynomial &b);
bool operator <= (const Fp_polynomial &a, const Fp_polynomial &b);
bool operator == (const Fp_polynomial& a, const Fp_polynomial& b);
bool operator != (const Fp_polynomial& a, const Fp_polynomial& b);

void negate (Fp_polynomial& x, const Fp_polynomial& a);
void add (Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b);
void add (Fp_polynomial & x, const Fp_polynomial& a, const bigint& b);
inline
void add (Fp_polynomial& x, const bigint& a, const Fp_polynomial& b)
{
    add(x, b, a);
}

void subtract (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_polynomial& b);
void subtract (Fp_polynomial & x, const Fp_polynomial& a, const bigint& b);
inline
void subtract (Fp_polynomial& x, const bigint& a, const Fp_polynomial& b)
{
    LiDIA::negate(x, b);
    add(x, x, a);
}


void multiply (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_polynomial& b);
inline
void multiply (Fp_polynomial & x, const Fp_polynomial& a, const bigint& b)
{
    multiply_by_scalar(x, a, b);
}
inline
void multiply (Fp_polynomial & x, const bigint& a, const Fp_polynomial& b)
{
    multiply_by_scalar(x, b, a);
}
void multiply_by_scalar (Fp_polynomial &c, const Fp_polynomial &a,
			 const bigint &b);

// These always use "classical" arithmetic
void plain_mul (Fp_polynomial& x, const Fp_polynomial& a,
		const Fp_polynomial& b);
void plain_sqr (Fp_polynomial& x, const Fp_polynomial& a);

// These always use FFT arithmetic
void fft_mul (Fp_polynomial& x, const Fp_polynomial& a,
	      const Fp_polynomial& b);
void fft_sqr (Fp_polynomial& x, const Fp_polynomial& a);

void square(Fp_polynomial& x, const Fp_polynomial& a);

// q = a/b, r = a%b
void div_rem (Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a,
	      const Fp_polynomial& b); 

// q = a/b
void divide (Fp_polynomial& q, const Fp_polynomial& a, const Fp_polynomial& b);
void divide (Fp_polynomial& q, const bigint& a, const Fp_polynomial& b);
void divide (Fp_polynomial& q, const Fp_polynomial& a, const bigint& b);

// r = a%b
void remainder (Fp_polynomial& r, const Fp_polynomial& a,
		const Fp_polynomial& b);


// computes c = a^{-1} % x^m
// constant term must be non-zero
void invert (Fp_polynomial& c, const Fp_polynomial& a, lidia_size_t m);

// These always use "classical" arithmetic
void plain_div_rem (Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a,
		    const Fp_polynomial& b);
void plain_div (Fp_polynomial& q, const Fp_polynomial& a,
		const Fp_polynomial& b);
void plain_rem (Fp_polynomial& r, const Fp_polynomial& a,
		const Fp_polynomial& b);

// These always use FFT arithmetic
void fft_div_rem (Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a,
		  const Fp_polynomial& b);
void fft_div (Fp_polynomial& q, const Fp_polynomial& a,
	      const Fp_polynomial& b);
void fft_rem (Fp_polynomial& r, const Fp_polynomial& a,
	      const Fp_polynomial& b);

// always uses "classical" algorithm
// ALIAS RESTRICTION: input may not alias output
void plain_inv (Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m);

// uses a Newton Iteration with the FFT.
// ALIAS RESTRICTION: input may not alias output
void newton_inv (Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m);


//***************************************************************
//
//	Miscellaneaous
//
//****************************************************************/

// swap x & y (only pointers are swapped)
void swap (Fp_polynomial& x, Fp_polynomial& y);

// c = a % x^m
void trunc (Fp_polynomial& c, const Fp_polynomial& a, lidia_size_t m);

// c = a/x^n
void shift_right (Fp_polynomial& c, const Fp_polynomial& a, lidia_size_t n);

// c = a*x^n
void shift_left (Fp_polynomial& c, const Fp_polynomial& a, lidia_size_t n);

// c = derivative of a
void derivative (Fp_polynomial& c, const Fp_polynomial& a);
inline
Fp_polynomial derivative (const Fp_polynomial& a)
{
    Fp_polynomial x;
    
    derivative(x, a);
    return x;
}

// f = g + s*x^n*h, n >= 0
void add_multiple (Fp_polynomial &f, const Fp_polynomial &g,
		   const bigint &s, lidia_size_t n, const Fp_polynomial &h);

// computes c = a mod x^m-1
void cyclic_reduce (Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m);

// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
// input may not alias output
void copy_reverse (Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t lo,
		   lidia_size_t hi);

// c = (a * x) mod f
void multiply_by_x_mod (Fp_polynomial& c, const Fp_polynomial& a,
			const Fp_polynomial& f);

//set_degree
void remainder (Fp_polynomial &, const Fp_polynomial &,
		const Fp_poly_modulus &);


//***************************************************************
//
//    			    gcd's
//
//***************************************************************

// x = gcd(a, b), x is always monic (or zero if a == b == 0).
void gcd (Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b);

inline Fp_polynomial
gcd(const Fp_polynomial& a, const Fp_polynomial& b)
{
	Fp_polynomial x;

	gcd(x, a, b);
	return x;
}



void
xgcd (Fp_polynomial& d, Fp_polynomial& s, Fp_polynomial& t,
      const Fp_polynomial& a, const Fp_polynomial& b);
// d = gcd(a, b), a s + b t = d

void
xgcd_left (Fp_polynomial& d, Fp_polynomial& s, const Fp_polynomial& a,
	   const Fp_polynomial& b);
// d = gcd(a, b), a s + b ?? = d

inline void
xgcd_right (Fp_polynomial& d, Fp_polynomial& t, const Fp_polynomial& a,
	    const Fp_polynomial& b)
{
	xgcd_left(d, t, b, a);
}



// d = gcd(a, b), a ?? + b t = d

void
resultant (bigint &, const Fp_polynomial &, const Fp_polynomial &);


inline bigint
resultant (const Fp_polynomial& a, const Fp_polynomial& b)
{
	bigint x;
	resultant (x, a, b);
	return x;
}



#ifndef HEADBANGER
// These always use "classical" arithmetic
void
plain_xgcd (Fp_polynomial& d, Fp_polynomial& s, Fp_polynomial& t,
	    const Fp_polynomial& a, const Fp_polynomial& b);
void
plain_xgcd_left (Fp_polynomial& d, Fp_polynomial& s,
		 const Fp_polynomial& a, const Fp_polynomial& b);
void
plain_gcd (Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b);
#endif	// HEADBANGER


//***************************************************************
//
//	Modular Arithmetic without pre-conditioning
//
//***************************************************************

#ifndef HEADBANGER
// arithmetic mod f.
// all inputs and outputs are polynomials of degree less than deg(f).
// ASSUMPTION: f is assumed monic, and deg(f) > 0.
// NOTE: if you want to do many computations with a fixed f,
//       use the Fp_poly_modulus data structure and associated routines below.

void
multiply_mod (Fp_polynomial& c, const Fp_polynomial& a, const Fp_polynomial& b, const Fp_polynomial& f);
// c = (a * b) % f

void
square_mod (Fp_polynomial& c, const Fp_polynomial& a, const Fp_polynomial& f);
// c = a^2 % f

void
multiply_by_x_mod (Fp_polynomial& c, const Fp_polynomial& a, const Fp_polynomial& f);
// c = (a * x) mod f

void
invert_mod (Fp_polynomial& c, const Fp_polynomial& a, const Fp_polynomial& f);
// c = a^{-1} % f, error if a is not invertible

bool
invert_mod_status(Fp_polynomial& c, const Fp_polynomial& a, const Fp_polynomial& f);
// if (a, f) = 1, returns 1 and sets c = a^{-1} % f
// otherwise, returns 0 and sets c = (a, f)

void
power_mod (Fp_polynomial& c, const Fp_polynomial& a, const bigint&  e, const Fp_polynomial& f);
// c = a^e % f
// WARNING: obsolete.  Use power_mod with Fp_poly_modulus (see below).

void
power_x_mod (Fp_polynomial& c, const bigint&  e, const Fp_polynomial& f);
//c = x^e mod f
// WARNING: obsolete.  Use power_mod with Fp_poly_modulus (see below).

void
power_x_plus_a_mod (Fp_polynomial& c, const bigint& a, const bigint&  e, const Fp_polynomial& f);
// c = (x + a)^e mod f
// WARNING: obsolete.  Use power_mod with Fp_poly_modulus (see below).

#endif 	// HEADBANGER



//***************************************************************
//
//    			    I/O
//
//***************************************************************

std::istream& operator >> (std::istream& s, Fp_polynomial& x);
std::ostream& operator << (std::ostream& s, const Fp_polynomial& a);



//***************************************************************
//
//			Miscellaneaous
//
//***************************************************************

void randomize (Fp_polynomial& x, const bigint& p, lidia_size_t n);
// generate a random polynomial of degree = n with coefficients in Z/pZ

void power (Fp_polynomial &x, const Fp_polynomial &a, lidia_size_t e);
// x = a^e, e >= 0



//*********************************************************************
//
//			miscellaneous
//
//*********************************************************************


// ************   min_poly, checked_min_poly, irred_poly   ************

void
prob_min_poly (Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
	       const Fp_poly_modulus& F);
// computes the monic minimal polynomial of (g mod f).
// m = a bound on the degree of the minimal polynomial.
// The algorithm is probabilistic, always returns a divisor of
// the minimal polynomial, and returns a proper divisor with
// probability at most m/p.


void min_poly (Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
	       const Fp_poly_modulus& F);
// same as above, but guarantees that result is correct


void irred_poly (Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
		 const Fp_poly_modulus& F);
// same as above, but assumes that f is irreducible
// (or at least that the minimal poly of g is itself irreducible).
// The algorithm is deterministic (and hence is always correct).



// ************   Multiplication of several polynomials   ************
//void multiply(Fp_polynomial& x, base_vector < Fp_polynomial > & a);
    // x = product of all a[i]'s, contents of a[i] are destroyed



// ************   Modular Composition   ************
// algorithms for computing g(h) mod f

void compose (Fp_polynomial& x, const Fp_polynomial& g,
	      const Fp_polynomial& h, const Fp_poly_modulus& F);
// x = g(h) mod f


void compose2 (Fp_polynomial& x1, Fp_polynomial& x2,
	       const Fp_polynomial& g1, const Fp_polynomial& g2,
	       const Fp_polynomial& h, const Fp_poly_modulus& F);
// xi = gi(h) mod f (i = 1, 2)
// ALIAS RESTRICTION:  xi may not alias gj, for i != j


void compose3 (Fp_polynomial& x1, Fp_polynomial& x2, Fp_polynomial& x3,
	       const Fp_polynomial& g1, const Fp_polynomial& g2, const Fp_polynomial& g3,
	       const Fp_polynomial& h, const Fp_poly_modulus& F);
// xi = gi(h) mod f (i = 1..3)
// ALIAS RESTRICTION:  xi may not alias gj, for i != j



// ************   update_map   ************
void update_map (base_vector< bigint > & x, const base_vector< bigint > & a,
		 const Fp_poly_multiplier& B, const Fp_poly_modulus& F);
// computes (a, b), (a, (b*X)%f), ..., (a, (b*X^{n-1})%f),
// where (,) denotes the vector inner product.
// This is really a "transposed" MulMod by B.

void update_map (base_vector< bigint > & xx, const base_vector< bigint > & a,
		 const Fp_polynomial& b, const Fp_polynomial& f);
// same as above, but uses only classical arithmetic



// ************   inner_product   ************
void inner_product (bigint& x, const base_vector< bigint > & a,
		    const Fp_polynomial &b, lidia_size_t offset = 0);



//*********************************************************************
//
//			factorization.cc
//
//*********************************************************************

void rec_find_roots (base_vector< bigint > & x, const Fp_polynomial& f);

base_vector< bigint > find_roots (const Fp_polynomial & f, int flag = 0);
// returns the list of roots of f (without multiplicities)
// if (flag != 0), f must be monic and the product of deg(f) distinct
//  linear factors
// otherwise, no assumptions on f are made

bigint find_root (const Fp_polynomial & ff);
// finds a single root of ff.
// assumes that ff is monic and splits into distinct linear factors

Fp_polynomial find_factor( const Fp_polynomial & ff, lidia_size_t d ); 
    // Finds an irreducible factor of ff of degree d 
    // Conditions:  ff splits into distinct irreducible factors of degree d.

Fp_polynomial find_factor_degree_d( const Fp_polynomial & ff, lidia_size_t d ); 
    // Finds a factor of ff of degree in [d, 2d]
    // Conditions:  ff splits into distinct irreducible factors of degree 1.

bool prob_irred_test (const Fp_polynomial & f, lidia_size_t iter = 1);
// performs a fast, probabilistic irreduciblity test
// the test can err only if f is reducible, and the
// error probability is bounded by p^{-iter}.
// Works for any p.


bool det_irred_test (const Fp_polynomial & f);
// performs a recursive deterministic irreducibility test


bool iter_irred_test (const Fp_polynomial & f);
// performs an iterative deterministic irreducibility test,
// based on DDF



void build_irred (Fp_polynomial& f, const bigint & p, lidia_size_t n);
// Build a monic irreducible poly of degree n.


void build_random_irred (Fp_polynomial & f, const Fp_polynomial & g);
// g is a monic irreducible polynomial.
// constructs a random monic irreducible polynomial f of the same degree.

lidia_size_t compute_degree (const Fp_polynomial & h, const Fp_poly_modulus & F,
			     lidia_size_t d = -1);
// f = F.f is assumed to be an "equal degree" polynomial
// h = x^p mod f
// the common degree of the irreducible factors of f is computed
// d is multiple of common degree, if -1, choose degree F.f
// This routine is useful in counting points on elliptic curves


lidia_size_t prob_compute_degree (const Fp_polynomial & h, const Fp_poly_modulus & F);
// same as above, but uses a slightly faster probabilistic algorithm
// the return value may be 0 or may be too big, but for large p
// (relative to n), this happens with very low probability.



void trace_map (Fp_polynomial & w, const Fp_polynomial & a, lidia_size_t d,
		const Fp_poly_modulus & F, const Fp_polynomial & b);
// w = a+a^q+...+^{q^{d-1}} mod f;
// it is assumed that d >= 0, and b = x^q mod f, q a power of p
// Space allocation can be controlled via ComposeBound



void power_compose (Fp_polynomial& w, const Fp_polynomial& b, lidia_size_t d,
		    const Fp_poly_modulus& F);
// w = x^{q^d} mod f;
// it is assumed that d >= 0, and b = x^q mod f, q a power of p
// Space allocation can be controlled via ComposeBound



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_FP_POLYNOMIAL

#include	"LiDIA/specialization/Fp_polynomial.special"



#endif	// LIDIA_FP_POLYNOMIAL_H_GUARD_
