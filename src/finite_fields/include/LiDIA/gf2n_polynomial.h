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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================

// Description  :  Polynomials over gf2n, basic operations
//                 This file essentially copies many functions directly
//                 from Fp_polynomial.


#ifndef LIDIA_GF2N_POLYNOMIAL_H_GUARD_
#define LIDIA_GF2N_POLYNOMIAL_H_GUARD_


#ifndef LIDIA_GF2N_H_GUARD_
# include	"LiDIA/gf2n.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif

#ifdef DEBUG
# include	<cassert>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf2n_polynomial
{
	friend class gf2n_poly_modulus;

private:

	static int default_size;
	static gf2n *gf2n_p_tmp;
	static int gf2n_p_top;
	static int gf2n_p_tmpsize;
	static int refc;

	int size;
	gf2n *coeff;
	int deg;


public:

	gf2n_polynomial(); // = 0

	gf2n_polynomial(unsigned int d); // = 1 * X^d

	gf2n_polynomial(const gf2n & f); // = f * X^0

	gf2n_polynomial(const gf2n_polynomial & p); // = p

	~gf2n_polynomial()
	{
		gf2n_polynomial::refc --;
		if (size > -1)
			delete[] coeff;

		if (!refc) {
			delete[] gf2n_p_tmp;
			gf2n_p_tmp = NULL;
			gf2n_p_top = gf2n_p_tmpsize = 0;
		}
	}

	//****  comparisons of polynomials ******************************


	friend bool operator == (const gf2n_polynomial &a,
				 const gf2n_polynomial &b);

	friend bool operator != (const gf2n_polynomial &a,
				 const gf2n_polynomial &b)
	{
		return (!(a == b));
	}

	//****** set defaultsize ***************************************

	static void set_defaultsize(unsigned int d)
	{
		if (d > 4)
			default_size = d;
		else
			default_size = 4;
	}

	void set_size(unsigned int d);

	//****** basic functions **************************************

	gf2n lead_coeff() const
	{
		if (deg >= 0)
			return coeff[deg];
		else return gf2n();
	}

	friend gf2n lead_coeff(const gf2n_polynomial & a)
	{
		return a.lead_coeff();
	}

	gf2n const_term() const
	{
		if (deg >= 0)
			return coeff[0];
		else
			return gf2n();
	}

	friend gf2n const_term(const gf2n_polynomial & a)
	{
		return a.const_term();
	}

	friend gf2n get_coefficient(const gf2n_polynomial & f, unsigned int i)
	{
		if (static_cast<int>(i) > f.deg) {
			lidia_error_handler("gf2n_polynomial", "i > degree(f)");
			return gf2n();
		}
		else
			return f.coeff[i];
	}

	void make_monic();

	friend int degree(const gf2n_polynomial &a)
	{
		return (a.deg);
	}

	unsigned int degree_of_definition() const;
	// return the minimal field degree over GF(2) that holds all coefficients

	int degree() const
	{
		return deg;
	}

	bool is_zero() const
	{
		return (deg < 0);
	}

	bool is_one() const
	{
		return (deg == 0 && coeff[0].is_one());
	}

	bool is_x() const
	{
		if (deg == 1)
			return (coeff[0].is_zero() && coeff[1].is_one());
		else
			return false;
	}

	bool is_monic() const
	{
		if (deg >= 0)
			return coeff[deg].is_one();
		else
			return false;
	}

	bool is_square() const;

	friend gf2n eval(const gf2n_polynomial & f, const gf2n & x);

	gf2n operator() (const gf2n & x) const
	{
		return eval(*this, x);
	}

	friend bool sqrt(gf2n_polynomial &f, const gf2n_polynomial &g);


	//************* assignments *******************************

	void assign(const gf2n_polynomial & f);

	void assign(const base_vector< gf2n > & v);

	void assign_zero()
	{
		deg = -1;
	}

	void assign_one()
	{
		if (size == -1)
			set_size(gf2n_polynomial::default_size);
		coeff[0].assign_one();
		deg = 0;
	}


	void assign_x()
	{
		(*this).assign(gf2n_polynomial(1));
	}

	void kill()
	{
		if (size > -1)
			delete[] coeff;
		size = -1;
		deg = -1;
	}

	void randomize(unsigned int degree);

	friend void randomize(gf2n_polynomial h, unsigned int degree)
	{
		h.randomize(degree);
	}

	gf2n_polynomial & operator = (const gf2n_polynomial &a)
	{
		assign(a);
		return *this;
	}

	void set_coefficient(const gf2n & x, unsigned int i);
	void set_coefficient(unsigned int);
	void get_coefficient(gf2n &, unsigned int i) const;

	gf2n operator [] (unsigned int i) const
	{
		gf2n x;

		get_coefficient(x, i);
		return x;
	}



	//************* operations *************************************

	// g = f * (X+a)
	friend void multiply_by_linear(gf2n_polynomial & g,
				       const gf2n_polynomial & f,
				       const gf2n & a);



	friend void multiply_by_scalar(gf2n_polynomial & g, const gf2n & b,
				       const gf2n_polynomial & a);


	friend void multiply_by_scalar(gf2n_polynomial & g,
				       const gf2n_polynomial & a,
				       const gf2n & b)
	{
		multiply_by_scalar(g, b, a);
	}


	friend void negate (gf2n_polynomial & c, const gf2n_polynomial & a)
		// for consistancy
	{
		c.assign(a);
	}

	friend void add (gf2n_polynomial &, const gf2n_polynomial &,
			 const gf2n_polynomial &);


	friend void subtract (gf2n_polynomial &c, const gf2n_polynomial &a,
			      const gf2n_polynomial &b)
		// for consistancy
	{
		add(c, a, b);
	}

	friend gf2n_polynomial operator + (const gf2n_polynomial & a,
					   const gf2n_polynomial & b)
	{
		gf2n_polynomial c;

		add(c, a, b);
		return(c);
	}


	friend gf2n_polynomial operator - (const gf2n_polynomial & a,
					   const gf2n_polynomial & b)
	{
		gf2n_polynomial c;

		add(c, a, b);
		return(c);
	}
	gf2n_polynomial & operator += (const gf2n_polynomial &a)
	{
		add(*this, *this, a);
		return *this;
	}

	gf2n_polynomial & operator -= (const gf2n_polynomial &a)
	{
		add(*this, *this, a);
		return *this;
	}



	friend void plain_mul(gf2n_polynomial &, const gf2n_polynomial &,
			      const gf2n_polynomial &);

	friend void multiply(gf2n_polynomial &, const gf2n_polynomial &,
			     const gf2n_polynomial &);



	friend void square(gf2n_polynomial&, const gf2n_polynomial&);

	friend gf2n_polynomial square (const gf2n_polynomial & a)
	{
		gf2n_polynomial c;

		square(c, a);
		return(c);
	}

	friend gf2n_polynomial operator * (const gf2n_polynomial & a,
					   const gf2n_polynomial & b)
	{
		gf2n_polynomial c;

		multiply(c, a, b);
		return(c);
	}

	friend gf2n_polynomial operator * (const gf2n_polynomial & a,
					   const gf2n & b)
	{
		gf2n_polynomial c;

		multiply_by_scalar(c, b, a);
		return c;
	}

	friend gf2n_polynomial operator * (const gf2n & b,
					   const gf2n_polynomial & a)
	{
		gf2n_polynomial c;

		multiply_by_scalar(c, b, a);
		return c;
	}


	gf2n_polynomial & operator *= (const gf2n_polynomial &a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	friend void div_rem(gf2n_polynomial & div, gf2n_polynomial & rem,
			    const gf2n_polynomial & a, const gf2n_polynomial & b)
	{
		if (&div != &b) {
			divide(div, a, b);
			add(rem, rem, div * b);
		}
		else {
			gf2n_polynomial t(b);
			divide(div, a, t);
			add(rem, rem, div * t);
		}
	}

	friend void plain_div_rem(gf2n_polynomial & div, gf2n_polynomial & rem,
				  const gf2n_polynomial &, const gf2n_polynomial &);

	friend void divide(gf2n_polynomial & div, const gf2n_polynomial & a,
			   const gf2n_polynomial & b);

	friend void plain_div(gf2n_polynomial & div, const gf2n_polynomial & a,
			      const gf2n_polynomial & b);

	friend void kara_div(gf2n_polynomial & q, const gf2n_polynomial & a,
			     const gf2n_polynomial & b);


	friend void remainder (gf2n_polynomial &r, const gf2n_polynomial &a,
			       const gf2n_polynomial &b);

	friend void plain_rem (gf2n_polynomial &, const gf2n_polynomial &,
			       const gf2n_polynomial &);

	friend void kara_rem(gf2n_polynomial & r, const gf2n_polynomial & a,
			     const gf2n_polynomial & b);

	friend gf2n_polynomial  operator / (const gf2n_polynomial & a,
					    const gf2n_polynomial & b)
	{
		gf2n_polynomial q; divide(q, a, b); return q;
	}

	gf2n_polynomial & operator /= (const gf2n_polynomial &a)
	{
		gf2n_polynomial h;

		div_rem(*this, h, *this, a);
		return *this;
	}

	friend gf2n_polynomial  operator % (const gf2n_polynomial & a,
					    const gf2n_polynomial & b)
	{
		gf2n_polynomial q; remainder(q, a, b); return q;
	}

	gf2n_polynomial & operator %= (const gf2n_polynomial &a)
	{
		remainder(*this, *this, a);
		return *this;
	}



private:

	friend void copy_reverse(gf2n_polynomial & x, const gf2n_polynomial & a,
				 lidia_size_t lo, lidia_size_t hi);

public:

	friend void invert(gf2n_polynomial & x, const gf2n_polynomial & a,
			   lidia_size_t m);

	friend void plain_inv(gf2n_polynomial& x, const gf2n_polynomial& a,
			      lidia_size_t m);

	friend void newton_inv(gf2n_polynomial& x, const gf2n_polynomial& a,
			       lidia_size_t m);


	//***************************************************************
	//
	//	Miscellaneaous
	//
	//***************************************************************

	// g = f * X^d
	friend void shift_left(gf2n_polynomial & g,
			       const gf2n_polynomial & f,
			       unsigned int d = 1);

	// g = f / X^d
	friend void shift_right(gf2n_polynomial & g,
				const gf2n_polynomial & f,
				unsigned int d = 1);

	// g = floor(f / X^d)
	friend void floor(gf2n_polynomial & g, const gf2n_polynomial & f,
			  unsigned int d = 1);

	friend void swap(gf2n_polynomial &, gf2n_polynomial &);

	friend void trunc(gf2n_polynomial& c, const gf2n_polynomial& a,
			  unsigned int m)    // c = a % x^m
	{
		c.assign(a);
		if (static_cast<int>(m-1) < a.deg)
			c.deg = m-1;
	}

	friend void derivative(gf2n_polynomial& c, const gf2n_polynomial& a);

	friend gf2n_polynomial derivative(const gf2n_polynomial& a)
	{
		gf2n_polynomial x;

		derivative(x, a);
		return x;
	}


	// memory is not freed, has to be done by calling application
	friend gf2n_polynomial *power_table(const gf2n_polynomial & a,
					    int d);


	//*********** higher level functions as gcd

	friend void gcd(gf2n_polynomial &, const gf2n_polynomial &,
			const gf2n_polynomial &);

	// d = gcd(a, b), a s + b t = d
	friend void xgcd(gf2n_polynomial& d, gf2n_polynomial& s, gf2n_polynomial& t,
			 const gf2n_polynomial& a, const gf2n_polynomial& b);

	friend gf2n_polynomial xgcd_left(gf2n_polynomial& s,
					 const gf2n_polynomial & a,
					 const gf2n_polynomial & b);

	gf2n_polynomial xgcd_right(gf2n_polynomial& t,
				   const gf2n_polynomial & a,
				   const gf2n_polynomial & b)
	{
		return xgcd_left(t, b, a);
	}

	void resultant(gf2n &, const gf2n_polynomial &,
		       const gf2n_polynomial &);

	gf2n resultant (const gf2n_polynomial& a, const gf2n_polynomial& b)
	{
		gf2n x;

		resultant (x, a, b);
		return x;
	}

	void randomize(gf2n_polynomial & x, lidia_size_t n);
	// generate a random polynomial of degree = n

	void power(gf2n_polynomial &x, const gf2n_polynomial &a, lidia_size_t e);
	// x = a^e, e >= 0


	//***************************************************************

	//  Modular Arithmetic without pre-conditioning

	//***************************************************************
	// arithmetic mod f.
	// all inputs and outputs are polynomials of degree less than deg(f).
	// ASSUMPTION: f is assumed monic, and deg(f) > 0.
	// NOTE: if you want to do many computations with a fixed f,
	//       use the gf2n_poly_modulus

	friend void multiply_mod(gf2n_polynomial & r, const gf2n_polynomial & a,
				 const gf2n_polynomial & b,
				 const gf2n_polynomial & f);

	friend void square_mod(gf2n_polynomial & r, const gf2n_polynomial & a,
			       const gf2n_polynomial & f);

	//if (a, f) = 1, returns 1 and sets c = a^{-1} % f, otherwise error
	friend void invert_mod(gf2n_polynomial & r, const gf2n_polynomial & a,
			       const gf2n_polynomial & f);

	// if (a, f) = 1, returns 1 and sets c = a^{-1} % f
	// otherwise, returns 0 and sets c = (a, f)

	friend bool invert_mod_status(gf2n_polynomial& c, const gf2n_polynomial& a,
				      const gf2n_polynomial& f);

	// c = (a * x) mod f
	friend void multiply_by_x_mod(gf2n_polynomial & c, const gf2n_polynomial & a,
				      const gf2n_polynomial & f)
	{
		shift_left(c, a, 1);
		remainder(c, c, f);
	}


	//********** functions for factoring polynomials ******************


	friend void find_all_roots (base_vector< gf2n > &, const gf2n_polynomial &,
				    unsigned int d = gf2n::extension_degree());

	friend gf2n find_root (const gf2n_polynomial &,
			       unsigned int d = gf2n::extension_degree());

	void inner_product(gf2n & x, const base_vector< gf2n > & a,
			   const gf2n_polynomial &b, lidia_size_t offset = 0);



	//************* input / output *********************************************

	friend std::ostream & operator << (std::ostream &, const gf2n_polynomial &);

	friend std::istream & operator >> (std::istream &, gf2n_polynomial &);


	//********* friend functions from gf2n_poly_modulus ***************

	friend void plain_rem(class gf2n_polynomial &, const class gf2n_polynomial &,
			      class gf2n_poly_modulus &);

	friend void kara_rem (class gf2n_polynomial &, const class gf2n_polynomial &,
			      class gf2n_poly_modulus &);

	friend void multiply(class gf2n_polynomial &, const class gf2n_polynomial &,
			     const class gf2n_polynomial &,
			     class gf2n_poly_modulus &);

	friend void square(class gf2n_polynomial &, const class gf2n_polynomial &,
			   class gf2n_poly_modulus &);

	friend void Xq(class gf2n_polynomial &, class gf2n_poly_modulus &,
		       unsigned int field_degree = gf2n::extension_degree());

	friend void compose(class gf2n_polynomial &, const class gf2n_polynomial &,
			    const class gf2n_polynomial &,
			    class gf2n_poly_modulus &);

	friend void trace_map(gf2n_polynomial & w, const gf2n_polynomial& a,
			      lidia_size_t d, gf2n_poly_modulus & F,
			      const gf2n_polynomial & b);

	friend void EDF_one_factor(gf2n_polynomial & w, const gf2n_polynomial & Xq,
				   gf2n_poly_modulus & Fpm, int d);

	friend lidia_size_t compute_degree(const class gf2n_polynomial &,
					   class gf2n_poly_modulus &, lidia_size_t);

private:

	static void resize_stack(int s = 0);
	static void resize_stack_for_degree(int d);
	friend void gf2n_p_mult(gf2n *r, const gf2n * a, const gf2n * b,
				int ad, int bd);
	friend void gf2n_p_mult32(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult33(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult42(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult43(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult44(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult52(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult53(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult54(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult55(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult62(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult63(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult64(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult65(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult66(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult72(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult73(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult74(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult75(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult76(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult77(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult82(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult83(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult84(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult85(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult86(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult87(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult88(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult92(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult93(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult94(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult95(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult96(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult97(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult98(gf2n * r, const gf2n * a, const gf2n * b);
	friend void gf2n_p_mult99(gf2n * r, const gf2n * a, const gf2n * b);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF2N_POLYNOMIAL_H_GUARD_
