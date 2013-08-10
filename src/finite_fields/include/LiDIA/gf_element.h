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
//	Author	: Detlef Anton (DA), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_GF_ELEMENT_H_GUARD_
#define LIDIA_GF_ELEMENT_H_GUARD_


#ifndef LIDIA_GALOIS_FIELD_H_GUARD_
# include	"LiDIA/galois_field.h"
#endif
#ifndef LIDIA_MULTI_BIGMOD_H_GUARD_
# include	"LiDIA/multi_bigmod.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class galois_field_rep;
class gf_element;

template <class T >
class polynomial;



class gf_element
{
	friend class galois_field_rep;


	// immediate access for class gf_polynomial:
	const Fp_polynomial& polynomial_rep2() const;
	friend void to_Kronecker(Fp_polynomial &g, const polynomial< gf_element > &f,
				 lidia_size_t lo, lidia_size_t hi);


	//
	// the C++ type we use to represent an Element
	// in a Polynomial Base of a Galois field
	//

	static galois_field uninitialized_field;
	static unsigned int output_format;
	enum {
		SHORT_OUTPUT = 0,
		VERBOSE_OUTPUT = 1
	};


	galois_field ff;
	void * rep;

	galois_field_rep const* get_ff_rep() const
	{
		return ff.rep;
	}
	void set_field(const galois_field&);
	gf_element promote(const bigint&) const;
	gf_element promote(const gf_element&) const;
	void power(const bigint&);

public:

	//
	// constructors and destructor
	//

	gf_element();
	gf_element(const gf_element&);
	gf_element(const galois_field&);
	//gf_element(const Fp_polynomial&);
	~gf_element();

	//
	// assignments
	//

	gf_element& operator = (const gf_element&);
	gf_element& operator = (const bigint&);
	void assign(const gf_element &);
	void assign(const bigint&);
	void assign_zero();
	void assign_zero(const galois_field&);
	void assign_one();
	void assign_one(const galois_field&);

	//
	// access functions
	//

	const galois_field& get_field() const;
	bigint characteristic () const;
	Fp_polynomial polynomial_rep() const;
	void set_polynomial_rep(const Fp_polynomial&);


	//
	// arithmetical operations
	//
	gf_element & operator += (const gf_element &);
	gf_element & operator -= (const gf_element &);
	gf_element & operator *= (const gf_element &);
	gf_element & operator /= (const gf_element &);

	gf_element & operator += (const bigint &);
	gf_element & operator -= (const bigint &);
	gf_element & operator *= (const bigint &);
	gf_element & operator /= (const bigint &);

	//
	// procedural versions
	//

	void negate();
	void multiply_by_2();
	void divide_by_2();
	void invert();

	friend void add     (gf_element&, const gf_element&, const gf_element&);
	friend void subtract(gf_element&, const gf_element&, const gf_element&);
	friend void multiply(gf_element&, const gf_element&, const gf_element&);
	friend void divide  (gf_element&, const gf_element&, const gf_element&);

	friend void add     (gf_element&, const bigint&, const gf_element&);
	friend void subtract(gf_element&, const bigint&, const gf_element&);
	friend void multiply(gf_element&, const bigint&, const gf_element&);
	friend void divide  (gf_element&, const bigint&, const gf_element&);

	friend void add     (gf_element&, const gf_element&, const bigint&);
	friend void subtract(gf_element&, const gf_element&, const bigint&);
	friend void multiply(gf_element&, const gf_element&, const bigint&);
	friend void divide  (gf_element&, const gf_element&, const bigint&);

	friend void negate  (gf_element&, const gf_element&);
	friend void invert  (gf_element&, const gf_element&);

	friend void power(gf_element &, const gf_element &, const bigint &);
	friend void pth_power(gf_element &, const gf_element &, lidia_size_t);

	friend void square(gf_element &, const gf_element &);
	friend gf_element sqrt(const gf_element &);

	//
	// comparisons
	//

	bool operator == (const gf_element &) const;
	bool operator != (const gf_element &) const;
	bool operator == (const bigint &) const;
	bool operator != (const bigint &) const;

	bool is_zero() const;
	bool is_one()  const;

	//
	// basic functions
	//

	void randomize()
	{
		randomize(ff.degree());
	}
	void randomize(int unsigned deg);
	friend void swap(gf_element &, gf_element &);

	//
	// high level functions
	///

	bigint       order() const;
	multi_bigmod trace() const;
	multi_bigmod norm()  const;
	bool is_primitive_element() const;
	bool is_free_element() const;
	bool is_square() const;
	
	void assign_primitive_element(const galois_field&);


	bool solve_quadratic(const gf_element&, const gf_element&);
	unsigned int relative_degree() const;
	unsigned int absolute_degree() const;


	bigint lift_to_Z() const;
	friend udigit hash(const gf_element &);
#if 0
	friend unsigned int in_subfield(const gf_element &);
#endif

	//
	// input / output
	//

	void input(std::istream&);
	void output(std::ostream&) const;
	bigint return_as_bigint() const;

	static void set_output_format(unsigned int);
	static unsigned int get_output_format();
};



gf_element operator - (const gf_element &);
gf_element operator + (const gf_element &, const gf_element &);
gf_element operator - (const gf_element &, const gf_element &);
gf_element operator * (const gf_element &, const gf_element &);
gf_element operator / (const gf_element &, const gf_element &);

gf_element operator + (const bigint &, const gf_element &);
gf_element operator - (const bigint &, const gf_element &);
gf_element operator * (const bigint &, const gf_element &);
gf_element operator / (const bigint &, const gf_element &);

gf_element operator + (const gf_element &, const bigint &);
gf_element operator - (const gf_element &, const bigint &);
gf_element operator * (const gf_element &, const bigint &);
gf_element operator / (const gf_element &, const bigint &);


gf_element inverse(const gf_element &);

void swap(gf_element &, gf_element &);

std::istream & operator >> (std::istream & in, gf_element &);
std::ostream & operator << (std::ostream & out, const gf_element &);

void add     (gf_element&, const gf_element&, const gf_element&);
void subtract(gf_element&, const gf_element&, const gf_element&);
void multiply(gf_element&, const gf_element&, const gf_element&);
void divide  (gf_element&, const gf_element&, const gf_element&);

void add     (gf_element&, const bigint&, const gf_element&);
void subtract(gf_element&, const bigint&, const gf_element&);
void multiply(gf_element&, const bigint&, const gf_element&);
void divide  (gf_element&, const bigint&, const gf_element&);

void add     (gf_element&, const gf_element&, const bigint&);
void subtract(gf_element&, const gf_element&, const bigint&);
void multiply(gf_element&, const gf_element&, const bigint&);
void divide  (gf_element&, const gf_element&, const bigint&);

void negate  (gf_element&, const gf_element&);
void invert  (gf_element&, const gf_element&);

void power(gf_element &, const gf_element &, const bigint &);
void pth_power(gf_element &, const gf_element &, lidia_size_t);

void square(gf_element &, const gf_element &);
gf_element sqrt(const gf_element &);

#define LIDIA_CLASS_GF_ELEMENT



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF_ELEMENT_H_GUARD_
