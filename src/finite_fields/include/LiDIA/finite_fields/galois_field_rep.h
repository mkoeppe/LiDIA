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


#ifndef LIDIA_GALOIS_FIELD_REP_H_GUARD_
#define LIDIA_GALOIS_FIELD_REP_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif
#ifndef LIDIA_LIDIA_REFERENCE_COUNTER_H_GUARD_
# include	"LiDIA/lidia_reference_counter.h"
#endif
#ifndef LIDIA_INFO_GF2N_H_GUARD_
# include	"LiDIA/finite_fields/info_gf2n.h"
#endif
#ifndef LIDIA_FP_POLY_MODULUS_H_GUARD_
# include	"LiDIA/Fp_poly_modulus.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf_element;

class gf_flags
{
public:
	enum gf_rep_enum {
		UNKNOWN = 0,
		GF2n = 1,	// fields of char. 2
		GFp = 2,	// prime fields of odd char.
		GFpn = 3	// fields of odd char. with deg > 1
	};
	gf_rep_enum gf_rep;
};



class galois_field_rep: public gf_flags
{
private:
  friend class galois_field;
  
  mutable lidia_reference_counter rc;

public:
	friend class gf_element;
	friend void add(gf_element&, const gf_element&, const gf_element&);
	friend void subtract(gf_element&, const gf_element&, const gf_element&);
	friend void multiply(gf_element&, const gf_element&, const gf_element&);
	friend void square(gf_element&, const gf_element&);
	friend udigit hash(const gf_element&);

	//-----------------------------------------------------------------------------
	// function pointers for elements of finite fields
	//-----------------------------------------------------------------------------
	void (galois_field_rep::*construct)(gf_element&) const;
	void (galois_field_rep::*copy)(gf_element&, const gf_element&) const;
	void (galois_field_rep::*destruct)(gf_element&) const;
	void (galois_field_rep::*as0)(gf_element&) const;
	void (galois_field_rep::*as1)(gf_element&) const;
	bool (galois_field_rep::*iseq)(const gf_element&, const gf_element&) const;
	bool (galois_field_rep::*is0)(const gf_element&) const;
	bool (galois_field_rep::*is1)(const gf_element&) const;
	void (galois_field_rep::*neg)(gf_element&) const;
	void (galois_field_rep::*add)(gf_element&, const gf_element&, const gf_element&) const;
	void (galois_field_rep::*sub)(gf_element&, const gf_element&, const gf_element&) const;
	void (galois_field_rep::*mul)(gf_element&, const gf_element&, const gf_element&) const;
	void (galois_field_rep::*inv)(gf_element&) const;
	void (galois_field_rep::*sqr)(gf_element&, const gf_element&) const;
	udigit (galois_field_rep::*hash)(const gf_element&) const;
	const Fp_polynomial& (galois_field_rep::*get_pol_rep)(const gf_element&) const;
	void (galois_field_rep::*set_pol_rep)(gf_element&, const Fp_polynomial&) const;
	//void (galois_field_rep::*promote)(gf_element&, const bigint&) const;
	//-----------------------------------------------------------------------------
	// functions for finite prime fields of odd characteristic
	//-----------------------------------------------------------------------------
	void construct_GFp(gf_element&) const;
	void copy_GFp(gf_element&, const gf_element&) const;
	void destruct_GFp(gf_element&) const;
	void as0_GFp(gf_element&) const;
	void as1_GFp(gf_element&) const;
	bool iseq_GFp(const gf_element&, const gf_element&) const;
	bool is0_GFp(const gf_element&) const;
	bool is1_GFp(const gf_element&) const;
	void neg_GFp(gf_element&) const;
	void add_GFp(gf_element&, const gf_element&, const gf_element&) const;
	void sub_GFp(gf_element&, const gf_element&, const gf_element&) const;
	void mul_GFp(gf_element&, const gf_element&, const gf_element&) const;
	void inv_GFp(gf_element&) const;
	void sqr_GFp(gf_element&, const gf_element&) const;
	udigit hash_GFp(const gf_element&) const;
	const Fp_polynomial& get_pol_rep_GFp(const gf_element&) const;
	void set_pol_rep_GFp(gf_element&, const Fp_polynomial&) const;
	//void promote_GFp(gf_element&, const bigint&) const;
	//-----------------------------------------------------------------------------
	// functions for finite fields of odd characteristic and degree > 1
	//-----------------------------------------------------------------------------
	void construct_GFpn(gf_element&) const;
	void copy_GFpn(gf_element&, const gf_element&) const;
	void destruct_GFpn(gf_element&) const;
	void as0_GFpn(gf_element&) const;
	void as1_GFpn(gf_element&) const;
	bool iseq_GFpn(const gf_element&, const gf_element&) const;
	bool is0_GFpn(const gf_element&) const;
	bool is1_GFpn(const gf_element&) const;
	void neg_GFpn(gf_element&) const;
	void add_GFpn(gf_element&, const gf_element&, const gf_element&) const;
	void sub_GFpn(gf_element&, const gf_element&, const gf_element&) const;
	void mul_GFpn(gf_element&, const gf_element&, const gf_element&) const;
	void inv_GFpn(gf_element&) const;
	void sqr_GFpn(gf_element&, const gf_element&) const;
	udigit hash_GFpn(const gf_element&) const;
	const Fp_polynomial& get_pol_rep_GFpn(const gf_element&) const;
	void set_pol_rep_GFpn(gf_element&, const Fp_polynomial&) const;
	//void promote_GFpn(gf_element&, const bigint&) const;
	//-----------------------------------------------------------------------------
	// functions for finite fields of characteristic 2
	//-----------------------------------------------------------------------------
	void construct_GF2n(gf_element&) const;
	void copy_GF2n(gf_element&, const gf_element&) const;
	void destruct_GF2n(gf_element&) const;
	void as0_GF2n(gf_element&) const;
	void as1_GF2n(gf_element&) const;
	bool iseq_GF2n(const gf_element&, const gf_element&) const;
	bool is0_GF2n(const gf_element&) const;
	bool is1_GF2n(const gf_element&) const;
	void neg_GF2n(gf_element&) const;
	void add_GF2n(gf_element&, const gf_element&, const gf_element&) const;
	void sub_GF2n(gf_element&, const gf_element&, const gf_element&) const;
	void mul_GF2n(gf_element&, const gf_element&, const gf_element&) const;
	void inv_GF2n(gf_element&) const;
	void sqr_GF2n(gf_element&, const gf_element&) const;
	udigit hash_GF2n(const gf_element&) const;
	const Fp_polynomial& get_pol_rep_GF2n(const gf_element&) const;
	void set_pol_rep_GF2n(gf_element&, const Fp_polynomial&) const;
	//void promote_GF2n(gf_element&, const bigint&) const;
	//-----------------------------------------------------------------------------

	info_gf2n I2;
	Fp_poly_modulus poly_mod;

	//-----------------------------------------------------------------------------

private:
	//
	// internal representation of a Galois Field
	//

	bigint                                p;
	unsigned int                          n;
	bigint                                p_pow_n;
	mutable rational_factorization        p_pow_n_minus_1;
        mutable gf_element const*             gen;

	static galois_field_rep *head; // this is a linked list
	galois_field_rep *next; // of galois_field_rep's

	static galois_field_rep* search(galois_field_rep*, const bigint &,
					unsigned int);

	void init_p_pow_n_minus_1() const;
	void init_fact(const rational_factorization&);
	void init_pol(const Fp_polynomial&);
	void init_functions();

private:
	galois_field_rep(const galois_field_rep&); //disable
	galois_field_rep& operator = (const galois_field_rep&); //disable
	galois_field_rep(const bigint& characteristic,
                         unsigned int degree,
                         const Fp_polynomial& pol,
                         const rational_factorization& fact);

public:
	// use the next two functions instead of constructors:
	static galois_field_rep* create(const bigint & characteristic,
					unsigned int degree,
					const rational_factorization & fact = rational_factorization());
	static galois_field_rep* create(const Fp_polynomial & pol,
					const rational_factorization & fact = rational_factorization());
	~galois_field_rep();

	const bigint& characteristic() const
	{
		return p;
	}
	unsigned int  degree() const
	{
		return n;
	}
	const bigint& number_of_elements() const
	{
		return p_pow_n;
	}
	const rational_factorization& factorization_of_mult_order() const;
	Fp_polynomial irred_polynomial() const;

        gf_element const& generator() const;

	//static void output(); // for debugging purposes
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GALOIS_FIELD_REP_H_GUARD_
