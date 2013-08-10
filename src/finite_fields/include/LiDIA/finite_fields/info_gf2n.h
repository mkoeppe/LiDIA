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
//	Author	: Franz-Dieter Berger (FDB), Patric Kirsch (PK),
//                Volker Mueller (VM), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_INFO_GF2N_H_GUARD_
#define LIDIA_INFO_GF2N_H_GUARD_


#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class galois_field_rep;


class info_gf2n
{
	friend class galois_field_rep;

	static udigit *B, *C, *F, *G; // temp. storage
	static unsigned int max_anzBI; // max. anzBI of all info_gf2n
	static void gen_tables();


	void mul1(udigit*, udigit*, udigit*) const;
	void mul2(udigit*, udigit*, udigit*) const;
	void mul3(udigit*, udigit*, udigit*) const;
	void mul4(udigit*, udigit*, udigit*) const;
	void mul5(udigit*, udigit*, udigit*) const;
	void mul6(udigit*, udigit*, udigit*) const;
	void mul7(udigit*, udigit*, udigit*) const;
	void mul8(udigit*, udigit*, udigit*) const;
	void karatsuba_mul(udigit*, udigit*, udigit*) const;
public:
	enum {
		FIXMUL = 9	// number of multiplication routines
	};
	static void (info_gf2n::*gf2nmul[]) (udigit*, udigit*, udigit*) const;

private:
	void tri_invert(udigit*, udigit*) const;
	void pent_invert(udigit*, udigit*) const;
	void general_invert(udigit*, udigit*) const;

	void tri_partial_reduce1(udigit*) const;
	void pent_partial_reduce1(udigit*) const;
	void general_partial_reduce1(udigit*) const;

	void tri_partial_reduce2(udigit*) const;
	void pent_partial_reduce2(udigit*) const;
	void general_partial_reduce2(udigit*) const;


	unsigned int degree;
	unsigned int anzBI; // number of udigits for each gf2n
	unsigned int *exponents; // array of exponents with 1 coefficients of
    				// generating polynomial (excl. degree and 0)
	unsigned int exp1, exp2, exp3; // for storing the exponents of 3/5-nomials
	unsigned int anz_exponents; // number of nonzero exponents


	void (info_gf2n::*mul_p) (udigit*, udigit*, udigit*) const;
	void (info_gf2n::*invert_p)(udigit*, udigit*) const;
	void (info_gf2n::*partial_reduce1_p) (udigit *) const;
	void (info_gf2n::*partial_reduce2_p) (udigit *) const;


public:

	info_gf2n();
	void init(unsigned int);
	void init(char*, unsigned int);

	void mul(udigit*c, udigit* a, udigit* b) const
	{
		(this->*mul_p)(c, a, b);
	}
	void inv(udigit* a, udigit* b) const
	{
		(this->*invert_p)(a, b);
	}
	void partial_reduce1(udigit* a) const
	{
		(this->*partial_reduce1_p)(a);
	}
	void partial_reduce2(udigit* a) const
	{
		(this->*partial_reduce2_p)(a);
	}

	void square(udigit*, udigit*) const;

	bool is_reduced(const udigit* a) const
	{
		return ((a[anzBI-1] >> (degree % BITS_PER_LONG)) == static_cast<udigit>(0));
	}

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_INFO_GF2N_H_GUARD_
