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
//	Author	: Nigel Smart (NS)
//                Adaption of John Cremona's code

//	Changes	:
//
//==============================================================================================


// This class implements  y^2 = quartic


#ifndef LIDIA_QUARTIC_H_GUARD_
#define LIDIA_QUARTIC_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include	"LiDIA/bigcomplex.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quartic
{
	bigint a, b, c, d, e; // coefficients
	bigcomplex* roots; // roots, array of 4 created in all constructors
	bigint ii, jj, disc;
	int type; // 1, 2 or 3

	// returns -1 for insoluble, 0 for undecided, +1 for soluble --- odd p
	int lemma6(const bigint& p, int nu, const bigint& x); // BSD-Lemma6
	// returns -1 for insoluble, 0 for undecided, +1 for soluble --- p = 2
	int lemma7(const bigint& p, int nu, const bigint& x); // BSD-Lemma7
	// Siksek's poly time test
	int samir_qp_soluble(const bigint& p);

public:

	// constructors
	quartic();
	quartic(const bigint& qa, const bigint& qb, const bigint& qc,
		const bigint& qd, const bigint& qe);
	~quartic()
	{
		delete[] roots;
	}

	quartic(const quartic& q);
	void assign(const bigint& qa, const bigint& qb, const bigint& qc,
		    const bigint& qd, const bigint& qe);
	void assign(const quartic& q);
	quartic & operator = (const quartic& q);

	// Access Functions
	bigint get_a() const
	{
		return a;
	}
	bigint get_b() const
	{
		return b;
	}
	bigint get_c() const
	{
		return c;
	}
	bigint get_d() const
	{
		return d;
	}
	bigint get_e() const
	{
		return e;
	}
	void get_coeffs(bigint& xa, bigint& xb, bigint& xc,
			bigint& xd, bigint& xe) const
	{
		xa = a;
		xb = b;
		xc = c;
		xd = d;
		xe = e;
	}
	int get_type() const
	{
		return type;
	}
	bigint get_I() const
	{
		return ii;
	}
	bigint get_J() const
	{
		return jj;
	}
	bigint get_H() const
	{
		return 8*a*c - 3*b*b;
	}
	bigint get_discriminant() const
	{
		return disc;
	}
	bigcomplex* get_roots() const
	{
		return roots;
	}

	// Advanced Stuff
	void doubleup()
	{
		b *= 2;
		c *= 4;
		d *= 8;
		e *= 16;
		ii *= 16;
		jj *= 64;
		disc *= 4096;
	}
	int trivial() const; // Checks for a rational root
	long no_roots_mod(long p) const;
	friend std::ostream& operator << (std::ostream& s, const quartic& q);

	// Checks for solublility in Zp with x = x0 (mod p^nu)
	int zp_soluble(const bigint& p, const bigint& x0, long nu);
	// Check Qp solublity
	int qp_soluble(const bigint& p);
	// Checks solubility at infinity and at plist (assumed to be a list of primes)
	int locally_soluble(const sort_vector< bigint > & plist);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUARTIC_H_GUARD_
