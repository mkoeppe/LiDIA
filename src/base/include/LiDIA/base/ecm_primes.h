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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ECM_PRIMES_H_GUARD_
#define LIDIA_ECM_PRIMES_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class ecm_primes;
class ec_point_M;
class ec_point_W;

class ec_curve
{
	friend class rational_factorization;
#if defined(_MSC_VER)
	friend template single_factor< bigint >;
#else
	friend class single_factor< bigint >;
#endif
	friend class ec_point_M;
	friend class ec_point_W;
	friend void mecgen16(ec_curve &, ec_point_M &, bigint &, long);
	friend bigint trans(ec_curve &, ec_point_W &, ec_point_M &);
	friend void mecfind(ec_curve &, ec_point_M &, bigint &, unsigned long,
			    unsigned long, ecm_primes &);
	friend void add (ec_point_W &, const ec_point_W &, const ec_point_W &,
			 bigint &);
	friend void multiply_by_2 (ec_point_W &, const ec_point_W &, bigint &);



private:

	bigmod a;
	bigmod b;



	//
	// c'tors and d'tor
	//

	ec_curve();
	ec_curve(const ec_curve & C);
	ec_curve(const bigmod & l_a, const bigmod & l_b);
	~ec_curve();



	//
	// accessors
	//

	const bigmod & A() const;
	const bigmod & B() const;



	//
	// assigners
	//

	void assign(const ec_curve & C);
	void assign(const bigmod & l_a, const bigmod & l_b);

	ec_curve & operator = (const ec_curve & C);

};



inline
ec_curve::ec_curve()
{
	// nothing to do
}



inline
ec_curve::ec_curve(const ec_curve & C)
	: a(C.a),
	  b(C.b)
{
	// nothing to do
}



inline
ec_curve::ec_curve(const bigmod & l_a, const bigmod & l_b)
	: a(l_a),
	  b(l_b)
{
	// nothing to do
}



inline
ec_curve::~ec_curve()
{
	// nothing to do
}



inline const bigmod &
ec_curve::A() const
{
	return a;
}



inline const bigmod &
ec_curve::B() const
{
	return b;
}



inline void
ec_curve::assign(const bigmod & l_a, const bigmod & l_b)
{
	a.assign(l_a);
	b.assign(l_b);
}



inline void
ec_curve::assign(const ec_curve & C)
{
	if (&C != this) {
		a.assign(C.a);
		b.assign(C.b);
	}
}



inline ec_curve &
ec_curve::operator = (const ec_curve & C)
{
	assign(C);
	return *this;
}






class ec_point_M
{
	friend class ec_curve;
	friend class rational_factorization;
#if defined(_MSC_VER)
	friend template single_factor< bigint >;
#else
	friend class single_factor< bigint >;
#endif
	friend void mecfind(ec_curve &, ec_point_M &, bigint &, unsigned long,
			    unsigned long, ecm_primes &);
	friend void mecgen16(ec_curve &, ec_point_M &, bigint &, long);
	friend void add(ec_point_M &, const ec_point_M &, const ec_point_M &, const
			ec_point_M &);
	friend void multiply_by_2(ec_point_M &, const ec_point_M &, const bigmod &);
	friend void multiply(ec_point_M &, const ec_point_M &, const bigmod &,
			     unsigned long, unsigned long);
	friend bigint trans(ec_curve &, ec_point_W &, ec_point_M &);


private:

	bigmod x;
	bigmod z;


	//
	// c'tors and d'tor
	//

	ec_point_M();
	ec_point_M(const ec_point_M & mp);
	~ec_point_M();



	//
	// accessors
	//

	const bigmod & X() const;
	const bigmod & Z() const;



	//
	// assigners
	//

	void assign(const ec_point_M & mp);
	void assign(const bigmod & l_x, const bigmod & l_y);

	ec_point_M & operator = (const ec_point_M & mp);

};



inline 
ec_point_M::ec_point_M()
{
	// nothing to do
}



inline
ec_point_M::ec_point_M(const ec_point_M & mp)
	: x(mp.x),
	  z(mp.z)
{
	// nothing to do
}



inline
ec_point_M::~ec_point_M()
{
	// nothing to do
}



inline const bigmod &
ec_point_M::X() const
{
	return x;
}



inline const bigmod &
ec_point_M::Z() const
{
	return z;
}



inline void
ec_point_M::assign(const ec_point_M & mp)
{
	if (&mp != this) {
		x.assign(mp.x);
		z.assign(mp.z);
	}
}



inline void
ec_point_M::assign(const bigmod & l_x, const bigmod & l_y)
{
	x.assign(l_x);
	z.assign(l_y);
}



inline ec_point_M &
ec_point_M::operator = (const ec_point_M & mp)
{
	assign(mp);
	return *this;
}






class ec_point_W
{
	friend class rational_factorization;
	friend class ec_curve;
#if defined(_MSC_VER)
	friend template single_factor< bigint >;
#else
	friend class single_factor< bigint >;
#endif
	friend void cont(ec_point_W &, bigint &, unsigned long, unsigned long,
			 ecm_primes &);
	friend void mecgen16(ec_curve &, ec_point_M &, bigint &, long);
	friend bigint trans(ec_curve &, ec_point_W &, ec_point_M &);

private:


	static const bigmod	zero;



	bigmod x;
	bigmod y;

	const ec_curve  * curve;

	bool is_0; //Flag if the point is infinite or not



	//
	// c'tors and d'tor
	//

	ec_point_W();
	ec_point_W(const ec_curve & we);
	ec_point_W(const ec_curve & we, const bigmod & l_x, const bigmod & l_y);
	ec_point_W(const ec_point_W & wp);
	~ec_point_W();



	//
	// accessors
	//

	const ec_curve * Curve() const;

	const bigmod & X() const;
	const bigmod & Y() const;



	//
	// assigners
	//

	void assign(const bigmod & l_x, const bigmod & l_y);
	void assign_one();
	void assign_zero();
	void assign(const ec_point_W & wp);

	ec_point_W & operator = (const ec_point_W & wp);


	//
	// predicates
	//

	bool is_negative(const ec_point_W & wp) const;



	friend bool operator == (const ec_point_W &, const ec_point_W &);
	friend void add(ec_point_W &, const ec_point_W &, const ec_point_W &,
			bigint &);
	friend void addPQ(ec_point_W &, const ec_point_W &, const ec_point_W &,
			  bigint &);
	friend void multiply_by_2(ec_point_W &, const ec_point_W &, bigint &);
	friend void multiply(ec_point_W &, const ec_point_W &, long, bigint &);

};



inline
ec_point_W::ec_point_W()
	: x(),
	  y(),
	  curve(NULL),
	  is_0(true)
{
	// nothing to do
}



inline
ec_point_W::ec_point_W(const ec_curve & we)
	: x(),
	  y(),
	  curve(&we),
	  is_0(true)
{
	// nothing to do
}



inline
ec_point_W::ec_point_W(const ec_curve & we, const bigmod & l_x, const bigmod & l_y)
	: x(l_x),
	  y(l_y),
	  curve(&we),
	  is_0(false)
{
	// nothing to do
}



inline
ec_point_W::ec_point_W(const ec_point_W & wp)
{
	if (is_0 != wp.is_0) {
		x = wp.x;
		y = wp.y;

		curve = wp.curve;
	}
}



inline
ec_point_W::~ec_point_W()
{
	// nothing to do
}



inline const ec_curve *
ec_point_W::Curve() const
{
	return curve;
}



inline const bigmod &
ec_point_W::X() const
{
	if (!is_0)
		return x;
	else {
		warning_handler("EC_POINT_W", "x-coordinate of infinite point");
		return zero;
	}
}



inline const bigmod &
ec_point_W::Y() const
{
	if (!is_0)
		return y;
	else {
		warning_handler("EC_POINT_W", "y-coordinate of infinite point");
		return zero;
	}
}



inline void
ec_point_W::assign(const bigmod & l_x, const bigmod & l_y)
{
	x.assign(l_x);
	y.assign(l_y);
	is_0 = false;
}



inline void
ec_point_W::assign_one()
{
	is_0 = false;
}



inline void
ec_point_W::assign_zero()
{
	is_0 = true;
}



inline void
ec_point_W::assign(const ec_point_W & wp)
{
	if (is_0 != wp.is_0) {
		x.assign(wp.x);
		y.assign(wp.y);
		curve = wp.curve;
	}
}



inline ec_point_W &
ec_point_W::operator = (const ec_point_W & wp)
{
	assign(wp);
	return *this;
}



inline bool
ec_point_W::is_negative(const ec_point_W & wp) const
{
	// added this to allow compiling on Mips with CC, TP
	bool i;
	bigmod h(wp.y);

	add(h, h, y);

	if ((is_0 == wp.is_0) && (is_0 || (x == wp.x) && h.is_zero()))
		i = true;
	else
		i = false;
	return i;
}






class ecm_primes
{
	friend class rational_factorization;
#if defined(_MSC_VER)
	friend template single_factor< bigint >;
#else
	friend class single_factor< bigint >;
#endif
	friend void mecfind(ec_curve &, ec_point_M &, bigint &, unsigned long,
			    unsigned long, ecm_primes &);
	friend void cont(ec_point_W &, bigint &, unsigned long, unsigned long,
			 ecm_primes &);

	friend void mpqs(lidia_size_t index, ecm_primes & prim);
	friend void mpqs_comp(lidia_size_t index);



private:

	unsigned int *table, // table to hold the actual primetable
		*table2; // table to hold the primes up to $sqrt$(max)
	unsigned int *over; // table to hold the overhead
	//
	unsigned int max, // upper limit for the primes
		p_max, // number of odd numbers < sqrt(max)
		last_prime, // last prime number in the actual primetable
		dim_ar, // length of the actual primetable
		space, // upper limit for memory

		feld, // pointers to hold the position in
		feld_abs, // the primetables
		pl;

	unsigned int  maske35, maske7, // patterns for sieving all
		maske11, maske13, // of the primes 3, 5, 7, 11, 13, 17,
		maske17, maske19, //               19, 23, 29, 31
		maske23, maske29,
		maske31, maske_i;

	char flag; // flag == 0 : table contains all primes
	// flag == 1 : we must sieve on several intervals



	void sieve_array(int);
	static const unsigned int  BIT;

	void initprimes(unsigned int, unsigned int, unsigned int);
	void killprimes();


public:
	ecm_primes(unsigned int us, unsigned int os, unsigned int sp)
	{
		initprimes(us, os, sp);
	}

	~ecm_primes()
	{
		killprimes();
	}

	unsigned int getprimes();
	void resetprimes(unsigned int);



	// inhibit:

private:

	ecm_primes();
	ecm_primes(const ecm_primes &);
	ecm_primes & operator = (const ecm_primes &);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ECM_PRIMES_H_GUARD_
