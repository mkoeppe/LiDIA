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
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


// The class ff2 realizes single precision arithmetic in a quadratic
// field extension GF(l^2).


#ifndef LIDIA_FF2_H_GUARD_
#define LIDIA_FF2_H_GUARD_



#ifndef LIDIA_FF1_H_GUARD_
# include	"LiDIA/ff1.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class ff2
{
private :

	ff1 a; // gen. pol. of field is X^2 - tau, tau non square
	ff1 b; // element is represented as b*X + a
	static ff1 tau;

public :

	ff2 ();
	ff2 (const ff2 &);
	ff2 (const ff1 &);
	ff2 (const ff1 & aa, const ff1 & bb);
	~ff2 ();

	static void  init_field (const udigit &);
	static udigit characteristic ();
	ff2 & operator = (const ff2 &);
	void assign (const ff2 &);
	void assign (const ff1 &aa, const ff1 & bb);
	void assign_one();
	void assign_zero();
	bool is_one() const;
	bool is_zero() const;

	bool as_udigit(udigit & x)
	{
		if(!b.is_zero()) {
			lidia_error_handler("ff2", "as_udigit::element not in prime field");
			return false;
		}
		else {
			x = a.as_udigit();
			return true;
		}
	}

	friend std::istream &  operator >> (std::istream &, ff2 &);
	friend std::ostream &  operator << (std::ostream &, const ff2 &);
	friend void swap (ff2 &, ff2 &);
	friend void add      (ff2 &, const ff2 &, const ff2 &);
	friend void subtract (ff2 &, const ff2 &, const ff2 &);
	friend void multiply (ff2 &, const ff2 &, const ff2 &);
	friend void divide   (ff2 &, const ff2 &, const ff2 &);
	friend void invert   (ff2 &, const ff2 &);
	friend void negate   (ff2 &, const ff2 &);
	friend void power (ff2 &, const ff2 &, long);

	udigit multiplicative_order () const;
	friend void nearly_all_of_order (udigit, base_vector< ff2 >&);
};

// friend functions of ff2

std::istream &  operator >> (std::istream &, ff2 &);
std::ostream &  operator << (std::ostream &, const ff2 &);
void swap (ff2 &, ff2 &);
void add      (ff2 &, const ff2 &, const ff2 &);
void subtract (ff2 &, const ff2 &, const ff2 &);
void multiply (ff2 &, const ff2 &, const ff2 &);
void divide   (ff2 &, const ff2 &, const ff2 &);
void invert   (ff2 &, const ff2 &);
void negate   (ff2 &, const ff2 &);
void power (ff2 &, const ff2 &, long);

void nearly_all_of_order (udigit, base_vector< ff2 >&);


inline void square  (ff2 & r, const ff2 & x)
{
	multiply(r, x, x);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FF2_H_GUARD_
