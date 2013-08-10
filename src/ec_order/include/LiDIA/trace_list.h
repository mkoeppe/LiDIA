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


// A trace_list is a list of trace_mod's. Note that only Atkin primes
// are stored in the list, Elkies primes are handled directly. Therefore
// (C3, M3) holds the collection of all Elkies primes.
//
//   template class: used for eco_prime (bigmod) and eco_gf2n (gf2n)


#ifndef LIDIA_TRACE_LIST_H_GUARD_
#define LIDIA_TRACE_LIST_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_LIDIA_VECTOR_H_GUARD_
# include	"LiDIA/lidia_vector.h"
#endif
#ifndef LIDIA_TRACE_MOD_H_GUARD_
# include	"LiDIA/trace_mod.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curve.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include	"LiDIA/point.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class trace_list
{
public:

	friend std::ostream & operator << (std::ostream & o, const trace_list & t);
	friend std::istream & operator >> (std::istream & i, trace_list & t);

private:

	elliptic_curve< gf_element > E;

	sort_vector< trace_mod > l; // l[i] are the possible values of the trace c
	// modulo the i-th Atkin prime
	bigint M1; 		 // babystep modulus
	bigint M2; // giantstep modulus
	bigint M3; 	         // c = c3 mod M3
	bigint four_sqrt_q; // size of interval
	bigint C3; // trace c mod M3
	lidia_size_t last_index; // use l[0, ..., last_index] in BG phase

	static int info;
	static int MAX_NOF_TRACES; // maximal number of candidates for the trace

	void tl_sort (lidia_size_t left, lidia_size_t right);

public:

	trace_list  ();
	~trace_list ();

	const bigint& get_M1() const;
	const bigint& get_M2() const;
	const bigint& get_M3() const;
	const bigint& get_C3() const;
	bigint get_absolute_smallest_C3() const;
	const sort_vector< trace_mod > & get_list() const;

public:
	static void set_info_mode(int i=0);
	static void set_max_nof_traces(int m);

	void clear();
	void set_curve(const elliptic_curve< gf_element > & e);


private:
	bigint number_of_combinations();

	bool split_baby_giant(sort_vector< bigint > & baby,
			       sort_vector< bigint > & giant);

	void transform_lists(sort_vector< bigint > & baby,
			      sort_vector< bigint > & giant);

public:
	bool append(const trace_mod &); // add a trace mod to the list, return
	// true iff BG can start

	bigint bg_search_for_order();
	bigint simple_search_for_order();
	bool baby_giant_lists_correct (const bigint & ec_order);


};

// friend functions of class trace_list

std::ostream & operator << (std::ostream & o, const trace_list & t);
std::istream & operator >> (std::istream & i, trace_list & t);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_TRACE_LIST_H_GUARD_
