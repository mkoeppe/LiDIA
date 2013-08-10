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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BF_LATTICE_BASIS_H_GUARD_
#define LIDIA_BF_LATTICE_BASIS_H_GUARD_


#ifndef LIDIA_BI_LATTICE_BASIS_H_GUARD_
# include	"LiDIA/lattices/bi_lattice_basis.h"
#endif
#ifndef LIDIA_BF_LATTICE_GENSYS_H_GUARD_
# include	"LiDIA/lattices/bf_lattice_gensys.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigfloat_lattice_basis:public bigfloat_lattice_gensys
{
public :

//
// Constructor / Destructors
//
	bigfloat_lattice_basis():bigfloat_lattice_gensys() {};
//    bigfloat_lattice_basis(lidia_size_t n):bigfloat_lattice_gensys(n) {};
	bigfloat_lattice_basis(lidia_size_t n, lidia_size_t m):bigfloat_lattice_gensys(n, m) {};
	bigfloat_lattice_basis(lidia_size_t n, lidia_size_t m, const bigfloat** abi):bigfloat_lattice_gensys(n, m, abi) {};
	bigfloat_lattice_basis(const math_matrix< bigfloat > & M):bigfloat_lattice_gensys(M) {};
	bigfloat_lattice_basis(const bigfloat_lattice_gensys& M):bigfloat_lattice_gensys(M) {};
	~bigfloat_lattice_basis() {};

//
// check - basis
//
	bool lll_check(double);
	bool lll_check(sdigit, sdigit);

	double lll_check_search();
	void lll_check_search(sdigit&, sdigit&);

//
// Conversion
//
	void extract_basis(const bigfloat_lattice_gensys&, lidia_size_t&);
	bool make_basis(const bigfloat_lattice_gensys&, lidia_size_t&, sdigit = 1);

//
// Tools
//
	void gram_schmidt_orth(math_matrix< bigfloat > &, math_matrix< bigfloat > &);

//
// Algorithms
//

//
// Schnorr - Euchner using doubles for parameter = 1
// and x*doubleprec for parameter = x
//
	void lll(sdigit = 1);
	void lll(const math_matrix< bigfloat > &, sdigit = 1);
        friend bigfloat_lattice_basis lll(const bigfloat_lattice_basis& L, sdigit factor);

	void lll_var(sdigit = 1);
	void lll_var(const math_matrix< bigfloat > &, sdigit = 1);
        friend bigfloat_lattice_basis lll_var(const bigfloat_lattice_basis& L, sdigit factor);

//    void lll_deep_insert(sdigit = 5); // not yet implemented
//    void lll_deep_insert(const math_matrix < bigfloat > &, sdigit = 5); // not yet implemented

	void lll_trans(math_matrix< bigint > &, sdigit = 1);
	void lll_trans(math_matrix< bigfloat > &, sdigit = 1);
	friend bigfloat_lattice_basis lll_trans(const bigfloat_lattice_basis& L,
						math_matrix< bigint > & T,
						sdigit factor);

	friend bigfloat_lattice_basis lll_trans(const bigfloat_lattice_basis& L,
						math_matrix< bigfloat > & T,
						sdigit factor);

	void lll_trans_var(math_matrix< bigint > &, sdigit = 1);
	void lll_trans_var(math_matrix< bigfloat > &, sdigit = 1);
	friend bigfloat_lattice_basis lll_trans_var(const bigfloat_lattice_basis& L,
						    math_matrix< bigint > & T,
						    sdigit factor);

	friend bigfloat_lattice_basis lll_trans_var(const bigfloat_lattice_basis& L,
						    math_matrix< bigfloat > & T,
						    sdigit factor);

protected :

//
// protected functions needed by algorithms
//
	bool Tr_gram_schmidt_orth(bigfloat**, bigfloat**);

	bool Tr_make_basis(lidia_size_t&, sdigit = 1);
	void Tr_extract_basis(lidia_size_t&);

	bool Tr_lll_check(double);
	void Tr_lll_check_search(sdigit&, sdigit&);

//
// Dimension Checking
//
	bool Tr_check_basis()
	{
		if (trans_flag == true) {
			if (rows < columns)
				return(false);
			else
				return(true);
		}
		else {
			if (columns < rows)
				return(false);
			else
				return(true);
		}
	}
};


// friend functions

inline
bigfloat_lattice_basis lll(const bigfloat_lattice_basis& L, sdigit factor = 1)
{
	bigfloat_lattice_basis TL(L);
	TL.lll(factor);
	return (TL);
}

inline
bigfloat_lattice_basis lll_var(const bigfloat_lattice_basis& L, sdigit factor = 1)
{
	bigfloat_lattice_basis TL(L);
	TL.lll_var(factor);
	return (TL);
}

inline
bigfloat_lattice_basis lll_trans(const bigfloat_lattice_basis& L,
				     math_matrix< bigint > & T,
				     sdigit factor = 1)
{
	bigfloat_lattice_basis LT(L);
	LT.lll_trans(T, factor);
	return (LT);
}

inline
bigfloat_lattice_basis lll_trans(const bigfloat_lattice_basis& L,
				     math_matrix< bigfloat > & T,
				     sdigit factor = 1)
{
	bigfloat_lattice_basis LT(L);
	LT.lll_trans(T, factor);
	return (LT);
}

inline
bigfloat_lattice_basis lll_trans_var(const bigfloat_lattice_basis& L,
					 math_matrix< bigint > & T,
					 sdigit factor = 1)
{
	bigfloat_lattice_basis LT(L);
	LT.lll_trans_var(T, factor);
	return (LT);
}

inline
bigfloat_lattice_basis lll_trans_var(const bigfloat_lattice_basis& L,
					 math_matrix< bigfloat > & T,
					 sdigit factor = 1)
{
	bigfloat_lattice_basis LT(L);
	LT.lll_trans_var(T, factor);
	return (LT);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BF_LATTICE_BASIS_H_GUARD_
