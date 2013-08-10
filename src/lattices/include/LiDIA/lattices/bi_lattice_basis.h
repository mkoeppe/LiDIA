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


#ifndef LIDIA_BI_LATTICE_BASIS_H_GUARD_
#define LIDIA_BI_LATTICE_BASIS_H_GUARD_


#ifndef LIDIA_BI_LATTICE_GENSYS_H_GUARD_
# include	"LiDIA/lattices/bi_lattice_gensys.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint_lattice_basis:public bigint_lattice_gensys
{
protected :
	friend          class bigfloat_lattice_gensys;
	friend          class bigfloat_lattice_basis;

public :

//
// Constructor / Destructors
//
	bigint_lattice_basis():bigint_lattice_gensys() {};
//    bigint_lattice_basis(lidia_size_t n):bigint_lattice_gensys(n) {};
	bigint_lattice_basis(lidia_size_t n, lidia_size_t m):bigint_lattice_gensys(n, m) {};
	bigint_lattice_basis(lidia_size_t n, lidia_size_t m, const bigint** abi):bigint_lattice_gensys(n, m, abi) {};
	bigint_lattice_basis(const math_matrix< bigint > & M):bigint_lattice_gensys(M) {};
	bigint_lattice_basis(const bigint_lattice_gensys& M):bigint_lattice_gensys(M) {};
	~bigint_lattice_basis() {};

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
	void extract_basis(const bigint_lattice_gensys&, lidia_size_t&);
	bool make_basis(const bigint_lattice_gensys&, lidia_size_t&, sdigit = 1);

//
// Tools
//
	void gram_schmidt_orth(math_matrix< bigfloat > &, math_matrix< bigfloat > &);

//
// Algorithms
//

//
// Schnorr - Euchner using doubles for parameter = 1
// and x*doubleprecision for parameter = x
//
	void lll(sdigit = 1);
	void lll(const math_matrix< bigint > &, sdigit = 1);
        friend bigint_lattice_basis lll(const bigint_lattice_basis& L,
					sdigit factor);

//    void lll_deep_insert(sdigit = 5); // not yet implemented
//    void lll_deep_insert(const math_matrix < bigint > &, sdigit = 5); // not yet implemented

	void lll_trans(math_matrix< bigint > &, sdigit = 1);
	friend bigint_lattice_basis lll_trans(const bigint_lattice_basis& L,
					      math_matrix< bigint > & T,
					      sdigit factor);

//
// Benne de Weger only using bigints
//
	void lll_benne_de_weger();
	void lll_benne_de_weger(const math_matrix< bigint > &);
        friend bigint_lattice_basis lll_benne_de_weger(const bigint_lattice_basis& L);

//    void lll_rand_benne_de_weger(); //   Warning :
//    void lll_rand_benne_de_weger(const math_matrix < bigint > &); //   Experimental Version

	void lll_trans_benne_de_weger(math_matrix< bigint > &);
	friend bigint_lattice_basis lll_trans_benne_de_weger(const bigint_lattice_basis& L,
							     math_matrix< bigint > & T);


protected :

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

//
// protected functions needed by algorithms
//
	bool Tr_gram_schmidt_orth_bdw(bigint_lattice_basis&, bigint_lattice_basis&, bigint*);
	bool Tr_gram_schmidt_orth(bigfloat**, bigfloat**);

	void Tr_extract_basis(lidia_size_t&);
	bool Tr_make_basis(lidia_size_t&, sdigit = 1);

	bool Tr_lll_check(double);
	void Tr_lll_check_search(sdigit&, sdigit&);
};


// friend functiuons

inline
bigint_lattice_basis lll(const bigint_lattice_basis& L, sdigit factor = 1)
{
	bigint_lattice_basis TL(L);
	TL.lll(factor);
	return (TL);
}

inline
bigint_lattice_basis lll_trans(const bigint_lattice_basis& L,
				   math_matrix< bigint > & T,
				   sdigit factor = 1)
{
	bigint_lattice_basis LT(L);
	LT.lll_trans(T, factor);
	return (LT);
}


//
// Benne de Weger only using bigints
//

inline
bigint_lattice_basis lll_benne_de_weger(const bigint_lattice_basis& L)
{
	bigint_lattice_basis TL(L);
	TL.lll_benne_de_weger();
	return (TL);
}

inline
bigint_lattice_basis lll_trans_benne_de_weger(const bigint_lattice_basis& L,
						  math_matrix< bigint > & T)
{
	bigint_lattice_basis LT(L);
	LT.lll_trans_benne_de_weger(T);
	return (LT);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BI_LATTICE_BASIS_H_GUARD_
