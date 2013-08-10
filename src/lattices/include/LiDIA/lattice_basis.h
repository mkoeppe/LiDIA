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


#ifndef LIDIA_LATTICE_BASIS_H_GUARD_
#define LIDIA_LATTICE_BASIS_H_GUARD_


#ifndef HEADBANGER

#ifndef LIDIA_LATTICE_GENSYS_H_GUARD_
# include	"LiDIA/lattice_gensys.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class lattice_basis: public lattice_gensys
{
public :

//
// Constructors / Destructor
//

	lattice_basis(lidia_size_t n, lidia_size_t m):lattice_gensys(n, m) { };
	lattice_basis(lidia_size_t n, lidia_size_t m, const double** ado):lattice_gensys(n, m, ado) { };
	lattice_basis(lidia_size_t n, lidia_size_t m, const bigint** abi):lattice_gensys(n, m, abi) { };
	lattice_basis(lidia_size_t n, lidia_size_t m, const bigfloat** abf):lattice_gensys(n, m, abf) { };
	lattice_basis(const lattice_basis& L):lattice_gensys(L) { };
	lattice_basis(const math_matrix< bigint > & M):lattice_gensys(M) { };
	lattice_basis(const math_matrix< bigfloat > & M):lattice_gensys(M) { };
	virtual ~lattice_basis() { };

//
// Conversion
//

	void extract_basis(const lattice_gensys&, lidia_size_t&);
	bool make_basis(const lattice_gensys&, lidia_size_t&);

//
// Algorithms
//
	void lll();
	void lll(const lattice_basis&);
        friend lattice_basis lll(const lattice_basis& L);

	void lll_trans(lattice_basis&);
        friend lattice_basis lll_trans(const lattice_basis& L,
				       lattice_basis& T);

	void gram_schmidt_orth(lattice_basis&, lattice_basis&);

	bool lll_check(double);
	bool lll_check(sdigit, sdigit);

	double lll_check_search();
	void lll_check_search(sdigit&, sdigit&);

	friend void lll(math_matrix< bigint > &, const math_matrix< bigint > &);
	friend void lll_trans(math_matrix< bigint > &, math_matrix< bigint > &);

	friend void lll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &);
	friend void lll_trans(math_matrix< bigint > &, math_matrix< bigfloat > &);
	friend void lll_trans(math_matrix< bigfloat > &, math_matrix< bigfloat > &);
};

// friend functions

inline
lattice_basis lll(const lattice_basis& L)
{
	lattice_basis TL(L);
	TL.lll();
	return (TL);
}

inline
lattice_basis lll_trans(const lattice_basis& L, lattice_basis& T)
{
	lattice_basis LT(L);
	LT.lll_trans(T);
	return (LT);
}

void lll(math_matrix< bigint > &, const math_matrix< bigint > &);
void lll_trans(math_matrix< bigint > &, math_matrix< bigint > &);

void lll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &);
void lll_trans(math_matrix< bigint > &, math_matrix< bigfloat > &);
void lll_trans(math_matrix< bigfloat > &, math_matrix< bigfloat > &);


#endif  // HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LATTICE_BASIS_H_GUARD_
