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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MODULAR_ARITHMETIC_H_GUARD_
#define LIDIA_MODULAR_ARITHMETIC_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class REP, class SINGLE_MODUL, class MULTI_MODUL >
class modular_arithmetic
{

public:

	const REP rep_modul;
	const SINGLE_MODUL long_modul;
	const MULTI_MODUL bigint_modul;

	//
	// constructor
	//

	modular_arithmetic() {}

	//
	// destructor
	//

	~modular_arithmetic() {}

	//
	// Chinese remaindering theorem
	//

	void chinrest(matrix< bigint > &, const matrix< bigint > *, const bigint *) const;
	void chinese_remainder(matrix< bigint > &, const bigint &,
			       const matrix< bigint > &, const bigint &) const;

	//
	// rank
	//

	lidia_size_t rank(const matrix< bigint > &, const bigint &) const;

	//
	// rank and linearly independent rows
	//

	lidia_size_t * lininr1(const matrix< bigint > &, const bigint &) const;
	lidia_size_t * lininr2(const matrix< bigint > &, const bigint &) const;

	//
	// rank linearly independent columns
	//

	lidia_size_t * lininc1(const matrix< bigint > &, const bigint &) const;
	lidia_size_t * lininc2(const matrix< bigint > &, const bigint &) const;

	//
	// adjoint matrix
	//

	void adj1(matrix< bigint > &, const matrix< bigint > &, const bigint &, const bigint &) const;
	void adj2(matrix< bigint > &, const matrix< bigint > &, const bigint &, const bigint &) const;
	void adj2(matrix< bigint > &, const matrix< bigint > &, const bigint &, const bigint &, int) const;

	//
	// lattice determinant
	//

	void latticedet1(const matrix< bigint > &, bigint &, const bigint &) const;
	void latticedet2(const matrix< bigint > &, bigint &, const bigint &) const;
	void latticedet3(const matrix< bigint > &, bigint &, const bigint &) const;

	void latticedet2(const matrix< bigint > &, bigint &, bigint &, int) const;
	void latticedet4(const matrix< bigint > &, bigint &, bigint &, int) const;
	lidia_size_t *latticedet5(matrix< bigint > &, bigint &, bigint &, int) const;

	void latticedet_special(const matrix< bigint > &, bigint &, const bigint &, int) const;

	//
	// determinant
	//

	void det(const matrix< bigint > &, bigint &, const bigint &) const;
	int det(const matrix< bigint > &RES, bigint & DET, const bigint &H, int num_same) const;

	//
	// characteristic polynomial
	//

	void charpoly(const matrix< bigint > &, bigint *, const bigint &) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/modular_arithmetic.cc"
#endif



#endif	// LIDIA_MODULAR_ARITHMETIC_H_GUARD_
