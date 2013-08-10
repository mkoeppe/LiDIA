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
//	Changes	: see CVS log
//
//==============================================================================================


#ifndef LIDIA_LATTICE_MODULES_H_GUARD_
#define LIDIA_LATTICE_MODULES_H_GUARD_


#ifndef LIDIA_LATTICE_DEFS_H_GUARD_
# include	"LiDIA/lattices/lattice_defs.h"
#endif
#ifndef LIDIA_P_VECTOR_H_GUARD_
# include	"LiDIA/lattices/p_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



const int Normal = 0;
const int VariationI = 1;
const int VariationII = 2;

const sdigit Magic_Precision_Factor = 2;



//
// class prec_modules < E, A >
// ~~~~~~~~~~~~~~~~~~~~~~~~
//
// a class  which cares about needed precision
// for both the template type E and the approximate type A
//
template< class E, class A >
class prec_modules
{
protected :
	sdigit tempsdt;
public :
	prec_modules() { };
	~prec_modules() { };

	void prec_startup(dense_alg< E > &);
	void prec_update(dense_alg< E > &);
	void prec_exact(const dense_alg< E > &);
	void prec_approx(const dense_alg< E > &);
	void prec_correct(const dense_alg< E > &);
	void prec_restore(const dense_alg< E > &);
};



//
// class base_modules < E, A, int >
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// a class  which implements parts of the lll -
// algorithm for both basis and gensys
//
template< class E, class A, const int V >
class base_modules : public prec_modules< E, A >
{
protected :
	sdigit   bi_bit_len;
	sdigit   bit_diff;
	bigint   mant;
	sdigit   expo;
	bigint   tempbin; // working on mantissa
	bigfloat tempbfl;

public :
	base_modules() : prec_modules< E, A > () { };
	~base_modules() { };

	bool A_convert_value_E_bound(dense_alg< E > &, E&, const A&, const A&);
	void E_convert_value_A(dense_alg< E > &, A&, const E&);
};



//
// class basis_modules < E, A, int >
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// a class  which implements specialized modules of the lll -
// algorithm for basis lattices
//
template< class E, class A, const int V >
class basis_modules : public base_modules< E, A, V >
{
protected :
	bool is_zero;
	void *tempP; // swapping rows

public :
	basis_modules() : base_modules< E, A, V > () { };
	~basis_modules() { };

	bool E_convert_vector_A(dense_alg< E > &, A**, const lidia_size_t);
	void E_convert_lattice_A(dense_alg< E > &, A**);
};



//
// class gensys_modules < E, A, int >
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// a class  which implements specialized modules of the lll -
// algorithm for gensys lattices
//
template< class E, class A, const int V >
class gensys_modules : public base_modules< E, A, V >
{
protected :
	bool is_zero;
	void *tempP; // swapping rows

public :
	gensys_modules() : base_modules< E, A, V > () { };
	~gensys_modules() { };

	bool E_convert_vector_A(dense_alg< E > &, A**, const lidia_size_t);
	void E_convert_lattice_A(dense_alg< E > &, A**);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#include	"LiDIA/lattices/lattice_modules.cc"



#endif	// LIDIA_LATTICE_MODULES_H_GUARD_
