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


#ifndef LIDIA_LATTICE_KERNEL_H_GUARD_
#define LIDIA_LATTICE_KERNEL_H_GUARD_


#ifndef LIDIA_XDOUBLE_H_GUARD_
# include	"LiDIA/xdouble.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif
#ifndef LIDIA_P_VECTOR_H_GUARD_
# include	"LiDIA/lattices/p_vector.h"
#endif
#ifndef LIDIA_LATTICE_DEFS_H_GUARD_
# include	"LiDIA/lattices/lattice_defs.h"
#endif
#ifndef LIDIA_LATTICE_MODULES_H_GUARD_
# include	"LiDIA/lattices/lattice_modules.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint_lattice;
class bigfloat_lattice;



//
// class lll_kernel_{op, fu}
// ~~~~~~~~~~~~~~~~~~~~~~~~
//
// These classes feature the real implementation of the
// variation of Schnorr - Euchner - lll we support for
// the class bigint/bigfloat_lattice. These algorithms are
// working row by row.
//
// Warning :
// ~~~~~~~~~
// These classes is for internal use only !!!
// There are many things you have to do before running
// an algorithm of these classes ! These things are done
// by the class bigint_lattice !
// For internal use only !
//
// This class is working with operators
//
template< class E, class A, class VectorOp, class AlgParts >
class lll_kernel_op
{
protected :
	VectorOp vector;
	AlgParts alg;

public :
	lll_kernel_op() { }
	~lll_kernel_op() { }

	VectorOp& vector_operations()
	{
		return (vector);
	}

//
// Dense Lattice :
// ~~~~~~~~~~~~~~~
// Schnorr - Euchner - lll algorithms using a datatype
// with fixed precision for approximation of the lattice
// The computation is done with implemented operators
// (like the double datatype)
// The datatype has to be functional compatible to
// the double datatype like the new xdouble datatype
// +, -, *, /, = , >, < , sqrt, fabs needed
//
	void lll(dense_alg< E > &, lattice_info&);
	void lll_var(dense_alg< E > &, lattice_info&);
//
// Versions modified for working with a Gram - matrix
//
	void lll_gram(dense_alg< E > &, lattice_info&);
	void lll_var_gram(dense_alg< E > &, lattice_info&);
};

//
// This class is working with function
//
template< class E, class A, class VectorOp, class AlgParts >
class lll_kernel_fu
{
protected :
	VectorOp vector;
	AlgParts alg;

public :
	lll_kernel_fu() { }
	~lll_kernel_fu() { }

	VectorOp& vector_operations()
	{
		return (vector);
	}

//
// Dense Lattice :
// ~~~~~~~~~~~~~~~
// Schnorr - Euchner - lll algorithms using a datatype
// with variable precision for approximation of the lattice
// The computation is done with implemented functions
// (like the bigfloat datatype)
//
	void lll(dense_alg< E > &, lattice_info&);
	void lll_var(dense_alg< E > &, lattice_info&);
//
// Versions modified for working with a Gram - matrix
//
	void lll_gram(dense_alg< E > &, lattice_info&);
	void lll_var_gram(dense_alg< E > &, lattice_info&);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/lattices/lattice_kernel.cc"
#endif



#endif	// LIDIA_LATTICE_KERNEL_H_GUARD_
