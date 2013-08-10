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


//
//  usage:
//
//	include all headers that are necessary to describe the types E_TYPE and A_TYPE
//
//	define the types E_TYPE and A_TYPE
//
//	define VAR (Normal or VariationI)
//
//	define MODULES (basis_modules or gensys_modules)
//
//	define VECTOR (vector_op or vector_op_SP)
//
//	include this file
//


#include	"LiDIA/bigint_lattice.h"
#include	"LiDIA/bigfloat_lattice.h"
#include	"LiDIA/lattices/lattice_kernel.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifndef FUNCTION
template class lll_kernel_op< E_TYPE, A_TYPE, VECTOR< E_TYPE, A_TYPE >,
	MODULES< E_TYPE, A_TYPE, VAR > >;
#else
template class lll_kernel_fu< E_TYPE, A_TYPE, VECTOR< E_TYPE, A_TYPE >,
	MODULES< E_TYPE, A_TYPE, VAR > >;
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
