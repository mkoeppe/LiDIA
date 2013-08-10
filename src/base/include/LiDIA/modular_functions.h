// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//	File    : modular_functions.h
//	Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_MODULAR_FUNCTIONS_H_GUARD_
#define LIDIA_MODULAR_FUNCTIONS_H_GUARD_


#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include "LiDIA/bigcomplex.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


bigcomplex
dedekind_eta( const bigcomplex &, bool is_fundamental = false );

bigcomplex
weber_f( const bigcomplex &, bool is_fundamental = false );

bigcomplex
weber_f1( const bigcomplex &, bool is_fundamental = false );

bigcomplex
weber_f2( const bigcomplex &, bool is_fundamental = false );

bigcomplex
gamma_2( const bigcomplex &, bool is_fundamental = false );

bigcomplex
modular_j( const bigcomplex &, bool is_fundamental = false );


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif	// LIDIA_MODULAR_FUNCTIONS_H_GUARD_
