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
//	Author	: Bruno Haible (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGINT_DEF_H_GUARD_
#define LIDIA_BIGINT_DEF_H_GUARD_


#include	<cln/dfloat.h>
#include	<cln/random.h>
#include	<cln/malloc.h>
#include	<cln/abort.h>
#include	<cln/integer_io.h>
#include	<cln/number.h>
#include	<cln/integer.h>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



typedef cln::cl_I	bigint_rep_t;
typedef uintD		base_digit;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_BIGINT_DEF_H_GUARD_
