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
//	Author	:Thomas Pfahler (TPf), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_UDIGIT_DEF_H_GUARD_
#define LIDIA_UDIGIT_DEF_H_GUARD_


// safety belt
#if !defined LIDIA_ARITHMETIC || LIDIA_ARITHMETIC != LIDIA_ARITH_CLN
# error No or wrong arithmetic kernel specified!
#endif


#include	<cln/types.h>
#include        <cln/integer.h>

#define mp_limb_t uintD
#define mp_ptr uintD*
#define mp_srcptr const uintD*
#define mp_size_t long



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



typedef uintD	udigit;

enum {
	UDIGIT_NBITS = sizeof(udigit) * CHAR_BIT
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_UDIGIT_DEF_H_GUARD_
