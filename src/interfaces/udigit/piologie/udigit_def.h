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
//	$Id: udigit_def.h,v 2.4 2002/06/24 09:41:57 lidiaadm Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_UDIGIT_DEF_H_GUARD_
#define LIDIA_UDIGIT_DEF_H_GUARD_


// safety belt
#if !defined LIDIA_ARITHMETIC || LIDIA_ARITHMETIC != LIDIA_ARITH_PIOLOGIE
# error No or wrong arithmetic kernel specified!
#endif


#include	<digit.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



typedef Digit udigit;

enum {
	UDIGIT_NBITS = BETA
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_UDIGIT_DEF_H_GUARD_
