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
//	$Id: bigint_def.h,v 2.3 2002/06/24 09:41:56 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGINT_DEF_H_GUARD_
#define LIDIA_BIGINT_DEF_H_GUARD_



#include	<integer.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



typedef Integer		bigint_rep_t;
typedef Digit		base_digit;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_BIGINT_DEF_H_GUARD_
