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
//	$Id: bigint_def.h,v 2.4 2002/06/24 09:41:55 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGINT_DEF_H_GUARD_
#define LIDIA_BIGINT_DEF_H_GUARD_



#include	<iint.h>
#include	<idigit.h>
#include	<imem.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



typedef Integer		bigint_rep_t;
typedef unsigned long	base_digit;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_BIGINT_DEF_H_GUARD_
