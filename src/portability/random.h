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
//	$Id: random.h,v 2.2 2002/06/24 09:42:10 lidiaadm Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_RANDOM_H_GUARD_
#define LIDIA_RANDOM_H_GUARD_


#include	<stddef.h>



#ifdef __cplusplus
extern "C" {
#endif


void	srandom (unsigned x);
char*	initstate (unsigned seed, char *arg_state, size_t n);
char*	setstate (const char *arg_state);
long	random (void);


#ifdef __cplusplus
}
#endif



#endif	// LIDIA_RANDOM_H_GUARD_
