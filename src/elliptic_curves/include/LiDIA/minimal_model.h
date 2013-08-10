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
//	$Id: minimal_model.h,v 2.3 2002/06/24 09:41:45 lidiaadm Exp $
//
//	Author	: Nigel Smart (NS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MINIMAL_MODEL_H_GUARD_
#define LIDIA_MINIMAL_MODEL_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// The following function computes the minimal model (over Z) of an
// elliptic_curve < bigrational >, returning also the isomorphism from
// the original to the minimal model so that points may be mapped
// between them.  For sample usage, see minimal_model_appl.c

void
minimal_model(elliptic_curve< bigint > & Emin, elliptic_curve< bigrational > & ER,
	      curve_isomorphism< bigrational, bigint > & iso);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MINIMAL_MODEL_H_GUARD_
