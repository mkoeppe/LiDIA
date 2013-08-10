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
//	$Id: interface_lib.h,v 2.5 2002/06/24 09:41:36 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_INTERFACE_LIB_H_GUARD_
#define LIDIA_INTERFACE_LIB_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef HEADBANGER

// file: single_precision.cc
sdigit       gcd (sdigit a, sdigit b);
sdigit       xgcd (sdigit & xv, sdigit & yv, sdigit a, sdigit b);
unsigned int power (unsigned int a, unsigned int b);
int          h_mul_mod(int, int, int);
int          mul_mod(int, int, int);
int          power_mod(int, int, int);
int          legendre(int, int);
int          invert(int, int);
int          jacobi(int, int);
int          ressol(int, int);

int integer_log (unsigned long x);

#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_INTERFACE_LIB_H_GUARD_
