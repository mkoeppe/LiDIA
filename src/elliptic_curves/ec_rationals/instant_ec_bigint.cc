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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#define HAS_NO_GET_FIELD

#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"

#include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.cc"
#include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//**** instantiation of classes *****************************

template class base_elliptic_curve_rep< bigint >;

template class base_elliptic_curve< bigint >;




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
