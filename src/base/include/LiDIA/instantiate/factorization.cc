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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


//
//  usage:
//
//	include all headers that are necessary to describe the type TYPE
//
//	define the type TYPE
//
//	include this file
//


#if defined (BASE_FACTOR)
# include	"LiDIA/base/base_factor.h"
# include	"LiDIA/base/base_factor.cc"
#endif

#if defined (SINGLE_FACTOR)
# include	"LiDIA/single_factor.h"
# include	"LiDIA/single_factor.cc"
#endif

#if defined (FACTORIZATION)
# include	"LiDIA/factorization.h"
# include	"LiDIA/factorization.cc"
# include	"LiDIA/base/ppair.h"
# include	"LiDIA/base/ppair.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if defined (BASE_FACTOR) || defined (SINGLE_FACTOR)

template class base_factor< TYPE >;

#endif


#if defined (SINGLE_FACTOR)

template class single_factor< TYPE >;

#endif


#if defined (FACTORIZATION)

template class factorization< TYPE >;

template class ppair< single_factor< TYPE >, lidia_size_t >;

#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
