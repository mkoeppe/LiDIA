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


#ifndef LIDIA_BASE_FACTOR_CC_GUARD_
#define LIDIA_BASE_FACTOR_CC_GUARD_


#ifndef LIDIA_BASE_FACTOR_H_GUARD_
# include	"LiDIA/base/base_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//****************************************************************
//			class base_factor< T >
//****************************************************************



template< class T >
base_factor< T > & base_factor< T >::operator = (const base_factor< T > &)
{
	lidia_error_handler("base_factor< T >",
			    "operator = (base_factor< T > &)::not implemented");
	return *this;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_FACTOR_CC_GUARD_
