// -*- C++ -*-
//=============================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//-----------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Christoph Ludwig (CL)
//	Changes	: See CVS log
//
//=============================================================================

#ifndef LIDIA_PARAMDESCRIPTION_CC_GUARD_
#define LIDIA_PARAMDESCRIPTION_CC_GUARD_


// LiDIA headers
#include "LiDIA/osstream.h"
#include "LiDIA/param_desc.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  namespace details {

    /**
     * class param_desc_impl
     */
    template <typename T>
    param_desc_impl<T>::
    param_desc_impl(const std::string& paramName,
			 const typename param_desc_impl::
			 value_type& value,
			 const std::string& paramCond)
      : param_desc(paramName, paramCond), value_(value) {
      //nothing to do
    }
    
    template <typename T>
    param_desc_impl<T>* param_desc_impl<T>::clone() const {
      return new param_desc_impl<T>(*this);
    }
    
    template <typename T>
    const std::string param_desc_impl<T>::valueStr() const {
      osstream oss;
      oss << value_;
      return extractString(oss);
    }
    
    template <typename T>
    const T& param_desc_impl<T>::value() const {
      return value_;
    }
    
  } // namespace details

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

#endif

