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

#ifndef LIDIA_PRECONDITIONERROR_CC_GUARD_
#define LIDIA_PRECONDITIONERROR_CC_GUARD_


// LiDIA headers
#include "LiDIA/osstream.h"
#include "LiDIA/precondition_error.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  /**
   * member templates of class precondition_error
   */
  template <typename T>
  void precondition_error::addParam(const std::string& paramName, 
				   const T& value,
				   const std::string& paramCond) {
    details::param_desc* ptr =
	  new details::param_desc_impl<T>(paramName, value, paramCond);
    params_.push_back(ptr);
  }
  
  template <typename T>
  const T& precondition_error::paramValue(size_type n) const {
    typedef details::param_desc_impl<T> ImplType;

    // cast down and access value
    const ImplType* paramImplPtr =
      dynamic_cast<const ImplType*>(access(n, "paramValue<T>"));
    if(paramImplPtr == 0) { // cast failed
      osstream oss;
      oss << "precondition_error::paramValue<T>(" << n << ")";
      error_handler(cast_error(extractString(oss), "dynamic_cast failed"));
    }
    return paramImplPtr->value();
  }
  

template< class T1 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  error_handler(except);
}



template< class T1, class T2 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  error_handler(except);
}



template< class T1, class T2, class T3 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4, class T5 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4, class T5, class T6 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4,
          class T5, class T6, class T7 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4,
          class T5, class T6, class T7, class T8 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const T8& para8, const char *name_8, const char *cond8,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  except.addParam(name_8, para8, cond8);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4,
          class T5, class T6, class T7, class T8, class T9 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const T8& para8, const char *name_8, const char *cond8,
		    const T9& para9, const char *name_9, const char *cond9,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  except.addParam(name_8, para8, cond8);
  except.addParam(name_9, para9, cond9);
  error_handler(except);
}



template< class T1, class T2, class T3, class T4,
          class T5, class T6, class T7, class T8, class T9, class T10 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const T8& para8, const char *name_8, const char *cond8,
		    const T9& para9, const char *name_9, const char *cond9,
		    const T10& para10, const char *name_10, const char *cond10,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  except.addParam(name_8, para8, cond8);
  except.addParam(name_9, para9, cond9);
  except.addParam(name_10, para10, cond10);
  error_handler(except);
}



template < class T1, class T2, class T3, class T4,
           class T5, class T6, class T7, class T8,
           class T9, class T10, class T11 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const T8& para8, const char *name_8, const char *cond8,
		    const T9& para9, const char *name_9, const char *cond9,
		    const T10& para10, const char *name_10, const char *cond10,
		    const T11& para11, const char *name_11, const char *cond11,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  except.addParam(name_8, para8, cond8);
  except.addParam(name_9, para9, cond9);
  except.addParam(name_10, para10, cond10);
  except.addParam(name_11, para11, cond11);
  error_handler(except);
}



template < class T1, class T2, class T3, class T4,
           class T5, class T6, class T7, class T8,
           class T9, class T10, class T11, class T12 >
void
precondition_error_handler (const T1& para1, const char *name_1, const char *cond1,
		    const T2& para2, const char *name_2, const char *cond2,
		    const T3& para3, const char *name_3, const char *cond3,
		    const T4& para4, const char *name_4, const char *cond4,
		    const T5& para5, const char *name_5, const char *cond5,
		    const T6& para6, const char *name_6, const char *cond6,
		    const T7& para7, const char *name_7, const char *cond7,
		    const T8& para8, const char *name_8, const char *cond8,
		    const T9& para9, const char *name_9, const char *cond9,
		    const T10& para10, const char *name_10, const char *cond10,
		    const T11& para11, const char *name_11, const char *cond11,
		    const T12& para12, const char *name_12, const char *cond12,
		    const char *proto, const char *file, const char *msg)
{
  precondition_error except(proto, file, msg);
  except.addParam(name_1, para1, cond1);
  except.addParam(name_2, para2, cond2);
  except.addParam(name_3, para3, cond3);
  except.addParam(name_4, para4, cond4);
  except.addParam(name_5, para5, cond5);
  except.addParam(name_6, para6, cond6);
  except.addParam(name_7, para7, cond7);
  except.addParam(name_8, para8, cond8);
  except.addParam(name_9, para9, cond9);
  except.addParam(name_10, para10, cond10);
  except.addParam(name_11, para11, cond11);
  except.addParam(name_12, para12, cond12);
  error_handler(except);
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif // LIDIA_PRECONDITIONERROR_CC_GUARD_

//  Local Variables:    --- for emacs
//  mode: c++           --- for emacs
//  tab-width: 4        --- for emacs
//  End:                --- for emacs
