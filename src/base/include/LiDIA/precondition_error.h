// -*- C++ -*-
//=============================================================================
//
//    This file is part of LiDIA --- a library for computational number theory
//
//    Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//    See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//-----------------------------------------------------------------------------
//
//    $Id$
//
//    Author	: Christoph Ludwig (CL)
//    Changes	: See CVS log
//
//=============================================================================

#ifndef LIDIA_PRECONDITIONERROR_H_GUARD_
#define LIDIA_PRECONDITIONERROR_H_GUARD_

// ISO C/C++ headers
#include <iosfwd>        // std::ostream
#include <string>        // std::string

// LiDIA headers
#include "LiDIA/error.h"
#include "LiDIA/param_desc.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  class precondition_error : public basic_error, std::logic_error {
  public:
    typedef details::param_container::size_type size_type;
    
    precondition_error(const std::string& proto,
		      const std::string& where,
		      const std::string& what_msg);
    precondition_error(const precondition_error& pce);
    precondition_error& operator=(const precondition_error& pce);
    virtual ~precondition_error() throw();
    virtual const char* what() const throw();
    
    virtual const std::string& offendingClass() const;
    
    template <typename T>
    void addParam(const std::string& paramName, const T& value,
		  const std::string& paramCond);
    
    size_type noParams() const;

    // indexOf() determines the index of the stored parameter
    // named paramName.
    // If indexOf(paramName) == noParams() then there is no
    // such parameter stored. 
    // Otherwise, paramNameStr(indexOf(paramName)) == paramName is
    // guaranteed to be true and not to throw.
    size_type indexOf(const std::string& paramName) const;

    // paramNameStr throws LiDIA::index_out_of_bounds_error
    const std::string& paramNameStr(size_type n) const;
    
    // paramValueStr throws LiDIA::index_out_of_bounds_error
    const std::string paramValueStr(size_type n) const;
    
    // paramConditionStr throws LiDIA::index_out_of_bounds_error
    const std::string& paramConditionStr(size_type n) const;
    
    // paramValue throws LiDIA::index_out_of_bounds_error, LiDIA::cast_error
    template <typename T>
    const T& paramValue(size_type n) const; 

  protected:
    virtual void traditional_error_handler_impl(std::string& classname,
						std::string& msg) const;
    
  private:
    // access throws LiDIA::index_out_of_bounds_error
    const details::param_desc* access(size_type n,
					    const std::string& caller) const;

    details::param_container params_;
    std::string prototype_;
    std::string location_;
    
  };

  template< class T1 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7, class T8 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const T8 & para8, const char *name_8, const char *cond8,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7, class T8, class T9 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const T8 & para8, const char *name_8, const char *cond8,
		     const T9 & para9, const char *name_9, const char *cond9,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7, class T8, class T9, class T10 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const T8 & para8, const char *name_8, const char *cond8,
		     const T9 & para9, const char *name_9, const char *cond9,
		     const T10 & para10, const char *name_10,
		     const char *cond10,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7, class T8, class T9, class T10, class T11 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const T8 & para8, const char *name_8, const char *cond8,
		     const T9 & para9, const char *name_9, const char *cond9,
		     const T10 & para10, const char *name_10,
		     const char *cond10,
		     const T11 & para11, const char *name_11,
		     const char *cond11,
		     const char *proto, const char *file, const char *msg);
  
  template< class T1, class T2, class T3, class T4, class T5, class T6,
    class T7, class T8, class T9, class T10, class T11, class T12 >
  void
  precondition_error_handler(const T1 & para1, const char *name_1, const char *cond1,
		     const T2 & para2, const char *name_2, const char *cond2,
		     const T3 & para3, const char *name_3, const char *cond3,
		     const T4 & para4, const char *name_4, const char *cond4,
		     const T5 & para5, const char *name_5, const char *cond5,
		     const T6 & para6, const char *name_6, const char *cond6,
		     const T7 & para7, const char *name_7, const char *cond7,
		     const T8 & para8, const char *name_8, const char *cond8,
		     const T9 & para9, const char *name_9, const char *cond9,
		     const T10 & para10, const char *name_10,
		     const char *cond10,
		     const T11 & para11, const char *name_11,
		     const char *cond11,
		     const T12 & para12, const char *name_12,
		     const char *cond12,
		     const char *proto, const char *file, const char *msg);
  
  
#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/precondition_error.cc"
#endif

#endif  // LIDIA_PRECONDITIONERROR_H_GUARD_


//  Local Variables:    --- for emacs
//  mode: c++           --- for emacs
//  tab-width: 8        --- for emacs
//  End:                --- for emacs
