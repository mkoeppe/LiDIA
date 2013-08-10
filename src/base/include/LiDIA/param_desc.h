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

#ifndef LIDIA_PARAMDESCRIPTION_H_GUARD_
#define LIDIA_PARAMDESCRIPTION_H_GUARD_

// LiDIA's headers
#include "LiDIA/LiDIA.h"
#include "LiDIA/base_vector.h"

// ISO C++ headers
#include <iosfwd>
#include <string>



#if defined(LIDIA_NAMESPACE)
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  class precondition_error; // defined in PreconditionError.h

  namespace details {

    class param_desc {
    public:
      param_desc(const std::string& paramName,
		       const std::string& paramCond);
      virtual ~param_desc();
      virtual param_desc* clone() const = 0;
      
      virtual const std::string valueStr() const = 0;
      const std::string& nameStr() const;
      const std::string& condStr() const;
       
    private:
      std::string name_;
      std::string cond_;
    };
    
    
    template <typename T>
    class param_desc_impl : public param_desc {
    public:
      typedef T value_type;
      
      param_desc_impl(const std::string& paramName,
			   const value_type& value,
			   const std::string& paramCond);
      virtual param_desc_impl* clone() const;
      virtual const std::string valueStr() const;
      const T& value() const;
      
    private:
      value_type  value_;
    };
    
    // Note:
    // I intended to use a container of refcounted smart pointers
    // (like boost::shared_ptr<T>). But it was too complicated to
    // track down all templates that need to be instantiated.
    // Since param_container is exclusively used in precondition_error and its
    // elements are never accessible from the outside, I decided
    // to use "dumb" pointers instead.
    class param_container {
    public:
      typedef lidia_size_t size_type;
      typedef param_desc* elem_type;

      virtual ~param_container();

      void push_back(elem_type pdp);
      const elem_type operator[](size_type n) const;
      size_type size() const;

      friend void swap(param_container&, param_container&);
      friend class LiDIA::precondition_error;

    private:
      explicit param_container(size_type n);
      param_container(const param_container&);
      param_container& operator=(const param_container&);
      
      base_vector<elem_type> container_;
    };

  }  // namespace details

    
  void swap(details::param_desc*& ptr1,
	    details::param_desc*& ptr2);
  std::ostream& operator<<(std::ostream& os,
			   const details::param_desc& pd);
  std::istream& operator>>(std::istream& is, details::param_desc*& ptr);

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#endif
