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


// LiDIA's headers
#include "LiDIA/precondition_error.h"
#include "LiDIA/osstream.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  
/**
 * class precondition_error
 */


  precondition_error::precondition_error(const std::string& proto,
				       const std::string& where,
				       const std::string& what_msg)
    : basic_error(what_msg.c_str()),
      std::logic_error(std::string("LiDIA: " + what_msg).c_str()),
	   params_(12), prototype_(proto), location_(where) {
    // Nothing to do
  }

  precondition_error::
  precondition_error(const precondition_error& pce)
	: basic_error(pce), std::logic_error(pce),
      params_(0), prototype_(pce.prototype_), location_(pce.location_) {
	// Create auxiliary container
	details::param_container tmpContainer(12);

	// Copy the elems from pce.params_ to tmpparam_container.
	// If an exception is thrown, the destructor of tmpparam_container
	// will free all acquired resources.
	size_type size = pce.params_.size();
	for(size_type index = 0; index < size; ++index) {
	  tmpContainer.push_back(pce.params_[index]->clone());
	}

	// swap tmpparam_container into params_. Won't throw.
	swap(params_, tmpContainer);
  }

  precondition_error& precondition_error::
  operator=(const precondition_error& pce) {
	precondition_error tmp(pce);

	// Hopefully, none of the remaining statements will throw.
	// (FIX ME: The code was written with the unbased assumption
	// that std::string::swap(std::string&) offers the nothrow-guarantee.
	// Unfortunately, it does not.)
	basic_error::operator=(pce);
	std::logic_error::operator=(pce);
	swap(params_, tmp.params_);
	// Note that some C++ headers have string swap as member function only
	prototype_.swap(tmp.prototype_);
	location_.swap(tmp.location_);
	return *this;
  }

  precondition_error::~precondition_error() throw() {
	// nothing to do
  }

  const char* precondition_error::what() const throw() {
    // this method needs only to disambiguitify which of the
    // inherited versions of what() is to be called.
    return basic_error::what();
  }

  const std::string& precondition_error::offendingClass() const {
    return location_;
  }

  precondition_error::size_type precondition_error::noParams() const {
    return params_.size();
  }

  precondition_error::size_type
  precondition_error::indexOf(const std::string& paramName) const {
    size_type size = params_.size();
    for(size_type index = 0; index < size; ++index) {
      if(paramNameStr(index) == paramName) {
	return index;
      }
    }
    return size;
  }

  const std::string&
  precondition_error::paramNameStr(precondition_error::size_type n) const {
    return access(n, "paramNameStr")->nameStr();
  }

  const std::string
  precondition_error::paramValueStr(precondition_error::size_type n) const {
    return access(n, "paramValueStr")->valueStr();
  }

  const std::string&
  precondition_error::paramConditionStr(precondition_error::size_type n) const {
    return access(n, "paramConditionStr")->condStr();
  }

  void precondition_error::
  traditional_error_handler_impl(std::string& classname,
				 std::string& msg) const {
    std::cerr << "FUNCTION:\n";
    std::cerr << prototype_ << "\n\n";
    std::cerr << "PARAMETERS:\n";

    size_type noElems = params_.size();
    for(size_type index = 0; index < noElems; ++index) {
      std::cerr << "actual parameter: " << *params_[index] << '\n';
    }
    std::cerr << '\n';
    std::cerr << "MESSAGE:\n";

    basic_error::traditional_error_handler_impl(classname, msg);
  }
  
  const details::param_desc*
  precondition_error::access(precondition_error::size_type n,
			    const std::string& caller) const {
    // We do not rely on the bounds check in basic_vector.member()
    // since we want to provide a helpfull error message.
    if(params_.size() <= n) {
      osstream oss;
      oss << "precondition_error::" + caller + '(' << n << ')';
      error_handler(index_out_of_bounds_error(extractString(oss), n));
    }

    return params_[n];
  }

    
#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

//  Local Variables:    --- for emacs
//  mode: c++           --- for emacs
//  tab-width: 4        --- for emacs
//  End:                --- for emacs
